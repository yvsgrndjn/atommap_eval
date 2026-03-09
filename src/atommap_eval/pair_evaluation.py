from concurrent.futures import ProcessPoolExecutor, TimeoutError, as_completed
from typing import List, Optional, Tuple, Union
import time

from wrapt_timeout_decorator import timeout

from .data_models import ReactionPair
from .evaluator import are_atom_maps_equivalent



def _unpack_pair(pair: Union[Tuple[str, str], "ReactionPair"]) -> Tuple[str, str]:
    """Extract (gt, pred) from a tuple or a ReactionPair."""
    if isinstance(pair, tuple):
        return pair
    elif hasattr(pair, "ground_truth") and hasattr(pair, "prediction"):
        return pair.ground_truth, pair.prediction
    else:
        raise TypeError(f"Unsupported pair type: {type(pair)}")


def _safe_equivalence_check(gt: str, pred: str, canonicalize: bool = False) -> bool:
    """Safely call are_atom_maps_equivalent in a subprocess."""
    return are_atom_maps_equivalent(gt_smi=gt, pred_smi=pred, canonicalize=canonicalize)


@timeout(dec_timeout=10.0, use_signals=False)
def timed_are_maps_equivalent(gt, pred, canonicalize):
    return are_atom_maps_equivalent(gt_smi=gt, pred_smi=pred, canonicalize=canonicalize)


def safe_timed_equivalence_check(gt: str, pred: str, canonicalize: bool=False):
    try:
        result = timed_are_maps_equivalent(gt, pred, canonicalize)
        return result, "ok"
    except TimeoutError:
        return None, "timeout"
    except Exception as e:
        return None, f"error:{type(e).__name__}"
    

def evaluate_pairs_batched(
    pairs: List[Union[Tuple[str, str], "ReactionPair"]],
    canonicalize: bool = False,
    num_workers: Optional[int] = None,
    batch_per_worker: int = 32,
) -> List[Tuple[Optional[bool], str]]:
    """
    Evaluate atom-map equivalence for many (gt, pred) pairs using batches and parallelism.

    Args:
        pairs: list of (gt_smiles, pred_smiles) tuples or ReactionPair objects
        canonicalize: whether to canonicalize before comparison
        num_workers: number of parallel workers
        batch_per_worker: target number of in-flight jobs per worker

    Returns:
        List of (equivalent, status) tuples.
    """
    results_all: List[Optional[Tuple[Optional[bool], str]]] = [None] * len(pairs)
    total_start = time.perf_counter()

    if not num_workers:
        num_workers = 1
    in_flight = batch_per_worker * num_workers

    print(
        f"[START] Evaluation of {len(pairs)} pairs with {num_workers} workers...",
        flush=True,
    )

    def submit_one(executor, idx, pair):
        try:
            gt, pred = _unpack_pair(pair)

            if (
                not isinstance(gt, str)
                or not isinstance(pred, str)
                or not gt.strip()
                or not pred.strip()
            ):
                results_all[idx] = (None, "invalid_input")
                return None

            return executor.submit(
                safe_timed_equivalence_check,
                gt,
                pred,
                canonicalize,
            )

        except Exception as e:
            results_all[idx] = (None, f"error:submit:{type(e).__name__}")
            return None

    def fill_pipeline(executor, iterator, future_to_idx, target_size):
        """
        Pull from iterator until future_to_idx reaches target_size
        or the iterator is exhausted.

        Invalid / bad inputs are consumed and written directly into results_all,
        but they do not contribute to future_to_idx size.
        """
        while len(future_to_idx) < target_size:
            try:
                idx, pair = next(iterator)
            except StopIteration:
                return

            fut = submit_one(executor, idx, pair)
            if fut is not None:
                future_to_idx[fut] = idx

    with ProcessPoolExecutor(max_workers=num_workers) as ex:
        it = iter(enumerate(pairs))
        future_to_idx = {}

        # initial fill
        fill_pipeline(ex, it, future_to_idx, in_flight)

        # drain one completed future at a time, then refill
        while future_to_idx:
            fut = next(as_completed(future_to_idx))
            idx = future_to_idx.pop(fut)
            try:
                result = fut.result()
                if result is None:
                    results_all[idx] = (None, "error:no_return")
                elif not isinstance(result, tuple) or len(result) != 2:
                    results_all[idx] = (None, f"error:bad_return:{type(result).__name__}")
                else:
                    results_all[idx] = result
            except Exception as e:
                results_all[idx] = (None, f"error:{type(e).__name__}")

            fill_pipeline(ex, it, future_to_idx, len(future_to_idx) + 1)

    bad = [i for i, x in enumerate(results_all) if x is None]
    if bad:
        raise RuntimeError(
            f"Unfilled result slots remain: first={bad[:10]}, count={len(bad)}"
        )

    total_elapsed = time.perf_counter() - total_start
    print(
        f"[DONE] Processed {len(pairs)} pairs in {total_elapsed:.2f}s using {num_workers} workers."
    )

    return results_all