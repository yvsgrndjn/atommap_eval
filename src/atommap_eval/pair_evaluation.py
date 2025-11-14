from concurrent.futures import ProcessPoolExecutor, TimeoutError, as_completed
from typing import List, Optional, Tuple, Union
from functools import partial
import time
import multiprocessing as mp

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
    batch_size: int = 1000,
) -> List[Tuple[Optional[bool], str]]:
    """
    Evaluate atom-map equivalence for many (gt, pred) pairs using batches, timeouts, and parallelism.

    Args:
        pairs: list of (gt_smiles, pred_smiles) tuples or ReactionPair objects
        canonicalize: whether to canonicalize before comparison
        num_workers: number of parallel workers (defaults to os.cpu_count())
        batch_size: number of reactions to process per batch
        timeout: per-reaction timeout in seconds

    Returns:
        List[Tuple[Optional[bool], str]]:
            Each element is (equivalent, status)
            where status is in {"ok", "timeout", "invalid_input", "error:<type>"}
    """
    results_all: List[Tuple[Optional[bool], str]] = []
    total_start = time.perf_counter()

    print(f"Starting evaluation of {len(pairs)} pairs with {num_workers or 'default'} workers...", flush=True)

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        for start in range(0, len(pairs), batch_size):
            end = min(start + batch_size, len(pairs))
            print(f"Submitting batch {start}-{end} / {len(pairs)}...", flush=True)

            batch = pairs[start:end]
            futures = []

            for pair in batch:
                gt, pred = _unpack_pair(pair)
                if not isinstance(gt, str) or not isinstance(pred, str) or not gt.strip() or not pred.strip():
                    results_all.append((None, "invalid_input"))
                    continue
                futures.append(executor.submit(safe_timed_equivalence_check, gt, pred, canonicalize))

            # collect results
            for future in as_completed(futures):
                try:
                    equivalent, status = future.result()
                except Exception as e:
                    equivalent, status = None, f"error:{type(e).__name__}"
                results_all.append((equivalent, status))

    total_elapsed = time.perf_counter() - total_start
    print(f"[DONE] Processed {len(pairs)} pairs in {total_elapsed:.2f}s using {num_workers or 'default'} workers.")
    return results_all
