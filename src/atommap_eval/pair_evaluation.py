from concurrent.futures import ProcessPoolExecutor, TimeoutError, as_completed
from typing import List, Optional, Tuple, Union
from functools import partial
import time

from .data_models import ReactionPair
from .evaluator import are_atom_maps_equivalent

"""
evalute_pairs_sequentially : great for debugging or small jobs
evaluate_pairs_in_parallel : multi-core batch execution
    (use chunksize > 1 for large datasets to reduce overhead)

# Example usage
```
from atommap_eval.parallel import evaluate_pairs_in_parallel
from atommap_eval.data_models import ReactionPair

pairs = [
    ReactionPair(
        "[CH3:1][OH:2]>>[CH3:1][O-:2]",
        "[Na+:3].[CH3:1][OH:2]>>[CH3:1][O-:2].[Na+:3]", id="rxn_1"
    ),
    (
        "[C:1](=[O:2])[O-:3]>>[C:1](=[O:2])[OH:3]",
        "[C:1](=[O:2])[OH:3]>>[C:1](=[O:2])[OH:3]"
    ),
]

results = evaluate_pairs_in_parallel(pairs, num_workers=4)
print(results)  # [True, True]

```

"""

# --- helpers ---

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


# --- main function ---

def evaluate_pairs_batched(
    pairs: List[Union[Tuple[str, str], "ReactionPair"]],
    canonicalize: bool = False,
    num_workers: Optional[int] = None,
    batch_size: int = 1000,
    timeout: float = 1.0,
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
            where status ∈ {"ok", "timeout", "invalid_input", "error:<type>"}
    """
    results_all: List[Tuple[Optional[bool], str]] = []
    total_start = time.perf_counter()

    # iterate over batches
    for start in range(0, len(pairs), batch_size):
        batch = pairs[start:start + batch_size]

        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            futures = []

            for pair in batch:
                gt, pred = _unpack_pair(pair)

                # validate quickly
                if not isinstance(gt, str) or not isinstance(pred, str) or not gt.strip() or not pred.strip():
                    results_all.append((None, "invalid_input"))
                    continue

                futures.append(executor.submit(_safe_equivalence_check, gt, pred, canonicalize))

            # collect as completed
            for future in as_completed(futures):
                try:
                    equivalent = future.result(timeout=timeout)
                    status = "ok"
                except TimeoutError:
                    equivalent, status = None, "timeout"
                except Exception as e:
                    equivalent, status = None, f"error:{type(e).__name__}"

                results_all.append((equivalent, status))

    total_elapsed = time.perf_counter() - total_start
    print(f"Processed {len(pairs)} pairs in {total_elapsed:.2f}s using {num_workers or 'default'} workers.")
    return results_all
