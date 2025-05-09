from concurrent.futures import ProcessPoolExecutor
from typing import List, Optional, Tuple, Union
from functools import partial

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


def _unpack_pair(pair: Union[Tuple[str, str], ReactionPair]) -> Tuple[str, str]:
    """
    Extracts (gt, pred) from a tuple or a ReactionPair.
    """
    if isinstance(pair, tuple):
        return pair
    elif isinstance(pair, ReactionPair):
        return pair.ground_truth, pair.prediction
    else:
        raise TypeError(f"Unsupported pair type: {type(pair)}")


def _evaluate_pair(pair, canonicalize: bool = True) -> bool:
    from .evaluator import are_atom_maps_equivalent
    from .parallel import _unpack_pair

    rxn1, rxn2 = _unpack_pair(pair)
    return are_atom_maps_equivalent(rxn1, rxn2, canonicalize=canonicalize)


def evaluate_pairs_sequentially(
    pairs: List[Union[Tuple[str, str], ReactionPair]],
    canonicalize: bool = False,
) -> List[bool]:
    """
    Evaluate atom-map equivalence for a list of reaction pairs
    (tuples or dataclass) sequentially.

    Args:
        pairs: List of (gt_smiles, pred_smiles) tuples or ReactionPair objects.

    Returns:
        List[bool]: True/False for each comparison.
    """
    return [are_atom_maps_equivalent(*_unpack_pair(pair)) for pair in pairs]


def evaluate_pairs_in_parallel(
    pairs: List[Union[Tuple[str, str], ReactionPair]],
    canonicalize: bool = False,
    num_workers: Optional[int] = None,
    chunksize: int = 1,
) -> List[bool]:
    """
    Evaluate atom-map equivalence for a list of reaction pairs (tuples or dataclass)
    in parallel.

    Args:
        pairs: List of (gt_smiles, pred_smiles) tuples or ReactionPair objects.
        num_workers: Number of parallel workers. Defaults to number of CPU cores.
        chunksize: How many items per worker task.

    Returns:
        List[bool]: True/False for each comparison.
    """
    bound_fn = partial(_evaluate_pair, canonicalize=canonicalize)

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        results = list(executor.map(bound_fn, pairs, chunksize=chunksize))
    return results
