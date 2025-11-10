from dataclasses import dataclass
from typing import Optional, Dict

"""
# Example

from atommap_eval.data_models import ReactionPair

rp = ReactionPair(
    ground_truth="[CH3:1][OH:2]>>[CH3:1][O-:2]",
    prediction="[CH3:1][OH:2]>>[CH3:1][O-:2]",
    id="rxn_42"
)

print(rp.ground_truth)
print(rp.id)

"""


@dataclass(frozen=True)
class ReactionPair:
    """
    Represents a ground-truth / predicted reaction SMILES pair.

    Attributes:
        ground_truth (str): Ground-truth atom-mapped reaction SMILES.
        prediction (str): Predicted atom-mapped reaction SMILES.
        id (Optional[str]): Optional identifier for tracking the reaction
            (e.g., filename, reaction ID).
    """

    ground_truth: str
    prediction: str
    id: Optional[str] = None


@dataclass
class PreprocessResult:
    """
    Represents the preprocessing structure of a ground-truth / predicted reaction SMILES pair.

    Attributes:
        reference: preprocessed ground-truth reaction. None if unfit for evaluation.
        prediction: preprocessed predicted reaction SMILES pair. None if unfit for evaluation.
        flags_ref: notify either why the preprocessing failed or what formatting problems are shown by the ground-truth.
        flags_pred: notify either why the preprocessing failed or what formatting problems are shown by the prediction
    """
    reference: Optional[str]
    prediction: Optional[str]
    flags_ref: Dict[str, bool]
    flags_pred: Dict[str, bool]