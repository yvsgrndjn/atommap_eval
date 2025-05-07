from dataclasses import dataclass
from typing import Optional

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
