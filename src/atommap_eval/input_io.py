import csv
import json
from typing import List

from .data_models import ReactionPair

"""
# Example CSV
id,ground_truth,prediction
rxn1,[CH3:1][OH:2]>>[CH3:1][O-:2],[Na+:3].[CH3:1][OH:2]>>[CH3:1][O-:2].[Na+:3]
rxn2,[C:1](=[O:2])[O-:3]>>[C:1](=[O:2])[OH:3],[C:1](=[O:2])[OH:3]>>[C:1](=[O:2])[OH:3]

# Example usage
from atommap_eval.input_io import load_pairs_from_csv
from atommap_eval.parallel import evaluate_pairs_in_parallel

pairs = load_pairs_from_csv("reactions.csv")
results = evaluate_pairs_in_parallel(pairs)

"""


def load_pairs_from_csv(
    path: str,
    gt_col: str = "ground_truth",
    pred_col: str = "prediction",
    id_col: str = "id",
) -> List[ReactionPair]:
    """
    Load reaction pairs from a CSV file and return a list of ReactionPair objects.

    Args:
        path (str): Path to CSV file.
        gt_col (str): Name of the column containing ground truth SMILES.
        pred_col (str): Name of the column containing predicted SMILES.
        id_col (str): Optional column for identifying reactions.

    Returns:
        List[ReactionPair]
    """
    pairs = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            gt = row[gt_col].strip()
            pred = row[pred_col].strip()
            rxn_id = row.get(id_col, None)
            pairs.append(ReactionPair(gt, pred, rxn_id))
    return pairs


def load_pairs_from_json(path: str) -> List[ReactionPair]:
    """
    Load reaction pairs from a JSON file containing a list of dicts with keys:
    "ground_truth", "prediction", and optionally "id".

    Args:
        path (str): Path to JSON file.

    Returns:
        List[ReactionPair]
    """
    with open(path, "r") as f:
        data = json.load(f)

    if not isinstance(data, list):
        raise ValueError("Expected a list of dictionaries in the JSON file.")

    pairs = []
    for entry in data:
        gt = entry["ground_truth"]
        pred = entry["prediction"]
        rxn_id = entry.get("id")
        pairs.append(ReactionPair(gt, pred, rxn_id))

    return pairs
