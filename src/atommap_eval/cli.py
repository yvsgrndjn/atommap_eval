import argparse
import csv
from typing import List

from .data_models import ReactionPair
from .input_io import load_pairs_from_csv, load_pairs_from_json
from .pair_evaluation import evaluate_pairs_in_parallel

"""
# Example usages
## CSV
python -m atommap_eval.cli reactions.csv -f csv -p 4 -o results.csv
## JSON
python -m atommap_eval.cli reactions.json -f json

"""


def parse_args():
    parser = argparse.ArgumentParser(
        description="Evaluate equivalence of atom-mapped reaction SMILES."
    )
    parser.add_argument("input", help="Path to input file (CSV or JSON)")
    parser.add_argument(
        "-f",
        "--format",
        choices=["csv", "json"],
        required=True,
        help="Input file format",
    )
    parser.add_argument("-o", "--output", help="Path to save output results as CSV")
    parser.add_argument(
        "-p",
        "--processes",
        type=int,
        default=None,
        help="Number of processes (default: all CPUs)",
    )
    parser.add_argument(
        "--gt-col", default="ground_truth", help="CSV column for ground truth SMILES"
    )
    parser.add_argument(
        "--pred-col", default="prediction", help="CSV column for predicted SMILES"
    )
    parser.add_argument(
        "--id-col", default="id", help="CSV column for reaction ID (optional)"
    )
    return parser.parse_args()


def load_pairs(
    input_path: str, fmt: str, gt_col: str, pred_col: str, id_col: str
) -> List[ReactionPair]:
    if fmt == "csv":
        return load_pairs_from_csv(
            input_path, gt_col=gt_col, pred_col=pred_col, id_col=id_col
        )
    elif fmt == "json":
        return load_pairs_from_json(input_path)
    else:
        raise ValueError("Unsupported format")


def write_results(pairs: List[ReactionPair], results: List[bool], output_path: str):
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["id", "ground_truth", "prediction", "equivalent"])
        for pair, result in zip(pairs, results):
            writer.writerow([pair.id or "", pair.ground_truth, pair.prediction, result])


def main():
    args = parse_args()

    print(f"Loading input pairs from: {args.input}")
    pairs = load_pairs(args.input, args.format, args.gt_col, args.pred_col, args.id_col)

    print(
        f"Running equivalence evaluation on {len(pairs)} pairs using "
        f"{args.processes or 'all available'} processes..."
    )
    results = evaluate_pairs_in_parallel(pairs, num_workers=args.processes)

    num_equiv = sum(results)
    print(f"Finished! {num_equiv} / {len(pairs)} reactions are equivalent.")

    if args.output:
        print(f"Saving results to {args.output}")
        write_results(pairs, results, args.output)


if __name__ == "__main__":
    main()
