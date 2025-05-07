import csv
import os
import subprocess
import tempfile


def test_cli_csv_roundtrip():
    with tempfile.TemporaryDirectory() as tmpdir:
        input_path = os.path.join(tmpdir, "input.csv")
        output_path = os.path.join(tmpdir, "output.csv")

        # Write a basic test input
        with open(input_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["id", "ground_truth", "prediction"])
            writer.writerow(
                ["rxn1", "[CH3:1][OH:2]>>[CH3:1][O-:2]", "[CH3:1][OH:2]>>[CH3:1][O-:2]"]
            )
            writer.writerow(
                ["rxn2", "[CH3:2][OH:1]>>[CH3:1][O-:2]", "[CH3:9][OH:8]>>[CH3:9][O-:8]"]
            )

        # Run CLI
        result = subprocess.run(
            [
                "python",
                "-m",
                "atommap_eval.cli",
                input_path,
                "-f",
                "csv",
                "-p",
                "2",
                "-o",
                output_path,
            ],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 0
        assert "2 reactions" in result.stdout
        assert os.path.exists(output_path)

        # Check output content
        with open(output_path, newline="") as f:
            rows = list(csv.DictReader(f))
            assert rows[0]["id"] == "rxn1"
            assert rows[0]["equivalent"] == "True"
            assert rows[1]["id"] == "rxn2"
            assert rows[1]["equivalent"] == "False"
