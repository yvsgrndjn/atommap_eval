# atommap-eval

**Evaluate the equivalence of atom-mapped reaction SMILES using graph-based isomorphism.**

---

## Overview

`atommap-eval` is a Python package for comparing two atom-mapped reactions and determining whether they are chemically equivalent, using graph theory and RDKit.

It supports:

- ✅ Canonicalization and standardization of SMILES
- ✅ Graph construction with atom-level attributes and mapping
- ✅ Graph isomorphism checks using `networkx`
- ✅ Evaluation over single or batch pairs (in parallel or sequentially)
- ✅ CLI support for evaluating CSV/JSON files
- ✅ Fast testing and linting with `pytest`, `ruff`, and `black`

---

## Installation

```bash
# Clone the repo and install in editable mode
git clone https://github.com/yvsgrndjn/atommap-eval.git
cd atommap-eval
pip install -e ".[dev]"
```

or with Conda:
```bash
conda create -n atommap-eval python=3.9
conda activate atommap-eval
pip install -e ".[dev]"
```

---

## Usage
### Python
```python
from atommap_eval.evaluator import are_atoms_equivalent

gt = ""
pred = ""
result = are_atom_maps_equivalent(gt, pred)
print(result) # True
```

### CLI
```bash
atommap-eval reactions.csv -f csv -p 4 -o results.csv
```

---

## Project structure
```bash
src/atommap_eval/
├── canonicalizer.py
├── cli.py
├── data_models.py
├── evaluator.py
├── input_io.py
├── parallel.py
├── rxn_graph.py
├── rxnmapper_utils.py
tests/
```

---

## Development
### Run tests:
```bash
make test
```

### Format code:
```bash
make format
```

### Lint:
```bash
make lint
```

---

## Test examples
Unit tests are located under `test/` and cover evaluator logic, CLI execution, and multiprocessing correctness.

## License
MIT License © 2025 Yves Grandjean