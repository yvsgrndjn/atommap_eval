# atommap_eval

**Evaluate the equivalence of atom-mapped reaction SMILES using graph-based isomorphism.**

## Overview

`atommap_eval` is a Python package for comparing two atom-mapped reactions and determining whether they are chemically equivalent, using their graph (`networkx`) representation and RDKit.

How it works:
- Optional preprocessing: "Canonicalization" and standardization of reaction SMILES to ensure all reactions are in the right format.
- Reactions graphs construction with atom-level / bond-level attributes and mapping
- Graph isomorphism checks using `networkx.is_isomorphic()`

It allows consistent evaluation of atom-mapping validity (e.g. against a ground truth atom-mapped reaction) by taking into account equivalence
of some atoms (i.e. all `CH3` in `t-Bu` are equivalent, any shuffling of atom-map indices should not impact correctness of the mapping)

*Warning:* tautomeric mappings are not considered equivalent even though from a chemist's perspective they are. Because template extraction
of the underlying reactivity would yield different results. Flags for tautomers will however be implemented in further implementations to better deal with this specific case.

---

## Installation

### Quick install for users (pip)
```bash
pip install atommap-eval
```

### For developers
```bash
# Clone the repo and install in editable mode
git clone https://github.com/yvsgrndjn/atommap_eval.git
cd atommap_eval
pip install -e ".[dev]"
```

or in case you want to create a new environment with Conda:
```bash
conda create -n atommap_eval python=3.9 -c conda-forge rdkit
conda activate atommap_eval
pip install -e ".[dev]"
```

---

## Usage
### Preprocessing
To preprocess data, either use the simple wrapper if it matches your needs:
```python
import atommap_eval.preprocess as preprocess

preprocess_df = preprocess.preprocess_dataset(df, path_to_save)
```

### Python
If you do not have too many examples, you can use the following:
```python
# simple case
from atommap_eval.evaluator import are_atom_maps_equivalent

gt = "[C:1](=[O:2])[O-:3].[H+:4]>>[C:1](=[O:2])[OH:3]"
pred = "[H+:4].[C:1](=[O:2])[O-:3]>>[C:1](=[O:2])[OH:3]"
result = are_atom_maps_equivalent(gt, pred)
print(result) # True
```

However, if you have more reactions to evaluate and what to compute them in parallel, use:
```python
from atommap_eval.parallel import evaluate_pairs_in_parallel

# pairs is either a list of tuples (rxn1, rxn2), or a ReactionPair object from atommap_eval.data_models
results = evaluate_pairs_in_parallel(pairs)
# results is a list of booleans
```

### CLI
```bash
atommap_eval reactions.csv -f csv -p 4 -o results.csv
```

---

## Project structure
```bash
src/atommap_eval/
├── preprocess.py
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