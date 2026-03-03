# atommap_eval

**Evaluate the equivalence of atom-mapped reaction SMILES using graph-based isomorphism.**

## Overview

`atommap_eval` is a Python package for comparing two atom-mapped reactions and determining whether they are chemically equivalent, using their graph (`networkx`) representation and RDKit.

How it works:
- Optional preprocessing: "Canonicalization" and standardization of reaction SMILES to ensure all reactions are in the right format.
- Reactions graphs construction with atom-level / bond-level attributes and mapping
- Graph isomorphism checks using `networkx.is_isomorphic()`

It allows consistent evaluation of atom-mapping equivalence (e.g. against a ground truth atom-mapped reaction) by taking into account equivalence (i.e. are both atom-mapped reactions describing the same reaction)
of some atoms (i.e. all `CH3` in `t-Bu` are equivalent, any shuffling of atom-map indices should not impact correctness of the mapping)

*Warning:* tautomeric mappings are not considered equivalent even though from a chemist's perspective they are. Because template extraction
of the underlying reactivity would yield different results. Flags for tautomers will however be implemented in further implementations to better deal with this specific case.

By default, if the isomorphism takes more than 10 seconds, it is interrupted and returns None with status "timeout".

**Please read the expected atom-mapping format and preprocessing sections before using this package**

### Coming updates:
- cleaning preprocessing based on edge cases
- adding sanitization only as a basic preprocessing option
- test CLI for >1.0.0
- update all tests >1.0.0
- define clearly how evaluation needs to be considered and what are edge cases examples
- fix `evaluate_pairs_batched` for batch_size>1 that shuffles rows among the batch
---

## Installation

### Quick install for users (pip)
(version 1.4.0)
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
### Expected atom-mapped reactions format
The ideal format enforced by the preprocessing is: all atoms on the product side are atom-mapped, and each of the mapping numbers need to have an equivalent on the reactant side of the reaction. 
RXNMapper builds atom-maps by traveling through product atoms and finding their predicted equivalent on the other side of the reaction iteratively.
If your data might contain many cases that would be removed by the preprocessing, simply don't use it.
Sanitize all (ground_truth, prediction) atom-mapped reaction pairs and run evaluation on them. (Soon available as an argument in the evaluation), in the meantime use: `from atommap_eval.preprocess import sanitize_reaction_smiles`.


### Preprocessing
Preprocessing helps format (ground_truth, prediction) atom-mapped reactions for a fair evaluation, removing pairs that are considered unfit atom-mapped reaction format (coming from the ground_truth side). 
It is split in 2 parts:
- canonicalization + sanitization : sorts reaction SMILES and atom-mapping indices deterministically. Sanitizes reactions. Returns None if one of the steps fails (associated with flags A, B, C, S)
- Format analyis : raises specific flags (D) if preprocessing worked but the reaction format
will lead to a negative evaluation.
Preprocessing removes rows that:
- had any flag for the ground truth reaction
- predicted reactions that raise flag B (an atom on the product side is not on the reactant side, which arises from the ground truth and will be detected as such in the future)

**Different flags**
Hard stop flags are associated with a None output for the preprocessed reaction:
- A: two product atoms are mapped with same index (could be solved in the future)
- B: one of the atoms in the product has no counterpart on the reactants' side
- C: impossible to canonicalize reaction SMILES, usually because some `->` characters are found within the string
- S: reaction could not be sanitized
Warning flags indicate a reaction out of preprocessing that will fail during evaluation (faulty ground-truth format):
- D: one of the atoms on the product side is not atom-mapped  

To preprocess data, either use the simple wrapper if it matches your needs:
```python
import atommap_eval.preprocess as preprocess

preprocessed_df = preprocess.preprocess_dataset(df, path_to_save)
```


### Python
If you have few examples, use the following:
```python
# simple case
from atommap_eval.evaluator import are_atom_maps_equivalent

gt = "[C:1](=[O:2])[O-:3].[H+:4]>>[C:1](=[O:2])[OH:3]"
pred = "[H+:4].[C:1](=[O:2])[O-:3]>>[C:1](=[O:2])[OH:3]"
result = are_atom_maps_equivalent(gt, pred)
print(result) # True
```

However, if you have more reactions to evaluate, use:
```python
from atommap_eval.pair_evaluation import evaluate_pairs_batched

# `pairs` is either a list of tuples (rxn1, rxn2) or ReactionPair objects from atommap_eval.data_models
# for example if you store reactions in `your_df` under columns "ground_truth_rxn" and "predicted_rxn":
pairs = [
    ReactionPair(row.ground_truth_rxn, row.predicted_rxn)
    for _, row in your_df.iterrows()
]

results = evaluate_pairs_batched(pairs, batch_size=1) #for now, a bigger batch_size shuffles the results,
# results is a list of tuples (result: bool, status: str) where status can be "ok", "timeout", "error:{e}"
```
When using `evaluate_pairs_batched`, a timeout of 10 seconds per pair is enforced. 


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
├── pair_evaluation.py
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