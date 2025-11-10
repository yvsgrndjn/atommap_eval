import re
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Dict, Iterable, List, Optional, Union

from atommap_eval.data_models import ReactionPair, PreprocessResult
from .rxnmapper_utils import (
    create_rxn_smiles_from_rxn_equation,
    generate_atom_mapped_reaction_atoms,
    process_reaction_with_product_maps_atoms,
)

# Level 0 - preprocessing
def canonicalize_rxn_smiles(rxn_smi: str) -> Optional[str]:
    """
    Canonicalizes a reaction SMILES string, reordering molecules
    and atom-mapping numbers deterministically.

    Args:
        rxn_smi (str): Atom-mapped reaction SMILES string.

    Returns:
        Optional[str]: Canonicalized atom-mapped reaction SMILES,
        or None if it fails.
    """
    try:
        rxn_smi_can_dict = process_reaction_with_product_maps_atoms(rxn_smi)

        if rxn_smi_can_dict is None or not isinstance(rxn_smi_can_dict, tuple):
            raise ValueError(
                f"Invalid output from "
                f"process_reaction_with_product_maps_atoms: {rxn_smi_can_dict}"
            )

        rxn_can, _ = generate_atom_mapped_reaction_atoms(
            *rxn_smi_can_dict, canonical=True
        )
        can_smi = create_rxn_smiles_from_rxn_equation(rxn_can)

        return can_smi
    except Exception as e:
        raise RuntimeError(
            f"Canonicalization failed for input SMILES: {rxn_smi}\nError: {e}"
        )


def sanitize_reaction_smiles(rxn_smi: str) -> Optional[str]:
    """
    Sanitize a reaction SMILES string.
    Returns an isomeric, canonical RDKit-consistent SMILES with proper aromaticity.

    Args:
        rxn_smi (str): Atom-mapped reaction SMILES string.
    
    Returns:
        Optional[str]: Sanitized isomeric, canonical SMILES,
        None if it fails.
    """
    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smi, useSmiles=True)
        if rxn is None:
            raise ValueError(f"Invalid reaction SMILES: {rxn_smi}")
    except Exception:
        return None

    def sanitize_mols(mols):
        sanitized = []
        for mol in mols:
            try:
                if mol is None:
                    continue
                mol = Chem.Mol(mol)
                Chem.SanitizeMol(mol)
                Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
                Chem.SetAromaticity(mol)
                sanitized.append(mol)
            except Exception:
                return None
        return sanitized

    reactants = sanitize_mols(rxn.GetReactants())
    products  = sanitize_mols(rxn.GetProducts())
    agents    = sanitize_mols(rxn.GetAgents())
    if reactants is None or products is None or agents is None:
        return None
    
    def to_smiles(mols):
        return '.'.join(sorted(
            Chem.MolToSmiles(m, canonical=True, isomericSmiles=True) for m in mols
        ))

    react_smiles   = to_smiles(reactants)
    agent_smiles   = to_smiles(agents)
    product_smiles = to_smiles(products)

    return f"{react_smiles}>{agent_smiles}>{product_smiles}"


# Level 0 - formatting diagnostics
def detect_atom_no_mapnum_in_preproc_prod(rxn_smi: str) -> bool:
    """
    Return True if the reaction contains at least one atom in the product side
    that does not have an atom-mapping number (after canonicalization, one atom-map number won't be found in the reactants side).
    Returns False otherwise.
    Returns None if the reaction SMILES cannot be parsed.
    """
    try:
        rxn = AllChem.ReactionFromSmarts(rxn_smi, useSmiles=True)
        if rxn is None:
            return None
    except Exception:
        return None
    set_prod_map = set()
    set_reac_map = set()
    for product in rxn.GetProducts():
        for atom in product.GetAtoms():
            amap = atom.GetAtomMapNum()
            if amap:
                set_prod_map.add(amap)
    for reactant in rxn.GetReactants():
        for atom in reactant.GetAtoms():
            amap = atom.GetAtomMapNum()
            if amap:
                set_reac_map.add(amap)
    
    missing = not set_prod_map.issubset(set_reac_map)
    return missing


# Level 1 - preprocessing
def preprocess_reaction_pair(ref_rxn: str, pred_rxn: str) -> PreprocessResult:
    """
    Run canonicalization + sanitization for a pair of reactions.
    Flags:
        A: two product atoms are mapped to the same precursor (i.e. two product atoms have the same mapping number)
        B: canonicalization failed because atom-mapping number not in reactants
        C: other canonicalization issue
        S: sanitization failed
    """
    base_flags = {"A": False, "B": False, "C": False, "S": False}

    def preprocess_one(rxn, flags):
        try: # canonicalization
            can_rxn = canonicalize_rxn_smiles(rxn)
        except Exception as e:
            msg = str(e)
            if "Two product atoms mapped to the same precursor atom" in msg:
                flags["A"] = True
            elif re.search(r"Error:\s*'(\d+)' is not in list", msg):
                flags["B"] = True
            else: # should be "Canonicalization failed for input SMILES:...."
                flags["C"] = True
            return None
        try: # sanitization
            sanitized_rxn = sanitize_reaction_smiles(can_rxn)
            if sanitized_rxn is None:
                flags["S"] = True
                return None
        except Exception:
            flags["S"] = True
            return None
        return sanitized_rxn
    
    flags_ref = base_flags.copy()
    flags_pred = base_flags.copy()

    ref_processed = preprocess_one(ref_rxn, flags_ref)
    pred_processed = preprocess_one(pred_rxn, flags_pred)
    
    return PreprocessResult(
        reference=ref_processed,
        prediction=pred_processed,
        flags_ref=flags_ref,
        flags_pred=flags_pred,
    )


# Level 1 - formatting diagnostics
def diagnostic_flags_reaction(rxn_smi: str) -> Dict[str, bool]:
    """
    Run validation checks on the canonicalized, sanitized reaction.
    Returns diagnostic flags that do not block evaluation
    """
    flags = {
        "D": False, # Unmapped atom in product
        #"E": False, # eventually future diagnostics to run
    }
    if detect_atom_no_mapnum_in_preproc_prod(rxn_smi):
        flags["D"] = True
    return flags


# Level 2 - preprocessing + diagnostics
def batch_preprocess_reaction_pairs(
    pairs: Iterable[Union[ReactionPair, tuple[str, str]]]
) -> List[PreprocessResult]:
    """
    Apply preprocessing and diagnostic flagging to a list of reaction pairs.

    Args:
        pairs: Iterable of either ReactionPair objects or (ground_truth, prediction) tuples.

    Returns:
        List of PreprocessResult objects (one per pair).
    """
    results = []

    for pair in pairs:
        if isinstance(pair, ReactionPair):
            ref_rxn = pair.ground_truth
            pred_rxn = pair.prediction
        else:
            ref_rxn, pred_rxn = pair  # tuple[str, str]

        # --- step 1: preprocessing ---
        result = preprocess_reaction_pair(ref_rxn=ref_rxn, pred_rxn=pred_rxn)

        # --- step 2: diagnostics ---
        if result.reference:
            diag_flag = diagnostic_flags_reaction(result.reference)
            result.flags_ref.update(diag_flag)
        if result.prediction:
            diag_flag = diagnostic_flags_reaction(result.prediction)
            result.flags_pred.update(diag_flag)

        results.append(result)

    return results


def convert_preproc_to_df_via_list(results: List[PreprocessResult]) -> pd.DataFrame:
    """
    Converts results of preprocessing (i.e. a list of PreprocessResult instances) into pd.DataFrame.

    Args:
        results: list of PreprocessResults instances, out from `batch_preprocess_reaction_pairs()`
    
    Returns:
        pd.DataFrame, same format as original file with preprocessed reactions and flags.
    """
    data = []
    for i, r in enumerate(results):
        ref_flags = [k for k,v in r.flags_ref.items() if v]
        pred_flags = [k for k,v in r.flags_pred.items() if v]

        row = {
            "pair_index": i,
            "preproc_ref": r.reference, 
            "preproc_pred": r.prediction,
            "flags_ref": ref_flags or None,
            "flags_pred": pred_flags or None,
        }
        data.append(row)
    return pd.DataFrame(data)


def convert_preproc_to_df(results: List[PreprocessResult]) -> pd.DataFrame:
    """
    Convert a list of PreprocessResult instances into a flat, analysis-ready DataFrame.

    Each flag is expanded into its own Boolean column (ref_A, ref_B, ... / pred_A, pred_B, ...),
    allowing for direct filtering and aggregation.

    Args:
        results: List of PreprocessResult objects (output of `batch_preprocess_reaction_pairs()`)

    Returns:
        pd.DataFrame with columns:
            - pair_index: integer index of the pair
            - preproc_ref: preprocessed reference reaction
            - preproc_pred: preprocessed predicted reaction
            - ref_<FLAG>: boolean column for each flag in reference
            - pred_<FLAG>: boolean column for each flag in prediction
    """
    rows = []
    all_flags = None
    for r in results:
        if r.reference is not None or r.prediction is not None:
            all_flags = sorted({*r.flags_ref.keys(), *r.flags_pred.keys()})
            break
    if all_flags is None:
        all_flags = []

    for i, r in enumerate(results):
        row = {
            "pair_index": i,
            "preproc_ref": r.reference,
            "preproc_pred": r.prediction,
        }

        # one column per flag: True/False for ref and pred sides
        for flag in all_flags:
            row[f"ref_{flag}"] = bool(r.flags_ref.get(flag, False))
            row[f"pred_{flag}"] = bool(r.flags_pred.get(flag, False))

        rows.append(row)

    return pd.DataFrame(rows)



def preprocess_dataset(df: pd.DataFrame, path_to_save: Optional[str] = None):
    """
    Simple placeholder wrapper of the whole preprocessing approach for pd.DataFrame inputs.

    Args:
        df: pd.DataFrame that contains colunms `ground_truth_rxn`, `predicted_rxn` 
        and `row_idx` to preprocess.
        path_to_save (optional): path to file where the preprocessed dataframe as CSV should be saved

    Returns:
        pd.DataFrame: DataFrame with same format as original, with added columns describing what happened.
        Data is ready to be sorted and evaluated.
    """
    # check that columns "ground_truth_rxn", "predicted_rxn", "row_idx" exist
    if not set(["ground_truth_rxn", "predicted_rxn"]).issubset(df.columns):
        raise ValueError(f"df does not contain 'ground_truth_rxn' and 'predicted_rxn' columns, only: {df.columns}")
    
    pairs = [
        ReactionPair(row.ground_truth_rxn, row.predicted_rxn, id=row.get("row_idx"))
        for _, row in df.iterrows()
    ]
    results = batch_preprocess_reaction_pairs(pairs)
    df_preproc = convert_preproc_to_df(results)
    merged = df.merge(df_preproc, how="left", left_index=True, right_on="pair_index").drop(columns=["pair_index"])
    if path_to_save:
        merged.to_csv(path_to_save, index=False)
    return merged
