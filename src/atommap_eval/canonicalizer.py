from typing import Optional

from .rxnmapper_utils import (
    create_rxn_smiles_from_rxn_equation,
    generate_atom_mapped_reaction_atoms,
    process_reaction_with_product_maps_atoms,
)


def canonicalize_rxn_smiles(rxn_smi: str) -> Optional[str]:
    """
    Canonicalizes a reaction SMILES string, reordering molecules
    and atom-maps deterministically.

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
