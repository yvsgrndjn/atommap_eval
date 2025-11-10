# Functions taken from RXNMapper original in rxnmapper/smiles_utils.py available here: https://github.com/rxn4chemistry/rxnmapper
import logging
from typing import List, Tuple

from rdkit import Chem
from rxn.chemutils.conversion import smiles_to_mol
from rxn.chemutils.exceptions import InvalidSmiles
from rxn.chemutils.reaction_equation import (
    ReactionEquation,
)
from rxn.chemutils.reaction_smiles import parse_any_reaction_smiles

LOGGER = logging.getLogger(__name__)


def create_rxn_smiles_from_rxn_equation(rxn_equation) -> str:
    """
    Converts a parsed reaction equation object into a SMILES string.

    Args:
        rxn_equation: Parsed reaction equation with reactants/products.

    Returns:
        str: Reaction SMILES string.
    """
    reactants = ".".join(rxn_equation.reactants)
    products = ".".join(rxn_equation.products)
    return f"{reactants}>>{products}"


def canonicalize_and_atom_map(smi: str, return_equivalent_atoms=False):
    """
    Remove atom mapping, canonicalize and return mapping numbers
    in order of canonicalization.

    Args:
        smi: reaction SMILES str
        return_equivalent_atoms

    Returns:

    """
    mol = smiles_to_mol(smi, sanitize=False)
    for atom in mol.GetAtoms():  # type: ignore[call-arg]
        if atom.HasProp("molAtomMapNumber"):
            atom_map = atom.GetAtomMapNum()
            atom.SetProp("atom_map", str(atom_map))
            atom.ClearProp("molAtomMapNumber")
        else:
            atom.SetProp("atom_map", str(0))
    can_smi = Chem.MolToSmiles(mol)
    order = list(
        mol.GetPropsAsDict(includePrivate=True, includeComputed=True)[
            "_smilesAtomOutputOrder"
        ]
    )

    atom_maps_canonical = [mol.GetAtoms()[idx].GetProp("atom_map") for idx in order]  # type: ignore[call-arg]

    if not return_equivalent_atoms:
        return (can_smi, atom_maps_canonical)

    raise NotImplementedError


def process_reaction_with_product_maps_atoms(rxn, skip_if_not_in_precursors=False):
    """
    Remove atom-mapping, move reagents to reactants and canonicalize reaction.
    If fragment group information is given, keep the groups together using
    the character defined with fragment_bond.

    Args:
        rxn: Reaction SMILES
        skip_if_not_in_precursors: accept unmapped atoms in the product (default: False)

    Returns: joined_precursors>>joined_products reaction SMILES
    """
    reactants, reagents, products = rxn.split(">")
    try:
        precursors = [canonicalize_and_atom_map(r) for r in reactants.split(".")]
        if len(reagents) > 0:
            precursors += [canonicalize_and_atom_map(r) for r in reagents.split(".")]
        products = [canonicalize_and_atom_map(p) for p in products.split(".")]
    except InvalidSmiles:
        return ""
    sorted_precursors = sorted(precursors, key=lambda x: x[0])
    sorted_products = sorted(products, key=lambda x: x[0])
    joined_precursors = ".".join([p[0] for p in sorted_precursors])
    joined_products = ".".join([p[0] for p in sorted_products])
    precursors_atom_maps = [
        i for p in sorted_precursors for i in p[1]
    ]  # could contain duplicate entries
    product_atom_maps = [
        i for p in sorted_products for i in p[1]
    ]  # could contain duplicate entries

    joined_rxn = f"{joined_precursors}>>{joined_products}"

    products_maps = []
    warnings = []

    for p_map in product_atom_maps:
        if skip_if_not_in_precursors and p_map not in precursors_atom_maps:
            products_maps.append(-1)
        elif int(p_map) == 0:
            products_maps.append(-1)
        else:
            corresponding_precursors_atom = precursors_atom_maps.index(p_map)
            if (
                corresponding_precursors_atom in products_maps
            ):  # handle equivalent atoms
                found_alternative = False
                for atom_idx, precursor_map in enumerate(precursors_atom_maps):
                    if (precursor_map == p_map) and atom_idx not in products_maps:
                        products_maps.append(atom_idx)
                        found_alternative = True
                        break
                if not found_alternative:
                    msg = f"Two product atoms mapped to the same precursor atom: {rxn}"
                    warnings.append(msg)
                    products_maps.append(corresponding_precursors_atom)
            else:
                products_maps.append(corresponding_precursors_atom)
    for w in list(set(warnings)):
        LOGGER.warning(w)
        raise ValueError(w)
    return joined_rxn, products_maps


def generate_atom_mapped_reaction_atoms(
    rxn: str, product_atom_maps, expected_atom_maps=None, canonical: bool = False
) -> Tuple[ReactionEquation, List[int]]:
    """
    Generate atom-mapped reaction from unmapped reaction and
    product-2-reactant atoms mapping vector.
    Args:
        rxn: unmapped reaction, in the format that the transformer model relies on.
        product_atom_maps: product to reactant atom maps.
        expected_atom_maps: if given, return the differences.
        canonical: whether to canonicalize the resulting SMILES.

    Returns: Atom-mapped reaction

    """

    reactants, agents, products = parse_any_reaction_smiles(rxn)
    precursors_mols = [smiles_to_mol(pr, sanitize=False) for pr in reactants + agents]
    products_mols = [smiles_to_mol(prod, sanitize=False) for prod in products]

    precursors_atom_maps = []

    differing_maps = []

    product_mapping_dict = {}

    i = -1
    atom_mapped_precursors_list = []
    for precursor_mol in precursors_mols:
        for atom in precursor_mol.GetAtoms():  # type: ignore[call-arg]
            i += 1
            if i in product_atom_maps:
                # atom maps start at an index of 1
                corresponding_product_atom_map = product_atom_maps.index(i) + 1
                precursors_atom_maps.append(corresponding_product_atom_map)
                atom.SetProp("molAtomMapNumber", str(corresponding_product_atom_map))

                indices = [idx for idx, x in enumerate(product_atom_maps) if x == i]

                if len(indices) > 1:
                    for idx in indices[1:]:
                        product_mapping_dict[idx] = corresponding_product_atom_map

                if expected_atom_maps is not None:
                    if (
                        i not in expected_atom_maps
                        or corresponding_product_atom_map
                        != expected_atom_maps.index(i) + 1
                    ):
                        differing_maps.append(corresponding_product_atom_map)
        atom_mapped_precursors_list.append(
            Chem.MolToSmiles(precursor_mol, canonical=canonical)
        )

    i = -1
    atom_mapped_products_list = []
    for products_mol in products_mols:
        for atom in products_mol.GetAtoms():  # type: ignore[call-arg]
            i += 1
            atom_map = product_mapping_dict.get(i, i + 1)
            atom.SetProp("molAtomMapNumber", str(atom_map))
        atom_mapped_products_list.append(
            Chem.MolToSmiles(products_mol, canonical=canonical)
        )

    atom_mapped_rxn = ReactionEquation(
        atom_mapped_precursors_list, [], atom_mapped_products_list
    )

    return atom_mapped_rxn, differing_maps
