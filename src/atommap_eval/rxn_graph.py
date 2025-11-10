from collections import defaultdict
from typing import List

import networkx as nx
from networkx.algorithms.isomorphism import (
    categorical_edge_match,
    categorical_node_match,
)
from rdkit import Chem
from rdkit.Chem import AllChem

from .preprocess import canonicalize_rxn_smiles


def build_atom_graph(mol: Chem.Mol) -> nx.Graph:
    """
    Convert an RDKit molecule to a NetworkX graph with key atom and bond features.
    """
    G = nx.Graph()
    for atom in mol.GetAtoms():
        G.add_node(
            atom.GetIdx(),
            atomic_num=atom.GetAtomicNum(),
            atom_map=atom.GetAtomMapNum(),
            formal_charge=atom.GetFormalCharge(),
            degree=atom.GetDegree(),
        )
    for bond in mol.GetBonds():
        G.add_edge(
            bond.GetBeginAtomIdx(),
            bond.GetEndAtomIdx(),
            bond_type=str(bond.GetBondType()),
        )
    return G


def get_mapped_graphs(mols: List[Chem.Mol]) -> List[nx.Graph]:
    """
    Keep only molecule graphs with at least one atom-mapped atom.
    """
    graphs = []
    for mol in mols:
        g = build_atom_graph(mol)
        if any(g.nodes[n]["atom_map"] > 0 for n in g.nodes):
            graphs.append(g)
    return graphs


def union_labeled_graphs(
    reactants: List[nx.Graph], products: List[nx.Graph]
) -> nx.Graph:
    """
    Create a union graph with uniquely labeled nodes from reactants and products.
    """
    unionG = nx.Graph()

    def relabel_and_add(graphs, prefix):
        for i, g in enumerate(graphs):
            mapping = {n: f"{prefix}{i}-{n}" for n in g.nodes}
            renamed = nx.relabel_nodes(g, mapping)
            unionG.update(renamed)

    relabel_and_add(reactants, "R-")
    relabel_and_add(products, "P-")

    return unionG


def add_mapping_edges(unionG: nx.Graph) -> nx.Graph:
    """
    Add special edges between reactant and product atoms with the same atom_map.
    """
    atom_map_dict = nx.get_node_attributes(unionG, "atom_map")
    amap_to_nodes = defaultdict(list)
    for node, amap in atom_map_dict.items():
        if amap > 0:
            amap_to_nodes[amap].append(node)

    for amap, nodes in amap_to_nodes.items():
        rs = [n for n in nodes if n.startswith("R-")]
        ps = [n for n in nodes if n.startswith("P-")]
        for r in rs:
            for p in ps:
                unionG.add_edge(r, p, label="mapping")

    return unionG


def build_reaction_graph(rxn_smi: str, canonicalize: bool = True) -> nx.Graph:
    """
    Build a reaction graph from a reaction SMILES, optionally canonicalizing it first.

    Args:
        rxn_smi: Atom-mapped reaction SMILES
        canonicalize: Whether to canonicalize the input before processing,
        canonicalizing both SMILES and atom-map order

    Returns:
        NetworkX graph of the reaction
    """

    if not isinstance(rxn_smi, str):
        raise TypeError(f"Expected string, got {type(rxn_smi)}")

    if canonicalize:
        rxn_smi = canonicalize_rxn_smiles(rxn_smi)
    rxn = AllChem.ReactionFromSmarts(rxn_smi, useSmiles=True)

    reactant_graphs = get_mapped_graphs(rxn.GetReactants())
    product_graphs = get_mapped_graphs(rxn.GetProducts())

    union = union_labeled_graphs(reactant_graphs, product_graphs)
    return add_mapping_edges(union)


def define_node_edge_matches():
    """
    Define node and edge matchers used for graph isomorphism testing.
    """
    node_match = categorical_node_match(
        ["atomic_num", "formal_charge", "degree"], [None, None, None]
    )
    edge_match = categorical_edge_match(["bond_type", "label"], [None, None])
    return node_match, edge_match
