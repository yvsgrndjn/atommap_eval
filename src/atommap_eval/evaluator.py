import networkx as nx

from .rxn_graph import build_reaction_graph, define_node_edge_matches

# Define node and edge match only once at import
_NODE_MATCH, _EDGE_MATCH = define_node_edge_matches()


def are_atom_maps_equivalent(gt_smi: str, pred_smi: str) -> bool:
    """
    Compare two atom-mapped reaction SMILES to determine whether
    the atom mappings are equivalent.

    This comparison is based on:
    - Canonicalizing the input SMILES
    - Constructing reaction graphs
    - Evaluating graph isomorphism with matching atom/bond features

    Args:
        gt_smi (str): Ground truth atom-mapped reaction SMILES
        pred_smi (str): Predicted atom-mapped reaction SMILES

    Returns:
        bool: True if the atom mappings are equivalent, False otherwise
    """
    rxn_gt_G = build_reaction_graph(gt_smi)
    rxn_pred_G = build_reaction_graph(pred_smi)

    return nx.is_isomorphic(
        rxn_gt_G, rxn_pred_G, node_match=_NODE_MATCH, edge_match=_EDGE_MATCH
    )
