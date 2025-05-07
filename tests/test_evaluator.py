import pytest

from atommap_eval.evaluator import are_atom_maps_equivalent


@pytest.mark.parametrize(
    "gt, pred",
    [
        # Simple identical mappings
        ("[CH3:1][OH:2]>>[CH3:1][O-:2]", "[CH3:1][OH:2]>>[CH3:1][O-:2]"),
        # Order of reactants/products doesn't matter
        (
            "[CH3:1][OH:2].[Na+:3]>>[CH3:1][O-:2].[Na+:3]",
            "[Na+:3].[CH3:1][OH:2]>>[Na+:3].[CH3:1][O-:2]",
        ),
        # More complex molecule — same mapping
        (
            "[C:1](=[O:2])[O-:3].[H+:4]>>[C:1](=[O:2])[OH:3]",
            "[H+:4].[C:1](=[O:2])[O-:3]>>[C:1](=[O:2])[OH:3]",
        ),
    ],
)
def test_equivalent_mappings(gt, pred):
    assert are_atom_maps_equivalent(gt, pred) is True


@pytest.mark.parametrize(
    "gt, pred",
    [
        # Indices mixup
        ("[CH3:2][OH:1]>>[CH3:1][O-:2]", "[CH3:9][OH:8]>>[CH3:9][O-:8]"),
        # Mapped vs unmapped
        ("[CH3:1][OH:2]>>[CH3:1][O-:2]", "CO>>C[O-]"),
        # Different structures (non-equivalent chemistry)
        ("[CH3:1][OH:2]>>[CH3:1][O-:2]", "[CH3:1][NH2:2]>>[CH3:1][NH:2]"),
    ],
)
def test_non_equivalent_mappings(gt, pred):
    assert are_atom_maps_equivalent(gt, pred) is False


def test_invalid_smiles():
    with pytest.raises(RuntimeError):
        are_atom_maps_equivalent(
            "INVALID>>[CH3:1][OH:2]", "[CH3:1][OH:2]>>[CH3:1][O-:2]"
        )
