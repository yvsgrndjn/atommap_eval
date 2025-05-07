from atommap_eval.data_models import ReactionPair
from atommap_eval.parallel import (
    evaluate_pairs_in_parallel,
    evaluate_pairs_sequentially,
)


def test_parallel_equivalence_on_tuples():
    pairs = [
        ("[CH3:1][OH:2]>>[CH3:1][O-:2]", "[CH3:1][OH:2]>>[CH3:1][O-:2]"),
        ("[CH3:2][OH:1]>>[CH3:1][O-:2]", "[CH3:9][OH:8]>>[CH3:9][O-:8]"),
    ]
    expected = [True, False]
    results = evaluate_pairs_in_parallel(pairs)
    assert results == expected


def test_parallel_equivalence_on_reaction_pairs():
    pairs = [
        ReactionPair(
            "[CH3:1][OH:2]>>[CH3:1][O-:2]", "[CH3:1][OH:2]>>[CH3:1][O-:2]", id="rxn1"
        ),
        ReactionPair(
            "[CH3:2][OH:1]>>[CH3:1][O-:2]", "[CH3:9][OH:8]>>[CH3:9][O-:8]", id="rxn2"
        ),
    ]
    expected = [True, False]
    results = evaluate_pairs_sequentially(pairs)
    assert results == expected
