import pytest

from atommap_eval.evaluator import are_atom_maps_equivalent


@pytest.mark.parametrize(
    "gt, pred",
    [
        # Identical reactions from index 102 in Schneider validation set
        # ground truth twice
        (
            "CC(=O)O.CS(C)=O.Cl[c:1]1[c:2]([O:3][CH:4]([C:5]([F:6])([F:7])[F:8])[c:9]2[cH:10][cH:11][cH:12][n:13][cH:14]2)[n:15][c:16]2[cH:17][cH:18][cH:19][cH:20][c:21]2[n:22]1.O=C([O-])[O-].[K+].[K+].[NH2:23][S:24](=[O:25])(=[O:26])[C:27]([F:28])([F:29])[F:30]>>[c:1]1([NH:23][S:24](=[O:25])(=[O:26])[C:27]([F:28])([F:29])[F:30])[c:2]([O:3][CH:4]([C:5]([F:6])([F:7])[F:8])[c:9]2[cH:10][cH:11][cH:12][n:13][cH:14]2)[n:15][c:16]2[cH:17][cH:18][cH:19][cH:20][c:21]2[n:22]1",
            "CC(=O)O.CS(C)=O.Cl[c:1]1[c:2]([O:3][CH:4]([C:5]([F:6])([F:7])[F:8])[c:9]2[cH:10][cH:11][cH:12][n:13][cH:14]2)[n:15][c:16]2[cH:17][cH:18][cH:19][cH:20][c:21]2[n:22]1.O=C([O-])[O-].[K+].[K+].[NH2:23][S:24](=[O:25])(=[O:26])[C:27]([F:28])([F:29])[F:30]>>[c:1]1([NH:23][S:24](=[O:25])(=[O:26])[C:27]([F:28])([F:29])[F:30])[c:2]([O:3][CH:4]([C:5]([F:6])([F:7])[F:8])[c:9]2[cH:10][cH:11][cH:12][n:13][cH:14]2)[n:15][c:16]2[cH:17][cH:18][cH:19][cH:20][c:21]2[n:22]1",
        ),
        # Identical mappings from index 102 in Schneider validation set
        # ground truth and prediction
        (
            "CC(=O)O.CS(C)=O.Cl[c:1]1[c:2]([O:3][CH:4]([C:5]([F:6])([F:7])[F:8])[c:9]2[cH:10][cH:11][cH:12][n:13][cH:14]2)[n:15][c:16]2[cH:17][cH:18][cH:19][cH:20][c:21]2[n:22]1.O=C([O-])[O-].[K+].[K+].[NH2:23][S:24](=[O:25])(=[O:26])[C:27]([F:28])([F:29])[F:30]>>[c:1]1([NH:23][S:24](=[O:25])(=[O:26])[C:27]([F:28])([F:29])[F:30])[c:2]([O:3][CH:4]([C:5]([F:6])([F:7])[F:8])[c:9]2[cH:10][cH:11][cH:12][n:13][cH:14]2)[n:15][c:16]2[cH:17][cH:18][cH:19][cH:20][c:21]2[n:22]1",
            "CC(=O)O.CS(C)=O.Cl[c:5]1[n:6][c:7]2[cH:8][cH:9][cH:10][cH:11][c:12]2[n:13][c:14]1[O:15][CH:16]([c:17]1[cH:18][cH:19][cH:20][n:21][cH:22]1)[C:23]([F:24])([F:25])[F:26].[O:1]=[S:2](=[O:3])([NH2:4])[C:27]([F:28])([F:29])[F:30].O=C([O-])[O-].[K+].[K+]>>[O:1]=[S:2](=[O:3])([NH:4][c:5]1[n:6][c:7]2[cH:8][cH:9][cH:10][cH:11][c:12]2[n:13][c:14]1[O:15][CH:16]([c:17]1[cH:18][cH:19][cH:20][n:21][cH:22]1)[C:23]([F:24])([F:25])[F:26])[C:27]([F:28])([F:29])[F:30]",
        ),
        # Identical mappings from index 102 in Schneider validation set
        # with modified Fluorine maps, ground truth and prediction
        (
            "CC(=O)O.CS(C)=O.Cl[c:1]1[c:2]([O:3][CH:4]([C:5]([F:6])([F:7])[F:8])[c:9]2[cH:10][cH:11][cH:12][n:13][cH:14]2)[n:15][c:16]2[cH:17][cH:18][cH:19][cH:20][c:21]2[n:22]1.O=C([O-])[O-].[K+].[K+].[NH2:23][S:24](=[O:25])(=[O:26])[C:27]([F:28])([F:29])[F:30]>>[c:1]1([NH:23][S:24](=[O:25])(=[O:26])[C:27]([F:28])([F:29])[F:30])[c:2]([O:3][CH:4]([C:5]([F:6])([F:7])[F:8])[c:9]2[cH:10][cH:11][cH:12][n:13][cH:14]2)[n:15][c:16]2[cH:17][cH:18][cH:19][cH:20][c:21]2[n:22]1",
            "CC(=O)O.CS(C)=O.Cl[c:5]1[n:6][c:7]2[cH:8][cH:9][cH:10][cH:11][c:12]2[n:13][c:14]1[O:15][CH:16]([c:17]1[cH:18][cH:19][cH:20][n:21][cH:22]1)[C:23]([F:25])([F:24])[F:26].[O:1]=[S:2](=[O:3])([NH2:4])[C:27]([F:29])([F:28])[F:30].O=C([O-])[O-].[K+].[K+]>>[O:1]=[S:2](=[O:3])([NH:4][c:5]1[n:6][c:7]2[cH:8][cH:9][cH:10][cH:11][c:12]2[n:13][c:14]1[O:15][CH:16]([c:17]1[cH:18][cH:19][cH:20][n:21][cH:22]1)[C:23]([F:24])([F:25])[F:26])[C:27]([F:28])([F:29])[F:30]",
        ),
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
        # Identical reactions from index 102 in Schneider validation set
        # ground truth twice
        (
            "CC(=O)O.CS(C)=O.Cl[c:1]1[c:2]([O:3][CH:4]([C:5]([F:6])([F:7])[F:8])[c:9]2[cH:10][cH:11][cH:12][n:13][cH:14]2)[n:15][c:16]2[cH:17][cH:18][cH:19][cH:20][c:21]2[n:22]1.O=C([O-])[O-].[K+].[K+].[NH2:23][S:24](=[O:25])(=[O:26])[C:27]([F:28])([F:29])[F:30]>>[c:1]1([NH:23][S:24](=[O:25])(=[O:26])[C:27]([F:28])([F:29])[F:30])[c:2]([O:3][CH:4]([C:5]([F:6])([F:7])[F:8])[c:9]2[cH:10][cH:11][cH:12][n:13][cH:14]2)[n:15][c:16]2[cH:17][cH:18][cH:19][cH:20][c:21]2[n:22]1",
            "CC(=O)O.CS(C)=O.Cl[c:1]1[c:2]([O:3][CH:4]([C:5]([F:6])([F:7])[F:8])[c:9]2[cH:10][cH:11][cH:12][n:13][cH:14]2)[n:15][c:16]2[cH:17][cH:18][cH:19][cH:20][c:21]2[n:22]1.O=C([O-])[O-].[K+].[K+].[NH2:23][S:24](=[O:25])(=[O:26])[C:27]([F:28])([F:29])[F:30]>>[c:1]1([NH:23][S:24](=[O:25])(=[O:26])[C:27]([F:28])([F:29])[F:30])[c:2]([O:3][CH:4]([C:5]([F:6])([F:7])[F:8])[c:9]2[cH:10][cH:11][cH:12][n:13][cH:14]2)[n:15][c:16]2[cH:17][cH:18][cH:19][cH:20][c:21]2[n:22]1",
        ),
        # Identical mappings from index 102 in Schneider validation set
        # ground truth and prediction
        (
            "CC(=O)O.CS(C)=O.Cl[c:1]1[c:2]([O:3][CH:4]([C:5]([F:6])([F:7])[F:8])[c:9]2[cH:10][cH:11][cH:12][n:13][cH:14]2)[n:15][c:16]2[cH:17][cH:18][cH:19][cH:20][c:21]2[n:22]1.O=C([O-])[O-].[K+].[K+].[NH2:23][S:24](=[O:25])(=[O:26])[C:27]([F:28])([F:29])[F:30]>>[c:1]1([NH:23][S:24](=[O:25])(=[O:26])[C:27]([F:28])([F:29])[F:30])[c:2]([O:3][CH:4]([C:5]([F:6])([F:7])[F:8])[c:9]2[cH:10][cH:11][cH:12][n:13][cH:14]2)[n:15][c:16]2[cH:17][cH:18][cH:19][cH:20][c:21]2[n:22]1",
            "CC(=O)O.CS(C)=O.Cl[c:5]1[n:6][c:7]2[cH:8][cH:9][cH:10][cH:11][c:12]2[n:13][c:14]1[O:15][CH:16]([c:17]1[cH:18][cH:19][cH:20][n:21][cH:22]1)[C:23]([F:24])([F:25])[F:26].[O:1]=[S:2](=[O:3])([NH2:4])[C:27]([F:28])([F:29])[F:30].O=C([O-])[O-].[K+].[K+]>>[O:1]=[S:2](=[O:3])([NH:4][c:5]1[n:6][c:7]2[cH:8][cH:9][cH:10][cH:11][c:12]2[n:13][c:14]1[O:15][CH:16]([c:17]1[cH:18][cH:19][cH:20][n:21][cH:22]1)[C:23]([F:24])([F:25])[F:26])[C:27]([F:28])([F:29])[F:30]",
        ),
        # Identical mappings from index 102 in Schneider validation set
        # with modified Fluorine maps, ground truth and prediction
        (
            "CC(=O)O.CS(C)=O.Cl[c:1]1[c:2]([O:3][CH:4]([C:5]([F:6])([F:7])[F:8])[c:9]2[cH:10][cH:11][cH:12][n:13][cH:14]2)[n:15][c:16]2[cH:17][cH:18][cH:19][cH:20][c:21]2[n:22]1.O=C([O-])[O-].[K+].[K+].[NH2:23][S:24](=[O:25])(=[O:26])[C:27]([F:28])([F:29])[F:30]>>[c:1]1([NH:23][S:24](=[O:25])(=[O:26])[C:27]([F:28])([F:29])[F:30])[c:2]([O:3][CH:4]([C:5]([F:6])([F:7])[F:8])[c:9]2[cH:10][cH:11][cH:12][n:13][cH:14]2)[n:15][c:16]2[cH:17][cH:18][cH:19][cH:20][c:21]2[n:22]1",
            "CC(=O)O.CS(C)=O.Cl[c:5]1[n:6][c:7]2[cH:8][cH:9][cH:10][cH:11][c:12]2[n:13][c:14]1[O:15][CH:16]([c:17]1[cH:18][cH:19][cH:20][n:21][cH:22]1)[C:23]([F:25])([F:24])[F:26].[O:1]=[S:2](=[O:3])([NH2:4])[C:27]([F:29])([F:28])[F:30].O=C([O-])[O-].[K+].[K+]>>[O:1]=[S:2](=[O:3])([NH:4][c:5]1[n:6][c:7]2[cH:8][cH:9][cH:10][cH:11][c:12]2[n:13][c:14]1[O:15][CH:16]([c:17]1[cH:18][cH:19][cH:20][n:21][cH:22]1)[C:23]([F:24])([F:25])[F:26])[C:27]([F:28])([F:29])[F:30]",
        ),
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
def test_equivalent_mappings_no_can(gt, pred):
    assert are_atom_maps_equivalent(gt, pred, canonicalize=False) is True


@pytest.mark.parametrize(
    "gt, pred",
    [
        # Indices mixup
        ("[CH3:2][OH:1]>>[CH3:1][O-:2]", "[CH3:9][OH:8]>>[CH3:9][O-:8]"),
        # Mapped vs unmapped
        ("[CH3:1][OH:2]>>[CH3:1][O-:2]", "CO>>C[O-]"),
        # Different structures (non-equivalent chemistry)
        ("[CH3:1][OH:2]>>[CH3:1][O-:2]", "[CH3:1][NH2:2]>>[CH3:1][NH:2]"),
        # Wrong mapping from index 110 in Schneider validation set
        (
            "C1CCOC1.CCN(CC)CC.Cl[C:1](=[O:2])[O:3][c:4]1[cH:5][cH:6][c:7]([N+:8](=[O:9])[O-:10])[cH:11][cH:12]1.[O:13]=[C:14]([O:15][C@@H:16]1[CH2:17][O:18][C@@H:19]2[C@H:20]([OH:21])[CH2:22][O:23][C@H:24]12)[N:25]1[CH2:26][CH2:27][CH:28]([O:29][N+:30](=[O:31])[O-:32])[CH2:33][CH2:34]1>>[C:1](=[O:2])([O:3][c:4]1[cH:5][cH:6][c:7]([N+:8](=[O:9])[O-:10])[cH:11][cH:12]1)[O:21][C@H:20]1[C@H:19]2[O:18][CH2:17][C@@H:16]([O:15][C:14](=[O:13])[N:25]3[CH2:26][CH2:27][CH:28]([O:29][N+:30](=[O:31])[O-:32])[CH2:33][CH2:34]3)[C@H:24]2[O:23][CH2:22]1",
            "C1CCOC1.CCN(CC)CC.Cl[C:2](=[O:1])[O:3][c:4]1[cH:5][cH:6][c:7]([N+:8](=[O:9])[O-:10])[cH:11][cH:12]1.[O:13]([C@@H:14]1[CH2:15][O:16][C@H:17]2[C@@H:18]1[O:19][CH2:20][C@H:21]2[OH:22])[C:23](=[O:24])[N:25]1[CH2:26][CH2:27][CH:28]([O:29][N+:30](=[O:31])[O-:32])[CH2:33][CH2:34]1>>[O:1]=[C:2]([O:3][c:4]1[cH:5][cH:6][c:7]([N+:8](=[O:9])[O-:10])[cH:11][cH:12]1)[O:13][C@@H:14]1[CH2:15][O:16][C@H:17]2[C@@H:18]1[O:19][CH2:20][C@H:21]2[O:22][C:23](=[O:24])[N:25]1[CH2:26][CH2:27][CH:28]([O:29][N+:30](=[O:31])[O-:32])[CH2:33][CH2:34]1",
        ),
        # Wrong mapping - index mixup in azide
        # from index 136 in Schneider validation set
        (
            "CN(C)C=O.[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[N:8]1[CH2:9][CH2:10][CH:11]([N:12]2[c:13]3[cH:14][cH:15][cH:16][cH:17][c:18]3[O:19][c:20]3[cH:21][c:22]([C:23]#[N:24])[cH:25][cH:26][c:27]32)[CH2:28][CH2:29]1.[Cl-].[N-:30]=[N+:31]=[N-:32].[NH4+].[Na+]>>[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[N:8]1[CH2:9][CH2:10][CH:11]([N:12]2[c:13]3[cH:14][cH:15][cH:16][cH:17][c:18]3[O:19][c:20]3[cH:21][c:22](-[c:23]4[n:24][n:30][n:31][nH:32]4)[cH:25][cH:26][c:27]32)[CH2:28][CH2:29]1",
            "[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[N:8]1[CH2:9][CH2:10][CH:11]([N:12]2[c:13]3[cH:14][cH:15][cH:16][cH:17][c:18]3[O:19][c:20]3[cH:21][c:22]([C:23]#[N:24])[cH:28][cH:29][c:30]32)[CH2:31][CH2:32]1.CN(C)C=O.[Cl-].[N+:25](=[N-:26])=[N-:27].[NH4+].[Na+]>>[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[N:8]1[CH2:9][CH2:10][CH:11]([N:12]2[c:13]3[cH:14][cH:15][cH:16][cH:17][c:18]3[O:19][c:20]3[cH:21][c:22](-[c:23]4[n:24][n:25][n:26][nH:27]4)[cH:28][cH:29][c:30]32)[CH2:31][CH2:32]1",
        ),
    ],
)
def test_non_equivalent_mappings(gt, pred):
    assert are_atom_maps_equivalent(gt, pred) is False


@pytest.mark.parametrize(
    "gt, pred",
    [
        # Indices mixup
        ("[CH3:2][OH:1]>>[CH3:1][O-:2]", "[CH3:9][OH:8]>>[CH3:9][O-:8]"),
        # Mapped vs unmapped
        ("[CH3:1][OH:2]>>[CH3:1][O-:2]", "CO>>C[O-]"),
        # Different structures (non-equivalent chemistry)
        ("[CH3:1][OH:2]>>[CH3:1][O-:2]", "[CH3:1][NH2:2]>>[CH3:1][NH:2]"),
        # Wrong mapping (index 110 in Schneider validation set)
        (
            "C1CCOC1.CCN(CC)CC.Cl[C:1](=[O:2])[O:3][c:4]1[cH:5][cH:6][c:7]([N+:8](=[O:9])[O-:10])[cH:11][cH:12]1.[O:13]=[C:14]([O:15][C@@H:16]1[CH2:17][O:18][C@@H:19]2[C@H:20]([OH:21])[CH2:22][O:23][C@H:24]12)[N:25]1[CH2:26][CH2:27][CH:28]([O:29][N+:30](=[O:31])[O-:32])[CH2:33][CH2:34]1>>[C:1](=[O:2])([O:3][c:4]1[cH:5][cH:6][c:7]([N+:8](=[O:9])[O-:10])[cH:11][cH:12]1)[O:21][C@H:20]1[C@H:19]2[O:18][CH2:17][C@@H:16]([O:15][C:14](=[O:13])[N:25]3[CH2:26][CH2:27][CH:28]([O:29][N+:30](=[O:31])[O-:32])[CH2:33][CH2:34]3)[C@H:24]2[O:23][CH2:22]1",
            "C1CCOC1.CCN(CC)CC.Cl[C:2](=[O:1])[O:3][c:4]1[cH:5][cH:6][c:7]([N+:8](=[O:9])[O-:10])[cH:11][cH:12]1.[O:13]([C@@H:14]1[CH2:15][O:16][C@H:17]2[C@@H:18]1[O:19][CH2:20][C@H:21]2[OH:22])[C:23](=[O:24])[N:25]1[CH2:26][CH2:27][CH:28]([O:29][N+:30](=[O:31])[O-:32])[CH2:33][CH2:34]1>>[O:1]=[C:2]([O:3][c:4]1[cH:5][cH:6][c:7]([N+:8](=[O:9])[O-:10])[cH:11][cH:12]1)[O:13][C@@H:14]1[CH2:15][O:16][C@H:17]2[C@@H:18]1[O:19][CH2:20][C@H:21]2[O:22][C:23](=[O:24])[N:25]1[CH2:26][CH2:27][CH:28]([O:29][N+:30](=[O:31])[O-:32])[CH2:33][CH2:34]1",
        ),
        # Wrong mapping - index mixup in azide
        # from index 136 in Schneider validation set
        (
            "CN(C)C=O.[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[N:8]1[CH2:9][CH2:10][CH:11]([N:12]2[c:13]3[cH:14][cH:15][cH:16][cH:17][c:18]3[O:19][c:20]3[cH:21][c:22]([C:23]#[N:24])[cH:25][cH:26][c:27]32)[CH2:28][CH2:29]1.[Cl-].[N-:30]=[N+:31]=[N-:32].[NH4+].[Na+]>>[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[N:8]1[CH2:9][CH2:10][CH:11]([N:12]2[c:13]3[cH:14][cH:15][cH:16][cH:17][c:18]3[O:19][c:20]3[cH:21][c:22](-[c:23]4[n:24][n:30][n:31][nH:32]4)[cH:25][cH:26][c:27]32)[CH2:28][CH2:29]1",
            "[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[N:8]1[CH2:9][CH2:10][CH:11]([N:12]2[c:13]3[cH:14][cH:15][cH:16][cH:17][c:18]3[O:19][c:20]3[cH:21][c:22]([C:23]#[N:24])[cH:28][cH:29][c:30]32)[CH2:31][CH2:32]1.CN(C)C=O.[Cl-].[N+:25](=[N-:26])=[N-:27].[NH4+].[Na+]>>[CH3:1][C:2]([CH3:3])([CH3:4])[O:5][C:6](=[O:7])[N:8]1[CH2:9][CH2:10][CH:11]([N:12]2[c:13]3[cH:14][cH:15][cH:16][cH:17][c:18]3[O:19][c:20]3[cH:21][c:22](-[c:23]4[n:24][n:25][n:26][nH:27]4)[cH:28][cH:29][c:30]32)[CH2:31][CH2:32]1",
        ),
    ],
)
def test_non_equivalent_mappings_no_can(gt, pred):
    assert are_atom_maps_equivalent(gt, pred, canonicalize=False) is False


def test_invalid_smiles():
    with pytest.raises(ValueError):
        are_atom_maps_equivalent(
            "INVALID>>[CH3:1][OH:2]", "[CH3:1][OH:2]>>[CH3:1][O-:2]"
        )
