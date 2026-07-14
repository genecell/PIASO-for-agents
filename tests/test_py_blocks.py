"""Phase 3.1 (Python) — functional tests + code-block compile checks.

Two layers:
  1. Functional: run the real PIASO/COSG/LARIS/Emergene/SCALAR/markerDB pipelines
     end-to-end on a PIASO-data fixture and assert the documented output keys.
     This is what proves the canonical code actually runs.
  2. Compile check: extract every ```python block from canonical/*.md and compile()
     it, catching syntax drift in the docs without needing every block's runtime state.

Run:  pytest tests/test_py_blocks.py -v
Network: functional markerDB + fixture download need internet (piaso.org, zenodo.org).
Skip the heavy functional test with:  PIASO_SKIP_FUNCTIONAL=1 pytest ...
"""
from __future__ import annotations
import os, re, warnings
from pathlib import Path
import pytest

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent.parent
CANON = ROOT / "canonical"
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fixtures  # noqa: E402


# --------------------------------------------------------------- compile check
def _python_blocks(md: Path) -> list[str]:
    return re.findall(r"```python\n(.*?)```", md.read_text(), re.DOTALL)


CANON_MD = sorted(CANON.rglob("*.md"))


@pytest.mark.parametrize("md", CANON_MD, ids=lambda p: str(p.relative_to(CANON)))
def test_python_blocks_compile(md: Path):
    """Every python block in canonical/ must at least be syntactically valid."""
    for i, block in enumerate(_python_blocks(md)):
        try:
            compile(block, f"{md.name}#py{i}", "exec")
        except SyntaxError as e:  # pragma: no cover
            pytest.fail(f"{md.relative_to(CANON)} block {i} does not compile: {e}\n{block}")


# ------------------------------------------------------------------ functional
skip_func = pytest.mark.skipif(
    os.environ.get("PIASO_SKIP_FUNCTIONAL") == "1", reason="functional tests disabled"
)


@pytest.fixture(scope="module")
def adata():
    import scanpy as sc
    import piaso, cosg  # noqa: F401
    h5 = fixtures.get("e18_v3_nuclei")
    ad = sc.read_10x_h5(str(h5))
    ad.var_names_make_unique()
    sc.pp.filter_cells(ad, min_genes=200)
    sc.pp.filter_genes(ad, min_cells=3)
    ad.layers["counts"] = ad.X.copy()
    piaso.tl.infog(ad, layer=None, n_top_genes=2000, key_added="infog")
    piaso.tl.runSVDLazy(ad, layer="infog", infog_layer="counts", n_components=50,
                        n_top_genes=2000, key_added="X_svd", random_state=1927)
    sc.pp.neighbors(ad, use_rep="X_svd", n_neighbors=15)
    sc.tl.leiden(ad, resolution=1.0, key_added="Leiden", flavor="igraph",
                 n_iterations=2, directed=False)
    ad.X = ad.layers["infog"]
    cosg.cosg(ad, groupby="Leiden", key_added="cosg", n_genes_user=30, mu=1.0)
    return ad


@skip_func
def test_infog_svd_cosg(adata):
    assert "infog" in adata.layers
    assert adata.obsm["X_svd"].shape[1] == 50
    assert "cosg" in adata.uns
    assert adata.obs["Leiden"].nunique() > 1


@skip_func
def test_gdr(adata):
    import piaso
    piaso.tl.runGDR(adata, groupby="Leiden", n_gene=30, mu=1.0, layer="infog",
                    score_layer="infog", scoring_method="scanpy", key_added="X_gdr")
    assert adata.obsm["X_gdr"].shape[0] == adata.n_obs


@skip_func
def test_score(adata):
    import piaso, numpy as np
    genes = [g for g in ["Gad1", "Gad2", "Slc17a7"] if g in adata.var_names] or list(adata.var_names[:3])
    piaso.tl.score(adata, gene_list=genes, layer="infog", key_added="myscore")
    assert "myscore" in adata.obs and np.isfinite(adata.obs["myscore"]).all()


@skip_func
def test_predict_marker(adata):
    import piaso
    if "X_gdr" not in adata.obsm:
        piaso.tl.runGDR(adata, groupby="Leiden", layer="infog", score_layer="infog",
                        scoring_method="scanpy", key_added="X_gdr")
    names = adata.uns["cosg"]["names"]
    mset = {cl: [names[cl][i] for i in range(len(names))] for cl in names.dtype.names}
    piaso.tl.predictCellTypeByMarker(adata, marker_gene_set=mset, score_layer="infog",
                                     use_rep="X_gdr", key_added="pred", smooth_prediction=True)
    assert "pred" in adata.obs


@skip_func
def test_scalar(adata):
    import piaso, numpy as np, pandas as pd
    ct = adata.obs["Leiden"].astype(str)
    cts = sorted(ct.unique(), key=int)
    X = adata.layers["infog"]
    means = np.vstack([np.asarray(X[(ct == c).values].mean(0)).ravel() for c in cts]).T
    spec = pd.DataFrame(means, index=adata.var_names, columns=cts)
    spec = spec.sub(spec.mean(1), axis=0).div(spec.std(1).replace(0, 1), axis=0)
    lr = pd.DataFrame([{"ligand": l, "receptor": r} for l, r in
                       [("Nrxn1", "Nlgn1"), ("Nrxn3", "Nlgn1"), ("Efna5", "Epha4")]
                       if l in adata.var_names and r in adata.var_names])
    res = piaso.tl.runSCALAR(adata, specificity_matrix=spec, lr_pairs=lr,
                             n_permutations=100, random_seed=42)
    assert {"interaction_score", "p_value_fdr", "sender", "receiver"} <= set(res.columns)


@skip_func
def test_markerdb_live_api():
    import piaso
    studies = piaso.tl.queryPIASOmarkerDB(list_studies=True)
    assert isinstance(studies, list) and len(studies) > 0
    df = piaso.tl.queryPIASOmarkerDB(cell_type="Microglia", limit=3)
    assert {"cell_type", "gene", "specificity_score"} <= set(df.columns)


@skip_func
def test_laris_spatial(adata):
    import laris, numpy as np, pandas as pd
    ad = adata[:1500].copy()
    ad.X = ad.layers["infog"]
    rng = np.random.default_rng(0)
    ad.obsm["X_spatial"] = rng.random((ad.n_obs, 2)) * 500
    ad.obs["CellTypes"] = pd.Categorical(ad.obs["Leiden"].astype(str))
    lrdb = laris.datasets.lrDatabase(species="mouse")
    present = set(ad.var_names)
    lrdb_f = lrdb[lrdb["ligand"].isin(present) & lrdb["receptor"].isin(present)].copy()
    lr_adata = laris.tl.prepareLRInteraction(ad, lr_df=lrdb_f, number_nearest_neighbors=10,
                                             use_rep_spatial="X_spatial")
    res = laris.tl.runLARIS(lr_adata, adata=ad, groupby="CellTypes",
                            n_permutations=50, n_top_lr=300, calculate_pvalues=True)
    res_df = res[0] if isinstance(res, tuple) else res
    assert len(res_df) > 0


@skip_func
def test_emergene_conditions(adata):
    import emergene, scanpy as sc, numpy as np, pandas as pd
    ad = adata.copy()
    ad.X = ad.layers["infog"]
    sc.pp.pca(ad, n_comps=30)
    ad.obs["condition"] = pd.Categorical((np.arange(ad.n_obs) % 2).astype(str))
    out = emergene.tl.runEMERGENE(ad, condition_key="condition", use_rep="X_pca",
                                  n_top_EG_genes=100, verbose=0)
    assert out is not None
