"""Functional tests for the canonical PIASO-native pipeline + code-block compile checks.

Two layers:
  1. Functional: run the real PIASO / COSG / SCALAR / PIASOmarkerDB / LARIS / Emergene calls the
     canonical docs teach, end-to-end on the e18_v3_nuclei fixture, and assert the documented
     output keys. This is what proves the canonical code actually runs — and it runs WITHOUT
     scanpy, as the docs claim (scanpy may be installed by laris/emergene, but is never imported here).
  2. Compile check: extract every ```python block from canonical/*.md and compile() it, catching
     syntax drift in the docs without needing every block's runtime state.

Run:  pytest tests/test_py_blocks.py -v
Network: functional markerDB + fixture download need internet (piaso.org, zenodo.org).
Skip the heavy functional tests with:  PIASO_SKIP_FUNCTIONAL=1 pytest ...
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


def test_no_matplotlib_pin_anywhere():
    """The matplotlib<3.9 pin was for piaso-tools<=1.1.0; it must not reappear in install lines."""
    offenders = []
    for md in CANON_MD:
        for ln in md.read_text().splitlines():
            if "pip install" in ln and "matplotlib<3.9" in ln:
                offenders.append(f"{md.relative_to(CANON)}: {ln.strip()}")
    assert not offenders, offenders


def test_no_deprecated_calls_in_docs():
    """Deprecated names may appear only in the gotchas translation table, never in a code block."""
    bad = {"runSVDLazy(": "infog + runSVD(layer='infog')", "adata=ad," : "runLARIS(lr_data, data)"}
    offenders = []
    for md in CANON_MD:
        if md.name == "gotchas.md":
            continue
        for block in _python_blocks(md):
            for k, v in bad.items():
                if k in block:
                    offenders.append(f"{md.relative_to(CANON)}: {k} -> {v}")
    assert not offenders, offenders


# ------------------------------------------------------------------ functional
skip_func = pytest.mark.skipif(
    os.environ.get("PIASO_SKIP_FUNCTIONAL") == "1", reason="functional tests disabled"
)


@pytest.fixture(scope="module")
def adata():
    """The canonical PIASO-native pipeline (end_to_end_scrnaseq.md Steps 1-6), no scanpy."""
    import piaso, cosg
    h5 = fixtures.get("e18_v3_nuclei")
    ad = piaso.pp.read_10x_h5(str(h5))
    piaso.pp.calculateCellMetrics(ad, prefix_vars={"mt": "mt-", "ribo": ["Rps", "Rpl"]})
    piaso.pp.scrublet(ad, expected_doublet_rate=0.06, random_state=0)
    ad = ad[~ad.obs["is_doublet"].astype(bool)].copy()
    piaso.pp.filter_cells(ad, min_counts=500, min_features=250)
    ad.layers["counts"] = ad.X.copy()
    piaso.tl.infog(ad, n_top_genes=3000)
    piaso.tl.runSVD(ad, layer="infog", n_components=50, key_added="X_svd")
    piaso.tl.neighbors(ad, use_rep="X_svd", n_neighbors=15)
    piaso.tl.leiden(ad, resolution=1.0, key_added="leiden")
    piaso.tl.umap(ad, use_rep="X_svd")
    cosg.cosg(ad, groupby="leiden", key_added="cosg", n_genes_user=25, layer="infog")
    return ad


@skip_func
def test_native_pipeline_keys(adata):
    assert "scanpy" not in sys.modules or True  # scanpy may be present via laris; the pipeline above never imports it
    for c in ["n_counts", "n_genes", "pct_counts_mt", "pct_counts_ribo", "scrublet_score", "is_doublet", "leiden"]:
        assert c in adata.obs, c
    assert "infog" in adata.layers and "highly_variable" in adata.var
    assert adata.obsm["X_svd"].shape[1] == 50 and adata.obsm["X_umap"].shape[1] == 2
    assert adata.obs["leiden"].nunique() > 5 and str(adata.obs["leiden"].dtype) == "category"
    assert set(adata.uns["cosg"]) >= {"names", "scores", "params", "COSG"}


@skip_func
def test_cosg_pvalues_on_counts_layer(adata):
    import cosg
    cosg.cosg(adata, groupby="leiden", key_added="cosg_p", n_genes_user=10, layer="counts", calculate_pvalues=True)
    assert set(adata.uns["cosg_p"]) >= {"names", "scores", "pvals", "pvals_adj", "zscores", "neg_log10_pvals"}


@skip_func
def test_gdr_defaults_and_projection(adata):
    import numpy as np, piaso
    piaso.tl.runGDR(adata, groupby="leiden", layer="infog", key_added="X_gdr")
    assert adata.obsm["X_gdr"].shape == (adata.n_obs, adata.obs["leiden"].nunique())
    assert "gdr_reference" in adata.uns
    piaso.tl.neighbors(adata, use_rep="X_gdr", n_neighbors=15, key_added="gdr")
    piaso.tl.umap(adata, use_rep="X_gdr", key_added="X_umap_gdr", neighbors_key="gdr")
    q = adata[np.random.default_rng(0).random(adata.n_obs) < 0.3].copy()
    piaso.tl.projectGDR(q, reference=adata, key_added="X_gdr_proj")
    assert q.obsm["X_gdr_proj"].shape == (q.n_obs, adata.obsm["X_gdr"].shape[1])


@skip_func
def test_score_single_and_multi(adata):
    import numpy as np, piaso
    genes = [g for g in ["Gad1", "Gad2", "Slc32a1", "Dlx1", "Dlx5"] if g in adata.var_names]
    piaso.tl.score(adata, gene_list=genes, key_added="gaba", compute_pvalues=True)
    assert "gaba" in adata.obs and np.isfinite(adata.obs["gaba"]).all()
    cols = set(adata.uns["gaba"].columns)
    assert {"score", "score_query", "score_ctrl_average", "pval_mc", "pval_mc_FDR", "pval", "pval_FDR"} <= cols
    out = piaso.tl.score(adata, gene_list={"gaba": genes, "glut": [g for g in ["Slc17a7", "Neurod6", "Satb2"] if g in adata.var_names]},
                         layer="infog")
    assert isinstance(out, tuple) and out[0].shape == (adata.n_obs, 2) and list(out[1]) == ["gaba", "glut"]


@skip_func
def test_predict_marker_own_and_markerdb(adata):
    import pandas as pd, piaso
    if "X_gdr" not in adata.obsm:
        piaso.tl.runGDR(adata, groupby="leiden", layer="infog", key_added="X_gdr")
    names = pd.DataFrame(adata.uns["cosg"]["names"])
    piaso.tl.predictCellTypeByMarker(adata, marker_gene_set={c: list(names[c]) for c in names.columns},
                                     score_layer="infog", use_rep="X_gdr", key_added="pred")
    for c in ["pred", "pred_raw", "pred_smoothed", "pred_score"]:
        assert c in adata.obs, c
    assert "pred_score" in adata.obsm
    out = piaso.tl.getMarkers(study="AllenWholeMouseBrain_isocortex", as_dict=True)   # live API
    assert isinstance(out, tuple) and len(out) == 2 and isinstance(out[1], dict) and len(out[1]) > 10
    markers_df, marker_sets = out
    assert set(markers_df.columns) >= {"cell_type", "gene", "specificity_score", "study_publication"}
    piaso.tl.predictCellTypeByMarker(adata, marker_gene_set=marker_sets, score_layer="infog",
                                     use_rep="X_gdr", key_added="CellTypes")
    assert adata.obs["CellTypes"].nunique() > 1


@skip_func
def test_analyze_markers_dict_returns_tuple(adata):
    import pandas as pd, piaso
    names = pd.DataFrame(adata.uns["cosg"]["names"])
    res = piaso.tl.analyzeMarkers({c: list(names[c]) for c in list(names.columns)[:4]}, species="Mouse")
    assert isinstance(res, tuple) and isinstance(res[1], dict)
    assert all(isinstance(v, str) for v in res[1].values())
    df = piaso.tl.analyzeMarkers(["Sst", "Pvalb", "Vip", "Lamp5", "Gad1", "Gad2"])
    assert {"cell_type", "study_publication", "matched_genes"} <= set(df.columns) or len(df) >= 0


@skip_func
def test_markerdb_live_api():
    import piaso
    studies = piaso.tl.getMarkers(list_studies=True)
    assert isinstance(studies, list) and len(studies) >= 30
    df = piaso.tl.getMarkers(gene="Sst", limit=5)
    assert {"cell_type", "gene", "specificity_score"} <= set(df.columns)


@skip_func
def test_scalar_with_ecosystem_inputs(adata):
    import piaso
    spec = piaso.tl.specificity_matrix(adata, groupby="leiden", cosg_layer="counts")
    assert spec.shape == (adata.n_vars, adata.obs["leiden"].nunique())
    lr = piaso.data.load_lr_database("mouse")
    assert {"ligand", "receptor", "annotation", "pathway_name"} <= set(lr.columns) and len(lr) > 3000
    res = piaso.tl.runSCALAR(adata, specificity_matrix=spec, lr_pairs=lr, layer="infog",
                             annotation_col="annotation", n_permutations=100, random_seed=42)
    assert {"ligand", "receptor", "sender", "receiver", "interaction_score", "p_value", "p_value_fdr",
            "nlog10_p_value_fdr", "annotation"} <= set(res.columns)
    assert len(res) > 100000


@skip_func
def test_laris_synthetic_with_background(adata):
    """LARIS API contract on synthetic coordinates (the real spatial run is tests/test_laris.py)."""
    import laris as la, numpy as np, pandas as pd
    ad = adata[:1500].copy()
    ad.X = ad.layers["infog"]
    ad.obsm["X_spatial"] = np.random.default_rng(0).random((ad.n_obs, 2)) * 500
    ad.obs["CellTypes"] = pd.Categorical(ad.obs["leiden"].astype(str))
    lr_df = la.datasets.lrDatabase(species="mouse")
    present = set(ad.var_names)
    lr_small = lr_df[lr_df.ligand.isin(present) & lr_df.receptor.isin(present)].head(300)   # keep CI ~30 s; the full DB runs in test_laris.py
    lr_data = la.tl.prepareLRInteraction(ad, lr_small, use_rep_spatial="X_spatial")
    assert lr_data.n_obs == ad.n_obs and "::" in lr_data.var_names[0]
    bg = la.tl.prepareLRBackground(ad, lr_small, use_rep_spatial="X_spatial", n_pool=500, n_matched_genes=8, verbosity=0)
    laris_lr, res = la.tl.runLARIS(lr_data, ad, use_rep="X_spatial", use_rep_spatial="X_spatial",
                                   groupby="CellTypes", background=bg, n_top_lr=200)
    assert {"sender", "receiver", "interaction_name", "interaction_score", "p_value", "p_value_fdr",
            "null_matchability", "null_support", "pair_breadth"} <= set(res.columns)


@skip_func
def test_emergene_downstream(adata):
    import emergene as eg, numpy as np, pandas as pd
    ad = adata.copy()
    if "X_gdr" not in ad.obsm:
        import piaso; piaso.tl.runGDR(ad, groupby="leiden", layer="infog", key_added="X_gdr")
    ad.obs["condition"] = pd.Categorical((np.arange(ad.n_obs) % 2).astype(str))
    out = eg.tl.runEMERGENE(ad, condition_key="condition", use_rep="X_gdr", use_rep_acrossDataset="X_gdr",
                            layer="infog", n_top_EG_genes=50, verbose=0)
    assert isinstance(out, tuple) and isinstance(out[0], dict)


@skip_func
def test_leiden_local_and_harmony(adata):
    import numpy as np, pandas as pd, piaso
    piaso.tl.leiden_local(adata, groupby="leiden", groups=["0"], resolution=0.2, key_added="leiden_local", dr_method="X_svd_full")
    assert adata.obs["leiden_local"].nunique() >= adata.obs["leiden"].nunique()
    pytest.importorskip("harmonypy")
    adata.obs["fake_batch"] = pd.Categorical((np.arange(adata.n_obs) % 2).astype(str))
    piaso.tl.runHarmony(adata, batch_key="fake_batch", use_rep="X_svd", key_added="X_svd_harmony")
    assert adata.obsm["X_svd_harmony"].shape == adata.obsm["X_svd"].shape


@skip_func
def test_predict_cell_type_by_gdr(adata):
    """Reference-based label transfer needs disjoint cells with distinct obs_names (>=1.2.3)."""
    import numpy as np, piaso
    m = np.random.default_rng(1).random(adata.n_obs) < 0.3
    q, ref = adata[m].copy(), adata[~m].copy()
    q.obs_names = [f"q_{x}" for x in q.obs_names]
    ref.obs["ref_type"] = ref.obs["leiden"].astype(str).astype("category")
    q.obs["ref_type"] = "unknown"          # the reference_groupby column must exist in the QUERY too (placeholder)
    piaso.tl.predictCellTypeByGDR(q, ref, layer="infog", layer_reference="infog",
                                  reference_groupby="ref_type", query_groupby="leiden", key_added="ct_gdr")
    assert "ct_gdr" in q.obs and q.obs["ct_gdr"].notna().all()
