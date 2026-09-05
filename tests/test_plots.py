"""Headless execution of the piaso.pp / piaso.pl calls the canonical docs use.

Plotting can't be asserted visually in CI, but it CAN be run under a non-interactive (Agg)
backend to prove the calls execute against a real object built with the PIASO-native pipeline.

Run:  pytest tests/test_plots.py -q     (needs internet for the fixture download)
      PIASO_SKIP_FUNCTIONAL=1 skips it.
"""
from __future__ import annotations
import os, warnings, sys
from pathlib import Path
import pytest
import matplotlib
matplotlib.use("Agg")  # headless — no display needed
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fixtures  # noqa: E402

pytestmark = pytest.mark.skipif(
    os.environ.get("PIASO_SKIP_FUNCTIONAL") == "1", reason="functional tests disabled"
)


@pytest.fixture(scope="module")
def adata():
    import numpy as np, pandas as pd, piaso, cosg
    piaso.settings.set_figure_params(style="cell")
    ad = piaso.pp.read_10x_h5(str(fixtures.get("e18_v3_nuclei")))
    piaso.pp.calculateCellMetrics(ad, prefix_vars={"mt": "mt-", "ribo": ["Rps", "Rpl"]})
    piaso.pp.filter_cells(ad, min_counts=500, min_features=250)
    piaso.tl.infog(ad, n_top_genes=2000)
    piaso.tl.runSVD(ad, layer="infog", n_components=30, key_added="X_svd")
    piaso.tl.neighbors(ad, use_rep="X_svd"); piaso.tl.leiden(ad, key_added="leiden"); piaso.tl.umap(ad, use_rep="X_svd")
    cosg.cosg(ad, groupby="leiden", key_added="cosg", n_genes_user=10, layer="infog")
    ad.obs["grp2"] = pd.Categorical((np.arange(ad.n_obs) % 3).astype(str))
    ad.obsm["X_spatial"] = np.random.default_rng(0).random((ad.n_obs, 2))
    return ad


def test_pp_helpers(adata):
    import piaso, pandas as pd
    assert isinstance(piaso.pp.table(adata.obs["leiden"], as_dataframe=True), pd.DataFrame)
    assert isinstance(piaso.pp.getCrossCategories(adata.obs, "leiden", "grp2"), pd.Categorical)
    piaso.pp.rotateSpatialCoordinates(adata, angle_degrees=90, spatial_key="X_spatial", inplace=True)
    assert adata.obsm["X_spatial"].shape[1] == 2
    gm = piaso.pp.calculateGroupMetrics(adata, groupby="leiden")
    assert "n_cells" in gm.columns


def test_pl_suite(adata):
    import piaso
    top = []
    for c in adata.obs["leiden"].cat.categories:
        top += [g for g in adata.uns["cosg"]["names"][c][:2] if g not in top]
    piaso.pl.embedding(adata, basis="X_umap", color="leiden", legend_loc="both", show=False); plt.close("all")
    piaso.pl.embedding(adata, basis="X_umap", color=["n_counts", "pct_counts_mt"], ncol=2, show=False); plt.close("all")
    piaso.pl.dotplot(adata, top[:16], groupby="leiden", standard_scale="var", show=False); plt.close("all")
    piaso.pl.violin(adata, ["n_genes", "pct_counts_mt"], groupby="leiden", show=False); plt.close("all")
    piaso.pl.scatter(adata, x="n_counts", y="n_genes", color="pct_counts_mt", logx=True, logy=True, marginals=True, show=False); plt.close("all")
    piaso.pl.stackedBarplot(adata, groupby="leiden", splitby="grp2", show=False); plt.close("all")
    piaso.pl.plot_features_violin(adata, feature_list=["n_genes"], groupby="leiden", show_figure=False); plt.close("all")
    piaso.pl.plot_embeddings_split(adata, color="leiden", splitby="grp2", basis="X_umap", show_figure=False); plt.close("all")
    out = piaso.pl.plotConfusionMatrix(adata, groupby_query="leiden", groupby_reference="grp2", return_objects=True); plt.close("all")
    assert isinstance(out, tuple)
    piaso.pl.sankey(adata, left="leiden", right="grp2", show=False); plt.close("all")
    gm = piaso.pp.calculateGroupMetrics(adata, groupby="leiden")
    piaso.pl.plotGroupMetrics(gm, data=adata, groupby="leiden", show=False); plt.close("all")
    assert piaso.pl.createCustomCmapFromHex(["#000000", "#ff0000", "#ffffff"]) is not None
