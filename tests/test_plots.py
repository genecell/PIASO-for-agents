"""Phase 3.1 — headless execution of the piaso.pp / piaso.pl helpers.

Plotting can't be asserted visually in CI, but it CAN be run under a non-interactive
(Agg) backend to prove the calls execute without error against a real object. This
covers the preprocessing utilities and the plotting suite documented in
components/piaso.md (which the main functional suite skips because they need a display).

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
    import scanpy as sc, numpy as np, pandas as pd
    import piaso  # noqa: F401
    ad = sc.read_10x_h5(str(fixtures.get("e18_v3_nuclei")))
    ad.var_names_make_unique()
    sc.pp.filter_cells(ad, min_genes=200)
    sc.pp.filter_genes(ad, min_cells=3)
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    sc.pp.pca(ad, n_comps=30)
    sc.pp.neighbors(ad)
    sc.tl.leiden(ad, key_added="Leiden", flavor="igraph", n_iterations=2, directed=False)
    sc.tl.umap(ad)
    ad.obs["grp2"] = pd.Categorical((np.arange(ad.n_obs) % 3).astype(str))
    ad.obsm["X_spatial"] = np.random.default_rng(0).random((ad.n_obs, 2))
    return ad


def test_pp_helpers(adata):
    import piaso, pandas as pd
    assert isinstance(piaso.pp.table(adata.obs["Leiden"], as_dataframe=True), pd.DataFrame)
    cc = piaso.pp.getCrossCategories(adata.obs, "Leiden", "grp2")
    assert isinstance(cc, pd.Categorical)
    piaso.pp.rotateSpatialCoordinates(adata, angle_degrees=90, spatial_key="X_spatial", inplace=True)
    assert adata.obsm["X_spatial"].shape[1] == 2


def test_pl_helpers(adata):
    import piaso
    gene = next(g for g in adata.var_names if adata[:, g].X.sum() > 0)
    piaso.pl.plot_embeddings_split(adata, color=gene, splitby="Leiden", basis="X_umap",
                                   show_figure=False)
    plt.close("all")
    piaso.pl.plot_features_violin(adata, feature_list=[gene], groupby="Leiden", show_figure=False)
    plt.close("all")
    out = piaso.pl.plotConfusionMatrix(adata, groupby_query="Leiden", groupby_reference="grp2",
                                       return_objects=True)
    plt.close("all")
    assert isinstance(out, tuple)
    cmap = piaso.pl.createCustomCmapFromHex(["#000000", "#ff0000", "#ffffff"])
    assert cmap is not None
