"""LARIS on real spatial data — the Slide-tags tonsil (241 MB, Zenodo 10.5281/zenodo.19981287).
Heavy (prepareLRBackground ~6 min): runs only with PIASO_HEAVY=1. The API contract on synthetic
coordinates is covered every run in tests/test_py_blocks.py::test_laris_synthetic_with_background.

Run:  PIASO_HEAVY=1 pytest tests/test_laris.py -v
"""
from __future__ import annotations
import os, sys, warnings
from pathlib import Path
import pytest

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fixtures  # noqa: E402

pytestmark = pytest.mark.skipif(os.environ.get("PIASO_HEAVY") != "1", reason="set PIASO_HEAVY=1 for the tonsil run")


def test_tonsil_end_to_end():
    import anndata as ad, laris as la
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    adata = ad.read_h5ad(str(fixtures.get("laris_tonsil")))
    assert "X_spatial" in adata.obsm and "cell_type" in adata.obs
    lr_df = la.datasets.lrDatabase(species="human")
    assert len(lr_df) == 2951
    lr_data = la.tl.prepareLRInteraction(adata, lr_df, use_rep_spatial="X_spatial")
    assert lr_data.n_obs == adata.n_obs and lr_data.n_vars > 1500          # 1985 pairs with both genes measured
    bg = la.tl.prepareLRBackground(adata, lr_df, use_rep_spatial="X_spatial", n_matched_genes=30)
    laris_lr, res = la.tl.runLARIS(lr_data, adata, use_rep="X_spatial", use_rep_spatial="X_spatial",
                                   groupby="cell_type", background=bg)
    assert set(laris_lr.columns) >= {"ligand", "receptor", "score"}
    assert {"sender", "receiver", "interaction_name", "interaction_score", "p_value", "p_value_fdr",
            "null_matchability", "null_support", "pair_breadth"} <= set(res.columns)
    sig = res[res.p_value_fdr < 0.05]
    assert 500 < len(sig) < 5000                                            # 1,372 on 2026-09-04 (tutorial: 1,345 at defaults)
    top = res.nlargest(1, "interaction_score").iloc[0]
    assert top["interaction_name"] == "FCER2::CR2"                          # published tonsil reference value
    la.pl.plotCCCHeatmap(res); plt.close("all")
    la.pl.plotCCCSpatial(lr_data, basis="X_spatial", interaction="C3::CR2", color_by="score"); plt.close("all")
    la.pl.plotCCCDotPlot(res, interactions_to_plot=["C3::CR2", "FCER2::CR2"], senders=["MRC", "B_naive"],
                         receivers=["FDC_LZDZ", "B_naive"]); plt.close("all")        # senders/receivers (or sender_receiver_pairs) are required
    la.pl.plotCCCNetwork(res, cell_type_of_interest=str(res.iloc[0]["sender"]), interaction_direction="sending",
                         data=adata, groupby="cell_type"); plt.close("all")
