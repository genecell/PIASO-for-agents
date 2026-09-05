"""cytorete — genome-free checks always; the real RNA regulon chain only with PIASO_HEAVY=1
(needs the SEA-AD cortex cytome 269 MB + JASPAR + hg38 .2bit ~800 MB, ~5 min).

Run:  pytest tests/test_cytorete.py -v        PIASO_HEAVY=1 pytest tests/test_cytorete.py -v
"""
from __future__ import annotations
import os, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import pytest

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fixtures  # noqa: E402

skip_func = pytest.mark.skipif(os.environ.get("PIASO_SKIP_FUNCTIONAL") == "1", reason="functional tests disabled")
heavy = pytest.mark.skipif(os.environ.get("PIASO_HEAVY") != "1", reason="set PIASO_HEAVY=1 for the real regulon run")


def test_import_surface_and_withheld_chains():
    import cytorete
    assert cytorete.tl is cytorete.tools
    for n in ["inferRegulon", "regulonActivity", "regulonSpecificity", "infer_regulon", "regulon_activity"]:
        assert callable(getattr(cytorete, n))
    for n in ["regulonActivity", "regulonEmbedding", "regulonNetwork", "regulonSpecificityScatter"]:
        assert callable(getattr(cytorete.pl, n))
    # the piaso.tl names are forwarders
    import piaso
    for n in ["inferRegulon", "inferGRN", "inferTFActivity", "regulonActivity", "regulonSpecificity"]:
        assert callable(getattr(piaso.tl, n))


@skip_func
def test_regulon_activity_from_supplied_dict():
    """The genome-free half: score a hand-made regulon dict (PIASOscore) + COSG specificity."""
    import piaso, cytorete
    ad = piaso.pp.read_10x_h5(str(fixtures.get("e18_v3_nuclei")))
    piaso.pp.filter_cells(ad, min_counts=500, min_features=250)
    piaso.tl.infog(ad, n_top_genes=2000)
    piaso.tl.runSVD(ad, layer="infog", n_components=30, key_added="X_svd")
    piaso.tl.neighbors(ad, use_rep="X_svd"); piaso.tl.leiden(ad, key_added="leiden")
    regs = {"GABA_TFs": ["Dlx1", "Dlx5", "Lhx6", "Sox6", "Gad1", "Gad2"], "Glut_TFs": ["Neurod2", "Neurod6", "Tbr1", "Satb2", "Slc17a7"]}
    regs = {k: [g for g in v if g in ad.var_names] for k, v in regs.items()}
    cytorete.regulonActivity(ad, regulons=regs, score_layer="infog", key_added="X_regulon", n_ctrl_set=50, verbose=0)
    assert ad.obsm["X_regulon"].shape == (ad.n_obs, 2) and "X_regulon_pval" in ad.obsm
    spec = cytorete.regulonSpecificity(ad, groupby="leiden", copy=True, verbose=0)
    assert isinstance(spec, pd.DataFrame) and {"cell_type", "regulon", "cosg_score"} <= set(spec.columns)


@skip_func
@heavy
def test_infer_regulon_real_run():
    """The executed block of components/cytorete.md (SEA-AD cortex cytome, 15 TFs)."""
    import piaso, cytorete
    pytest.importorskip("py2bit")
    jaspar = piaso.data.fetch_jaspar(); twobit = piaso.data.fetch_2bit("hg38"); piaso.data.fetch_genome("hg38")
    ds = piaso.data.load_dataset("sea_ad_mtg_20k_cytome", return_type="cytome")
    TFS = ["SOX9", "OLIG2", "SPI1", "IRF8", "MEF2C", "NEUROD2", "DLX2", "LHX6", "TCF7L2", "NFIB", "RORB", "CUX2", "FEZF2", "PAX6", "EGR1"]
    cytorete.inferRegulon(ds, "hg38", "Subclass", jaspar_path=jaspar, twobit_path=twobit, tf_list=TFS)
    keys = list(ds.embeddings.keys())
    assert "X_regulon" in keys and "X_regulon_pval" in keys
    md = ds.metadata["regulon"]
    assert {"names", "regulons", "specificity"} <= set(md) and len(md["names"]) >= 10
    A = np.asarray(ds.embeddings["X_regulon"]); assert A.shape[0] == ds.n_cells
    spec = cytorete.regulonSpecificity(ds, groupby="Subclass", copy=True)
    assert {"cell_type", "regulon", "cosg_score", "rank"} <= set(spec.columns)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    cytorete.pl.regulonActivity(ds, groupby="Subclass", show=False); plt.close("all")
    cytorete.pl.regulonEmbedding(ds, regulons=md["names"][:4], basis="X_umap", show=False); plt.close("all")
    cytorete.pl.regulonNetwork(ds, tf=md["names"][0], max_targets=10, show=False); plt.close("all")
    ds.close()
