"""cytome — the streaming path: build a .cytome from the fixture, run the same PIASO / COSG calls
on it, read results back. Covers workflows/streaming_large_data.md and components/cytome.md at
zero download cost (the 19 MB fixture).

Run:  pytest tests/test_cytome.py -v      PIASO_SKIP_FUNCTIONAL=1 skips.
"""
from __future__ import annotations
import os, sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
import pytest

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fixtures  # noqa: E402

pytestmark = pytest.mark.skipif(os.environ.get("PIASO_SKIP_FUNCTIONAL") == "1", reason="functional tests disabled")


@pytest.fixture(scope="module")
def workdir(tmp_path_factory):
    return tmp_path_factory.mktemp("cytome")


def test_from_10x_h5_streamed_pipeline(workdir):
    import cytome, piaso, cosg
    path = workdir / "e18_raw.cytome"
    ds = cytome.from_10x_h5(str(fixtures.get("e18_v3_nuclei")), str(path), sample_name="E18")   # use the RETURNED dataset
    assert ds.n_cells > 5000 and "RNA" in ds.modalities and "RNA_counts" in ds.list_matrices()
    piaso.pp.calculateCellMetrics(ds, modality="RNA", prefix_vars={"mt": "mt-", "ribo": ["Rps", "Rpl"]})
    assert {"n_counts", "n_genes", "pct_counts_mt"} <= set(ds.cells.columns)
    piaso.tl.infog(ds, modality="RNA", n_top_genes=3000, save_layer=True)
    assert "RNA_infog" in ds.list_matrices()
    piaso.tl.runSVD(ds, modality="RNA", layer="infog", n_components=30, key_added="X_svd")
    piaso.tl.neighbors(ds, use_rep="X_svd", n_neighbors=15)
    piaso.tl.leiden(ds, resolution=1.0, key_added="leiden")
    piaso.tl.umap(ds, use_rep="X_svd")
    assert "leiden" in ds.cells.columns and "X_umap" in list(ds.embeddings.keys())
    res = cosg.cosg(ds, groupby="leiden", modality="RNA", n_genes_user=25, layer="infog")   # RETURNS a dict
    assert isinstance(res, dict) and {"names", "scores", "groups_order"} <= set(res)
    names = pd.DataFrame(res["names"], columns=list(res["groups_order"]))
    assert names.shape[0] == 25
    # cells table SQL + streaming reads
    sub = ds.cells.query(f"leiden = '{res['groups_order'][0]}'")
    assert len(sub) > 0
    n = sum(1 for _ in ds.RNA.counts.iter_rows())
    assert n > 1 and ds.RNA.counts[:10, :10].shape == (10, 10)
    # subset is a copy
    mask = ds.cells.query_mask(f"leiden IN ('{res['groups_order'][0]}','{res['groups_order'][1]}')")
    piaso.pp.filter_cells(ds, mask=np.asarray(mask), inplace=False, output=str(workdir / "sub.cytome"), overwrite=True)
    sub_ds = cytome.open(str(workdir / "sub.cytome"))
    assert 0 < sub_ds.n_cells < ds.n_cells
    sub_ds.close(); ds.close()


def test_from_anndata_round_trip_and_analysis(workdir):
    import cytome, piaso, cosg
    ad = piaso.pp.read_10x_h5(str(fixtures.get("e18_v3_nuclei")))
    piaso.pp.filter_cells(ad, min_counts=500, min_features=250)
    ad.layers["counts"] = ad.X.copy()
    piaso.tl.infog(ad, n_top_genes=2000)
    piaso.tl.runSVD(ad, layer="infog", n_components=30, key_added="X_svd")
    piaso.tl.neighbors(ad, use_rep="X_svd"); piaso.tl.leiden(ad, key_added="leiden"); piaso.tl.umap(ad, use_rep="X_svd")
    path = workdir / "e18.cytome"
    ds = cytome.from_anndata(ad, modality="RNA", output=str(path))
    assert ds.n_cells == ad.n_obs
    assert "RNA_counts" in ds.list_matrices() and "RNA_infog" in ds.list_matrices()   # integer .X -> RNA_counts
    assert {"RNA_umap", "RNA_svd"} <= set(ds.embeddings.keys())                      # obsm renamed
    genes = [g for g in ["Gad1", "Gad2", "Slc32a1"] if g in ad.var_names]
    piaso.tl.score(ds, gene_list=genes, key_added="gaba", modality="RNA", layer="infog")
    assert "gaba" in ds.cells.columns
    piaso.tl.runGDR(ds, groupby="leiden", layer="infog", modality="RNA", key_added="X_gdr")
    assert "X_gdr" in list(ds.embeddings.keys())
    res = cosg.cosg(ds, groupby="leiden", modality="RNA", layer="infog", n_genes_user=10)
    assert isinstance(res, dict)
    ds.set_categories("leiden", order=sorted(ad.obs["leiden"].astype(str).unique(), key=int))
    back = ds.to_anndata(modality="RNA")
    assert back.shape == ad.shape and "X_umap" in back.obsm and "infog" in back.layers
    ds.close()
    # backed h5ad conversion + merge
    h5ad = workdir / "e18.h5ad"; ad.write_h5ad(str(h5ad))
    dsb = cytome.from_h5ad(str(h5ad), output=str(workdir / "e18_backed.cytome"), modality="RNA", backed=True, verbose=False)
    assert dsb.n_cells == ad.n_obs; dsb.close()
    mg = cytome.merge([str(path), str(workdir / "e18_backed.cytome")], output=str(workdir / "merged.cytome"),
                      batch_key="sample_id", force=True)
    assert mg.n_cells >= ad.n_obs; mg.close()
