"""cosg-only environment (pip install cosg): no piaso, no scanpy — the docs claim COSG is
standalone since 1.1.0, so this reads the 10x h5 with h5py + anndata and normalizes by hand."""
import warnings; warnings.filterwarnings("ignore")
import sys
for mod in ("piaso", "scanpy"):
    try:
        __import__(mod); print(f"FAIL: {mod} present in cosg-only env"); sys.exit(1)
    except ImportError:
        pass
import h5py, numpy as np, scipy.sparse as sp, anndata as ad, cosg

path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/piaso-data-cache/e18_v3_nuclei.h5"
with h5py.File(path) as f:                                        # 10x h5 -> AnnData without scanpy
    m = f["matrix"]
    X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(m["shape"][:])).T.tocsr()
    genes = m["features"]["name"][:].astype(str); cells = m["barcodes"][:].astype(str)
adata = ad.AnnData(X=X.astype(np.float32)); adata.var_names = genes; adata.obs_names = cells
adata.var_names_make_unique()
adata = adata[np.asarray(adata.X.sum(1)).ravel() >= 500].copy()
adata = adata[:, np.asarray((adata.X > 0).sum(0)).ravel() >= 3].copy()
sf = 1e4 / np.asarray(adata.X.sum(1)).ravel()                      # normalize_total + log1p by hand
adata.layers["log1p"] = sp.diags(sf) @ adata.X; adata.layers["log1p"].data = np.log1p(adata.layers["log1p"].data)
rng = np.random.default_rng(0)
adata.obs["group"] = rng.choice([f"g{i}" for i in range(6)], adata.n_obs)   # any label column works for the API contract
cosg.cosg(adata, groupby="group", key_added="cosg", n_genes_user=30, mu=1.0, layer="log1p")
assert "cosg" in adata.uns and set(adata.uns["cosg"]) >= {"names", "scores"}
adata.layers["counts"] = adata.X.copy()
cosg.cosg(adata, groupby="group", key_added="cosg_p", n_genes_user=10, layer="counts", calculate_pvalues=True)
assert "pvals" in adata.uns["cosg_p"]
print("COSG-only OK: markers + p-values written (no piaso, no scanpy)")
