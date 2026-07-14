import warnings; warnings.filterwarnings("ignore")
import sys
try:
    import piaso; print("FAIL: piaso present"); sys.exit(1)
except ImportError: pass
import scanpy as sc, numpy as np, pandas as pd, emergene
adata = sc.read_10x_h5(sys.argv[1] if len(sys.argv)>1 else "/tmp/piaso-data-cache/e18_v3_nuclei.h5"); adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_genes=200); sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=30)
adata.obs["condition"] = pd.Categorical((np.arange(adata.n_obs)%2).astype(str))
out = emergene.tl.runEMERGENE(adata, condition_key="condition", use_rep="X_pca", n_top_EG_genes=100, verbose=0)
print("EMERGENE-only OK:", type(out).__name__, "(no piaso)")
