import warnings; warnings.filterwarnings("ignore")
import sys
try:
    import piaso; print("FAIL: piaso present"); sys.exit(1)
except ImportError: pass
import scanpy as sc, numpy as np, pandas as pd, laris
adata = sc.read_10x_h5(sys.argv[1] if len(sys.argv)>1 else "/tmp/piaso-data-cache/e18_v3_nuclei.h5"); adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_genes=200); sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=30); sc.pp.neighbors(adata)
sc.tl.leiden(adata, key_added="Leiden", flavor="igraph", n_iterations=2, directed=False)
adata = adata[:1500].copy()
adata.obsm["X_spatial"] = np.random.default_rng(0).random((adata.n_obs,2))*500
adata.obs["CellTypes"] = pd.Categorical(adata.obs["Leiden"].astype(str))
lrdb = laris.datasets.lrDatabase(species="mouse")
present = set(adata.var_names)
lrdb_f = lrdb[lrdb["ligand"].isin(present) & lrdb["receptor"].isin(present)].copy()
la = laris.tl.prepareLRInteraction(adata, lr_df=lrdb_f, number_nearest_neighbors=10, use_rep_spatial="X_spatial")
res = laris.tl.runLARIS(la, adata=adata, groupby="CellTypes", n_permutations=50, n_top_lr=300, calculate_pvalues=True)
res = res[0] if isinstance(res, tuple) else res
print("LARIS-only OK: rows", len(res), "(no piaso)")
