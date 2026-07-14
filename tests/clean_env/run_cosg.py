import warnings; warnings.filterwarnings("ignore")
import sys
try:
    import piaso; print("FAIL: piaso present in cosg-only env"); sys.exit(1)
except ImportError: pass
import scanpy as sc, cosg
adata = sc.read_10x_h5(sys.argv[1] if len(sys.argv)>1 else "/tmp/piaso-data-cache/e18_v3_nuclei.h5"); adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_genes=200); sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
sc.pp.pca(adata, n_comps=50); sc.pp.neighbors(adata)
sc.tl.leiden(adata, key_added="Leiden", flavor="igraph", n_iterations=2, directed=False)
cosg.cosg(adata, groupby="Leiden", key_added="cosg", n_genes_user=30, mu=1.0)
assert "cosg" in adata.uns
print("COSG-only OK: clusters", adata.obs["Leiden"].nunique(), "| markers written (no piaso)")
