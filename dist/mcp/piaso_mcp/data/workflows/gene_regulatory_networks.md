# Workflow — gene regulatory networks / regulons (cytorete)

From an annotated dataset (AnnData or `.cytome`) to TF → target **regulons** with per-cell
activity and per-type specificity, on RNA alone. Executed on the SEA-AD cortex cytome (20,000
nuclei, `Subclass` labels) with cytorete 0.1.1; mirrors the executed piaso.org tutorial
*cytorete: RNA regulon inference*. The genome-free parts (scoring a supplied regulon dict, building
a cistrome from sequences) also run on the small `e18_v3_nuclei` fixture in the hub's tests.

## Install

```bash
pip install 'cytorete[motif]'      # cytorete + py2bit (needed to read the .2bit genome); pulls piaso-tools, cosg, cytome
```

## Decision rule — which tool

"Which **genes** define this cell type" → COSG (`components/cosg.md`). "Which **TFs** drive it" →
cytorete. cytorete needs: a cell-type / cluster column, a genome (`hg38` or `mm10`: `.2bit`
sequence + TSS BED via `piaso.data`), a motif database (JASPAR / CIS-BP), and — on a targeted
spatial panel — the TF *and* its targets measured. It does **not** need ATAC (that is the withheld
`inferGRN` chain).

## Step 1 — Inputs, all fetchable

```python
import numpy as np, pandas as pd
import piaso, cytorete
jaspar = piaso.data.fetch_jaspar()            # JASPAR CORE vertebrates (MEME); or fetch_cisbp()
twobit = piaso.data.fetch_2bit("hg38")        # ~800 MB, cached under ~/.piaso/data
piaso.data.fetch_genome("hg38")               # TSS BED (+ promoters, cCREs)
ds = piaso.data.load_dataset("sea_ad_mtg_20k_cytome", return_type="cytome")   # 20k nuclei, 'Subclass'; any annotated AnnData works too
```
**Out:** paths + an open dataset. Species must match the genome (`"mm10"` for mouse; gene symbols
are matched case-sensitively).

## Step 2 — One call

```python
TFS = ["SOX9", "OLIG2", "SPI1", "IRF8", "MEF2C", "NEUROD2", "DLX2", "LHX6",
       "TCF7L2", "NFIB", "RORB", "CUX2", "FEZF2", "PAX6", "EGR1"]      # omit tf_list to run every expressed TF with a motif
cytorete.inferRegulon(ds, "hg38", "Subclass", jaspar_path=jaspar, twobit_path=twobit, tf_list=TFS)
```
**Out (written onto the object):** embeddings `X_regulon` (cells × regulons, activity) and
`X_regulon_pval`; `ds.metadata["regulon"]` (or `adata.uns["regulon"]`) with `names`, `regulons`
(TF → targets), `weights`, `per_celltype`, `edges`, `specificity`, `params`. Run log on the
cortex: 2552 genes with promoters; cistrome 15 TFs × 2552 genes, 3411 edges; 15 global regulons
(median 45 targets); per-cell-type regulons for 23 types; ~50 s of compute after loading.

What decided the edges: strand-aware promoter windows (−1000/+500 around each TSS) scanned with
a per-motif background (`cistrome_method="motif_bg"`; `"nes"` for RcisTarget-style NES); a
motif-supported TF→gene pair survives only if the TF's and target's COSG specificity profiles
agree across the `groupby` types (positive-sign trans co-specificity); TFs keep ≥ `min_targets`
targets. Activity is `piaso.tl.score` against size- and expression-matched control sets, hence
the per-cell p-value.

## Step 3 — Read the biology

```python
md = ds.metadata["regulon"]; names = md["names"]                          # ['CUX2', 'DLX2', 'EGR1', 'FEZF2', 'IRF8', 'LHX6', ...]
A = np.asarray(ds.embeddings["X_regulon"]); P = np.asarray(ds.embeddings["X_regulon_pval"])   # (20000, 15)
cytorete.pl.regulonActivity(ds, groupby="Subclass")                        # regulon x cell-type heatmap
cytorete.pl.regulonEmbedding(ds, regulons=["CUX2", "FEZF2", "LHX6", "IRF8"], basis="X_umap")   # where each is active
cytorete.pl.regulonNetwork(ds, tf="CUX2", max_targets=20)                 # what is in it — check the targets before believing it
cytorete.pl.regulonSpecificityScatter(ds, groupby="Subclass", key="X_regulon")   # is the TF's own expression driving it? (default key is X_grn — pass key=)
spec = cytorete.regulonSpecificity(ds, groupby="Subclass", copy=True)     # long form: cell_type, regulon, cosg_score, rank
spec.pivot(index="regulon", columns="cell_type", values="cosg_score")     # pivot before a heatmap
```
**Out:** figures; `spec` (219 rows on the cortex run). Expected controls: SPI1/IRF8 → microglia,
OLIG2 → OPC/oligodendrocytes, SOX9/PAX6 → astrocytes, LHX6 → MGE interneurons, FEZF2 → deep
layers, CUX2 → L2/3 IT. Read **activity and p-value together**: −log10 p saturates at 3.0 with
1,000 control sets; a small, sharp regulon (IRF8) is significant in few cells and that is the
correct answer, not a weak one.

## Step 4 — Per-cell values for any other plot

```python
for tf in ["CUX2", "IRF8"]:
    j = names.index(tf)
    ds.cells[f"act_{tf}"] = A[:, j]
    ds.cells[f"nlp_{tf}"] = -np.log10(np.clip(P[:, j], 1e-300, 1))
ds.flush()
piaso.pl.plotEmbedding(ds, color="act_CUX2", basis="X_umap", vmin_pct=5, vmax_pct=95)
piaso.pl.plotEmbedding(ds, color="nlp_CUX2", basis="X_umap", vmin=0, vmax=3)
```
On an AnnData: `adata.obs[f"act_{tf}"] = adata.obsm["X_regulon"][:, j]`.

## Variant — score regulons you already have (no genome needed)

```python
regulons = {"GABA_TFs": ["Dlx1", "Dlx5", "Lhx6", "Sox6"], "Glut_TFs": ["Neurod2", "Neurod6", "Tbr1", "Satb2"]}
cytorete.regulonActivity(adata, regulons={k: [g for g in v if g in adata.var_names] for k, v in regulons.items()},
                         score_layer="infog", key_added="X_regulon")        # PIASOscore per set + p-values -> obsm['X_regulon'], ['X_regulon_pval']
cytorete.regulonSpecificity(adata, groupby="leiden")                       # COSG on the activity matrix
```
This is `piaso.tl.score` with regulon bookkeeping — useful for published regulon lists.

## Variant — spatial sections and development

Same `inferRegulon` on a spatial cytome (`"mm10"`, `groupby="annotation"`, `cosg_layer="infog",
score_layer="infog"`), then `regulonEmbedding(ds, regulons=[...], basis="spatial")` draws
activity on the tissue; across stages, compare `X_regulon` per stage. Check panel coverage first
(`spatial_transcriptomics.md`).

## Not in this release

`cytorete.inferGRN` / `inferGRN_consensus` (multiome RNA + ATAC), `inferTFActivity` (ATAC) and
`pp.build_peak_cistrome` raise `ImportError: ... not part of this distribution` at call time.
