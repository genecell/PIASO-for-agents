# Workflow — ligand–receptor / cell–cell communication

Two tools in this ecosystem infer ligand–receptor interactions; **pick by whether the data is
spatial** before writing any code. Both runnable paths below are executed and passing.

## Decision rule — read this first

Key off **spatial coordinates**:

- **Spatial coordinates present** — Visium / Slide-seq / MERFISH / Xenium / Stereo-seq; coordinates
  in `adata.obsm['spatial']` or `adata.obsm['X_spatial']`; the request mentions spatial neighborhoods
  or tissue location → **LARIS** (`laris.tl.runLARIS`). LARIS bundles CellChatDB.
- **No coordinates — dissociated single-cell RNA-seq** (expression + cell-type labels only) →
  **PIASO SCALAR** (`piaso.tl.runSCALAR`). SCALAR needs a **user-supplied** LR-pair list and a
  specificity matrix; there is no bundled DB.

Routing a dissociated dataset to LARIS (its spatial kNN graph is meaningless without coordinates) or
a spatial one to SCALAR is the failure mode to avoid. When both are possible, coordinates win → LARIS.

## Install

```bash
pip install piaso-tools laris cosg "matplotlib<3.9"
```

`laris` pulls `cosg`; install `piaso-tools` for the SCALAR path. Both paths assume a clustered
AnnData with cell-type / cluster labels (e.g. `adata.obs["Leiden"]` from
`end_to_end_scrnaseq.md`) and a normalized layer such as `adata.layers["infog"]`.

---

## Path A — SCALAR (dissociated single-cell)

`piaso.tl.runSCALAR` scores interaction = (ligand specificity in sender) × (receptor specificity in
receiver), with an expression-matched permutation null. You supply **two things**: a
`specificity_matrix` (a genes × cell-types DataFrame) and an `lr_pairs` DataFrame with
`ligand`/`receptor` columns whose genes exist in the data.

- **specificity_matrix:** the concrete, tested recipe is a z-scored per-cluster mean (below). You
  *can* build it from COSG scores, but note `cosg.cosg` writes only the **top-N genes per cluster**
  into `adata.uns['cosg']` (a recarray), so you must first assemble those into a full
  genes × cell-types matrix yourself — the z-scored mean is simpler and covers all genes.
- **lr_pairs:** there is no bundled DB, but you don't have to hand-write one — if `laris` is
  installed you can reuse its CellChatDB as a ready-made list:
  `lr = laris.datasets.lrDatabase(species="mouse")[["ligand","receptor"]]`, then filter to genes
  present in `adata` (as below).

```python
import numpy as np, pandas as pd, piaso
# specificity_matrix: genes x cell types — z-scored per-cluster mean (tested recipe)
ct = adata.obs["Leiden"].astype(str)
cell_types = sorted(ct.unique(), key=int)
X = adata.layers["infog"]
means = np.vstack([np.asarray(X[(ct==c).values].mean(0)).ravel() for c in cell_types]).T
spec = pd.DataFrame(means, index=adata.var_names, columns=cell_types)
spec = spec.sub(spec.mean(1), axis=0).div(spec.std(1).replace(0,1), axis=0)
# lr_pairs: DataFrame with ligand/receptor columns, genes present in adata
lr = pd.DataFrame([{"ligand":l,"receptor":r} for l,r in
                   [("Nrxn1","Nlgn1"),("Nrxn3","Nlgn1"),("Efna5","Epha4")]
                   if l in adata.var_names and r in adata.var_names])
res = piaso.tl.runSCALAR(adata, specificity_matrix=spec, lr_pairs=lr, n_permutations=200, random_seed=42)
# DataFrame: ligand, receptor, sender, receiver, interaction_score, p_value, p_value_fdr, nlog10_p_value_fdr
```
**Out:** a DataFrame of sender→receiver interactions with scores and BH-FDR p-values (no AnnData
mutation). Visualize with `piaso.pl.plotLigandReceptorInteraction` / `plotLigandReceptorLollipop`.

**Before/after:** requires cluster labels + a normalized layer; produces a ranked interaction table.
Permutation cost scales with (pairs × `n_permutations`).

---

## Path B — LARIS (spatial)

`laris.tl.runLARIS` diffuses ligand × receptor products over a spatial kNN graph. It has two setup
gotchas that this block handles explicitly.

**Precondition:** `adata.obsm["X_spatial"]` must hold spatial coordinates. LARIS uses the key
`X_spatial` (not scanpy's `spatial`) — set or pass it accordingly. The tutorial fixtures in this hub
are dissociated (no coordinates), so use a real spatial object here.

```python
import laris, pandas as pd
lrdb = laris.datasets.lrDatabase(species="mouse")            # bundled CellChatDB (3105 mouse pairs)
# GOTCHA 1: filter LR-DB to pairs whose ligand AND receptor are both present in the data
present = set(adata.var_names)
lrdb_f = lrdb[lrdb["ligand"].isin(present) & lrdb["receptor"].isin(present)].copy()
# GOTCHA 2: the groupby column must be a categorical dtype
adata.obs["CellTypes"] = pd.Categorical(adata.obs["Leiden"].astype(str))
lr_adata = laris.tl.prepareLRInteraction(adata, lr_df=lrdb_f, number_nearest_neighbors=10,
                                         use_rep_spatial="X_spatial")
res = laris.tl.runLARIS(lr_adata, adata=adata, groupby="CellTypes",
                        n_permutations=100, n_top_lr=500, calculate_pvalues=True)
# returns DataFrame (or tuple); cols: ligand, receptor, score, Rank
```
**Out:** an interaction table (with `by_celltype=True`, a tuple whose cell-type result has
`sender, receiver, ligand, receptor, interaction_score, p_value, p_value_fdr`). Visualize with
`laris.pl.plotCCCHeatmap` / `plotCCCNetwork` / `plotCCCDotPlot`.

**Before/after:** requires spatial coordinates in `X_spatial`, cluster labels as a **categorical**,
and an LR-DB **filtered to present genes** (both gotchas above; skipping either breaks the run). Pick
`species="mouse"` vs `"human"` so DB gene symbols match `adata.var_names`.

---

## Summary

| | SCALAR | LARIS |
|---|---|---|
| Data | dissociated scRNA-seq | spatial transcriptomics |
| Coordinates | none | required (`obsm['X_spatial']`) |
| LR pairs | user-supplied | bundled CellChatDB |
| Function | `piaso.tl.runSCALAR` | `laris.tl.runLARIS` |
