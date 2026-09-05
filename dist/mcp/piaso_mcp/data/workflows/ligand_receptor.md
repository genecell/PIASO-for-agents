# Workflow — ligand–receptor / cell–cell communication

Two tools in this ecosystem infer ligand–receptor interactions; **pick by whether the data is
spatial** before writing any code. They share the database. Path A is executed on the
`e18_v3_nuclei` fixture, Path B on the LARIS Slide-tags tonsil object (see `data.md`).

## Decision rule — read this first

Key off **spatial coordinates**:

- **Coordinates present** (Visium / Visium HD / Xenium / MERFISH / Stereo-seq / Slide-tags;
  `adata.obsm['spatial']` or `['X_spatial']`; the request mentions neighbourhoods, tissue,
  where) → **LARIS**. Interaction is scored **per cell against its spatial neighbours**, so a pair
  only scores where the partners are adjacent; the sender→receiver summary is derived from that
  and can be drawn on the section.
- **No coordinates — dissociated scRNA / snRNA-seq** (expression + cell-type labels) → **SCALAR**
  (`piaso.tl.runSCALAR`). Every ordered cell-type pair is scored as if contact were possible.
- **Same database**: CellChatDB via `piaso.data.load_lr_database(species)` or
  `laris.datasets.lrDatabase(species)` (human 2951 / mouse 3105 pairs, `annotation` column =
  mechanism). A pair list curated for one transfers to the other.
- **Natural pairing**: SCALAR on the dissociated reference to find which interactions exist in the
  tissue at all, LARIS on the spatial section to ask where they happen.
- **Targeted panel? Count first**: an LR pair needs *both* genes measured (313-gene Xenium: 11 of
  2951; 5,006-gene: 1,637; whole-transcriptome in situ: 2,831; the tonsil: 1,985).

Routing a dissociated dataset to LARIS (its spatial graph is meaningless without coordinates)
or a spatial one to SCALAR is the failure mode to avoid. When both are possible, coordinates win.

## Install

```bash
pip install piaso-tools laris      # SCALAR lives in piaso-tools; laris pulls cosg + scanpy
```

---

## Path A — SCALAR (dissociated single-cell) *(executed)*

`piaso.tl.runSCALAR` scores `interaction = specificity[ligand, sender] × specificity[receptor,
receiver]` with a permutation null from expression-matched control genes
(`n_nearest_neighbors=30`, `n_permutations=1000`) and BH FDR **per sender–receiver pair**. Both
inputs come from the ecosystem: the specificity matrix from COSG (all genes) via
`piaso.tl.specificity_matrix`, the pair list from CellChatDB via `piaso.data`.

Prerequisites: an AnnData with a normalized layer (`layers["infog"]`), a raw-counts layer (here
`layers["counts"]`; the fixture has counts in `.X`), and a cell-type / cluster column.

```python
import numpy as np, pandas as pd
import piaso, cosg
# adata from end_to_end_scrnaseq.md: layers['infog'], obs['leiden']
adata.layers["counts"] = adata.X.copy()                                   # specificity_matrix reads raw counts from cosg_layer
spec = piaso.tl.specificity_matrix(adata, groupby="leiden", cosg_layer="counts")   # (n_genes, n_groups) COSG lambda, all genes
lr = piaso.data.load_lr_database("mouse")                                 # CellChatDB (3105, 28); "human" -> 2951; annotation= to slice a mechanism
res = piaso.tl.runSCALAR(adata, specificity_matrix=spec, lr_pairs=lr, layer="infog",
                         annotation_col="annotation", n_permutations=1000, random_seed=42)
# DataFrame (no adata mutation): ligand, receptor, sender, receiver, interaction_score, p_value, p_value_fdr, nlog10_p_value_fdr, annotation
sig = res[res["p_value_fdr"] < 0.05]
sig.groupby(["sender", "receiver"]).size().sort_values(ascending=False).head()
```
**Out:** a ranked interaction table (1.47 M rows on the fixture: ~3,000 usable pairs × 22² ordered
cluster pairs). With Leiden clusters on 6k nuclei nothing reaches FDR < 0.05 — SCALAR wants real
cell types with sharp markers; on the 20k-nucleus SEA-AD cortex the tutorial finds VLMC → astrocyte
collagen signalling, ICAM2 → ITGAM/ITGB2 endothelium → microglia, and 1,102 significant
interneuron → excitatory interactions. **Counts track abundance and marker sharpness** — compare
pairs of comparable size before saying one type "talks more".

Plot one pair with the two specificities that produced it, and split by mechanism:

```python
sig = sig.copy()
sig["CellTypeXCellType"] = sig["sender"] + "@" + sig["receiver"]
sig["ligandXreceptor"] = sig["ligand"] + "-->" + sig["receptor"]
sig["ligand_specificity"] = [spec.at[r.ligand, r.sender] for r in sig.itertuples()]
sig["receptor_specificity"] = [spec.at[r.receptor, r.receiver] for r in sig.itertuples()]
piaso.pl.plotLigandReceptorInteraction(interactions_df=sig, specificity_df=spec,
                                       cell_type_pairs=["VLMC@Astrocyte"], ligand_receptor_sep="-->",
                                       top_n=30, y_max=float(np.ceil(sig["interaction_score"].max() * 10) / 10),
                                       heatmap_cmap="Purples", shared_legend=True)
ecm = sig[sig["annotation"] == "ECM-Receptor"]                            # same sender, one mechanism -> a controlled comparison
# piaso.pl.plotLigandReceptorLollipop(sig, cell_type_pairs=[...], col_cell_type_pair="CellTypeXCellType", sort_by_category=True)
```
**Set `y_max` from the data** (default 10; scores top out below 1). The specificity matrix is
the raw COSG λ: scores compare *within* a sender–receiver pair, not across pairs.

---

## Path B — LARIS (spatial) *(executed on the tonsil)*

Three calls: diffuse the pairs over the spatial graph, build the (reusable) matched-gene null, run.

**Preconditions:** coordinates in `adata.obsm["X_spatial"]` (LARIS's key — copy from `spatial` if
needed), a cell-type column with no missing values, log-normalized `.X`.

```python
import anndata as ad, laris as la
adata = ad.read_h5ad("adata_tonsil.h5ad")                                 # 5,695 cells, obsm['X_spatial'], obs['cell_type']
# adata.obsm["X_spatial"] = np.asarray(adata.obsm["spatial"], dtype=np.float64)   # when coordinates are under scanpy's key
lr_df = la.datasets.lrDatabase(species="human")                           # 2951 pairs; pairs with absent genes are dropped with a warning
lr_data = la.tl.prepareLRInteraction(adata, lr_df, use_rep_spatial="X_spatial")      # new AnnData (5695 x 1985 LR pairs), var 'ligand::receptor'
bg = la.tl.prepareLRBackground(adata, lr_df, use_rep_spatial="X_spatial")            # the expensive step (minutes); reusable across groupby / sweeps -> pickle it
laris_lr, res = la.tl.runLARIS(lr_data, adata, use_rep="X_spatial", use_rep_spatial="X_spatial",
                               groupby="cell_type", background=bg)                   # returns a TUPLE
res[res.p_value_fdr < 0.05].nsmallest(6, "p_value")[["sender", "receiver", "interaction_name", "interaction_score", "p_value_fdr"]]
```
**Out:** `laris_lr` (spatially specific LR pairs, scored) and `res`, one row per
sender–receiver–pair: `interaction_score, p_value, p_value_fdr, nlog10_p_value_fdr` +
diagnostics `null_matchability`, `null_support`, `pair_breadth`. P-values are exact tail counts —
deterministic, floor `1/(n_matched_genes² + 1)` (1e-4 at defaults). Executed here on the tonsil:
389,060 combinations, **1,372 significant at FDR < 0.05** (`n_matched_genes=30`; the tutorial's default run
reports 1,345), with
`FCER2::CR2` B_naive→B_naive on top by score and `C3::CR2` MRC→FDC among the most significant.

```python
la.pl.plotCCCHeatmap(res)                                                            # sender x receiver: significant pairs
la.pl.plotCCCSpatial(lr_data, basis="X_spatial", interaction="C3::CR2", color_by="score")   # the per-cell score on the section — what LARIS is for
la.pl.plotCCCNetwork(res, cell_type_of_interest="FDC_LZDZ", interaction_direction="sending", data=adata, groupby="cell_type")
la.pl.plotCCCDotPlot(res, interactions_to_plot=["C3::CR2", "FCER2::CR2"], senders=["MRC", "B_naive"], receivers=["FDC_LZDZ", "B_naive"])   # senders/receivers or sender_receiver_pairs= required
```

Group by **clusters or coarse types** (LARIS tests every sender–receiver pair). Several samples
or conditions → per-sample `res` tables → `la.tl.compareLARIS(...)` (subject-level statistics,
volcano via `plotCompareLARIS`) or `compareLARISMatched` after `buildJointEmbedding` (compare at
matched cell states). Same calls on a `.cytome`. Details, version notes (p-values changed in 0.12
and 0.13; `n_permutations` is legacy) and the cross-condition family: `components/laris.md`.

---

## Summary

| | SCALAR | LARIS |
|---|---|---|
| Data | dissociated scRNA/snRNA-seq | spatial transcriptomics |
| Coordinates | none | required (`obsm['X_spatial']`) |
| Unit of the answer | (pair, sender, receiver) | per cell **and** (pair, sender, receiver) |
| What constrains a hit | specificity in both partners | specificity **and** physical proximity |
| Inputs | `piaso.tl.specificity_matrix` + `piaso.data.load_lr_database` | `laris.datasets.lrDatabase` (same CellChatDB) |
| Null | permutation over matched control genes (`n_permutations`) | exact enumeration over matched pseudo-pairs (`prepareLRBackground`) |
| Function | `piaso.tl.runSCALAR` | `laris.tl.prepareLRInteraction` → `prepareLRBackground` → `runLARIS` |
| Output can be | ranked, filtered, plotted | all of that, plus mapped onto the tissue (`plotCCCSpatial`) |
