# Emergene — component reference (self-sufficient)

Emergene performs **individual-cell differential transcriptomics across conditions** (disease vs
control, developmental stages, perturbations): instead of one fold change per gene per cluster,
it scores the condition signature **per cell**, so a response confined to 20 % of a population
stays visible instead of averaging into a weak response everywhere. Works on scRNA-seq and
spatial data. Standalone package (does not depend on `piaso-tools`), but designed to run
**downstream of a PIASO analysis**: it consumes the embedding and normalized layer you already
have. Tutorial: https://piaso.org/tutorials/emergene/ · docs: https://genecell.github.io/Emergene/.

## Install

```bash
pip install emergene     # pins annoy<1.17.0 (BBKNN segfault guard); pulls scanpy + bbknn
```

Executed against **emergene 1.0.2** (unchanged since 2026-02). `runEMERGENE` raises an
`ImportError` with fix instructions if it detects annoy ≥ 1.17.

## Import / public surface

`import emergene` (or `import emergene as eg`). `emergene.tl` = `emergene.tools`, plus `.pp`, `.pl`.

- `eg.tl`: `runEMERGENE`, `runMarkG`, `score`, `identifyGeneModule`
- `eg.pp`: `infog`, `convertTopGeneDictToDF`

## What `runEMERGENE` computes

`runEMERGENE(adata, condition_key='Sample', use_rep='X_pca', use_rep_acrossDataset='X_pca',
layer=None, n_top_EG_genes=500, ...)` builds cross-condition connectivity with **BBKNN** (batching
on the condition). Per condition it computes a target specificity (cosine of expression vs its
within-condition diffusion), a shuffled-graph random background, and a cross-condition background
(diffusion from the *other* conditions); the Emergene score is `GSP − mu·random_GSP −
beta·condition_GSP`. It selects the top `n_top_EG_genes` per condition and writes per-cell local
fold changes.

Required input state:
- An embedding in **`adata.obsm[use_rep]`** — and the same (or another) key in
  **`use_rep_acrossDataset`**, which defaults to `'X_pca'` separately; **pass both** when your
  embedding is `X_gdr` / `X_svd`, or you get `'X_pca' not found`.
- Condition labels in `adata.obs[condition_key]` with ≥ 2 conditions (warns and points to
  `runMarkG` if only one).
- Expression in `adata.X` or `layer=` (INFOG or log-normalized).
- Side effect: always writes `adata.layers['localFC']`, even with `inplace=False`.

## Executed block (downstream of the PIASO pipeline)

```python
import emergene as eg
# adata: PIASO output with layers['infog'], obsm['X_gdr'], obs['condition'] (>= 2 levels)
out = eg.tl.runEMERGENE(adata, condition_key="condition", use_rep="X_gdr", use_rep_acrossDataset="X_gdr",
                        layer="infog", n_top_EG_genes=100)
# returns (dict per condition, DataFrame); single condition -> eg.tl.runMarkG instead
```

**Why `use_rep="X_gdr"`:** cells are matched to comparable cells across conditions *in the
embedding*. A variance-driven embedding is often dominated by the perturbation itself, so
matching in it is circular; GDR's marker-guided space is less so. `layer="infog"` keeps the
normalization consistent with the rest of the analysis. Two modes: **per-condition** (each
condition vs control — cleaner signatures, and one comparable gene-weight vector per condition,
so cosine similarity between conditions gives a data-driven hierarchy of perturbations) and
**all-in-one** (each condition vs all others).

## Decision rule — `runEMERGENE` vs `runMarkG`

- **≥ 2 conditions** to contrast → `runEMERGENE` (needs `condition_key`; BBKNN cross-condition diffusion).
- **Single condition** (marker / spatially variable genes, no comparison) → `runMarkG(adata, use_rep=, layer=)`.

## Citation

Emergene has **no paper of its own** — cite Wu, Dai *et al.*, *Nature* (2026):

> Wu, S.J., Dai, M. et al. Pyramidal neurons proportionately alter cortical interneuron
> subtypes. *Nature* (2026). DOI: 10.1038/s41586-025-09996-8
