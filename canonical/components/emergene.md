# Emergene — component reference (self-sufficient)

Emergene performs **individual-cell differential transcriptomics across conditions**
(e.g. disease vs control, developmental stages): given multiple conditions it
identifies which individual cells and genes change most between them, using graph
diffusion + cosine similarity + cross-condition background correction. It works on
scRNA-seq and spatial data. This file assumes nothing about PIASO being installed;
Emergene is a standalone package (it does not depend on `piaso-tools`).

## Install

```bash
pip install emergene     # pins annoy<1.17.0 (BBKNN segfault guard)
```

The `annoy<1.17.0` pin is a hard requirement: annoy ≥ 1.17.0 causes BBKNN segfaults,
and `runEMERGENE` raises an ImportError with fix instructions if it detects a newer
annoy.

## Import / public surface

`import emergene` (or `import emergene as eg`). Short submodule aliases:
`emergene.tl` = `emergene.tools`, `emergene.pp`, `emergene.pl`.

- `eg.tl`: `runEMERGENE`, `runMarkG`, `score`, `identifyGeneModule`
- `eg.pp`: `infog`, `convertTopGeneDictToDF`

## What `runEMERGENE` computes

`runEMERGENE(adata, condition_key='Sample', use_rep='X_pca', n_top_EG_genes=500, ...)`
builds cross-dataset connectivity with **BBKNN** (batching on the condition). Per
condition it computes a target specificity (cosine similarity of expression vs a
within-condition diffused version), a shuffled-graph random background, and a
cross-condition background (diffusion from the *other* conditions); the final
Emergene score is `GSP − mu*random_GSP − beta*condition_GSP`. It selects the top
`n_top_EG_genes` per condition and also writes per-cell local fold changes.

Required input state:

- A low-dim embedding in **`adata.obsm[use_rep]`** (default `'X_pca'`).
- Condition labels in **`adata.obs[condition_key]`** (default `'Sample'`) with ≥2
  conditions (it warns and suggests `runMarkG` if only one is found).
- Expression in `adata.X` or `adata.layers[layer]` (log-normalized or INFOG
  recommended).

Side effect: `runEMERGENE` always writes `adata.layers['localFC']`, even with
`inplace=False`.

## Verified block (verbatim)

```python
import emergene, scanpy as sc
sc.pp.pca(adata, n_comps=30)
adata.obs["condition"] = pd.Categorical(adata.obs["condition"].astype(str))  # >=2 conditions
out = emergene.tl.runEMERGENE(adata, condition_key="condition", use_rep="X_pca", n_top_EG_genes=100)
# returns (dict, DataFrame). Single-condition marker analysis -> emergene.tl.runMarkG instead.
```

## Decision rule — runEMERGENE vs runMarkG

Both live in Emergene:

- **≥2 conditions** to contrast → **`runEMERGENE`** (needs `condition_key`; uses
  BBKNN cross-condition diffusion).
- **Single condition** (just want marker / spatially-variable genes, no condition
  comparison, no BBKNN) → **`runMarkG`**.

## Citation

Emergene has **no paper of its own** — cite Wu, Dai *et al.*, *Nature* (2026):

> Wu, S.J., Dai, M. et al. Pyramidal neurons proportionately alter cortical
> interneuron subtypes. *Nature* (2026). DOI: 10.1038/s41586-025-09996-8
