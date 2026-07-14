# PIASO-data — fixtures for every code block

Every runnable code block in this hub loads from **PIASO-data**, the ecosystem's data repository
(`github.com/genecell/PIASO-data`). It has two halves: **tutorial datasets** hosted on Zenodo, and
**genome reference files** committed directly in the repo. Code blocks should use the small tutorial
fixtures below so they stay cheap to run.

## Zenodo record

- Record: <https://zenodo.org/records/19699639>
- DOI: **10.5281/zenodo.19699639**

Each tutorial dataset is a single file on that record. Fetch it by the direct content URL pattern:

```
https://zenodo.org/api/records/19699639/files/<filename>/content
```

You can download programmatically (`piaso.data.load_dataset(id)` / `fetch_dataset(id)`, cached under
`~/.piaso/data/datasets/`) or just fetch the URL directly with `curl`/`wget`/`requests`.

## Smallest fixtures — use these in code blocks

| Purpose | id | filename | size | reference |
|---|---|---|---|---|
| Loadable AnnData (real scRNA object) | `e18_v3_nuclei` | `SC3_v3_NextGem_DI_Nuclei_5K_SC3_v3_NextGem_DI_Nuclei_5K_count_sample_feature_bc_matrix.h5` | **20.25 MB** (20,250,624 B) | 10x Genomics public (E18 mouse brain nuclei, 5K, v3.1) |
| PIASOmarkerDB marker CSV | `piaso_markerdb_allen_immune` | `PIASOmarkerDB_AllenHumanImmuneHealthAtlas_L2_251219.csv` | **117 KB** (117,350 B) | Gong et al. Nature 648, 696–706 (2025) |

`e18_v3_nuclei` is a 10x `.h5` (~5,000 cells) — the smallest real expression matrix — and is what
the workflow blocks load. Fetch and load it:

```bash
curl -L -o e18_v3_nuclei.h5 \
  "https://zenodo.org/api/records/19699639/files/SC3_v3_NextGem_DI_Nuclei_5K_SC3_v3_NextGem_DI_Nuclei_5K_count_sample_feature_bc_matrix.h5/content"
```
```python
import scanpy as sc
adata = sc.read_10x_h5("e18_v3_nuclei.h5")   # PIASO-data fixture, ~5k cells
adata.var_names_make_unique()
```

`piaso_markerdb_allen_immune` is a 117 KB CSV — the smallest fixture overall — for offline marker
work. Note that the **live PIASOmarkerDB REST API** (used by `queryPIASOmarkerDB` / `analyzeMarkers`)
is a separate remote service and needs internet; this CSV is a published static slice, not the API.

md5 checksums (for verification): `e18_v3_nuclei` = `81a6ceb41e2def93ac0d0f824a610849`;
`piaso_markerdb_allen_immune` = `d4177960c47f995562ad572bb8a5f9f7`.

## Larger datasets — exist, but avoid in tests

The record also holds full atlases that are **too large for routine code blocks / CI** — do not use
them in tests:

- `sea_ad_mtg_20k` — SEA-AD MTG 20K human snRNA `.h5ad`, **1.92 GB** (Gabitto et al. Nat Neurosci 2024).
- `adult_cortex_multiome_rna` — Adult Mouse Cortex Multiome RNA `.h5ad`, **2.66 GB** (Bravo
  González-Blas et al. Nat Methods 2023).

Mid-size 10x `.h5` options (all mouse/human scRNA, tens of MB) also exist —
`mouse_brain_10k_gemx` (68.7 MB), `e18_v3_cell` (47.6 MB), `e18_v4_cell` (67.6 MB),
`pbmc_multiome_san1` (76.7 MB), `pbmc_multiome_san2` (87.8 MB) — but `e18_v3_nuclei` at 20 MB is the
lightest and is the default fixture here. All tutorial data is scRNA/snRNA: **no spatial and no
standalone ATAC matrices are shipped**, so spatial workflows (LARIS) must supply their own
coordinates.

## Genome references (committed in-repo)

Separately, `hg38/` and `mm10/` directories hold ENCODE/UCSC-derived annotation supports (gene
bodies, promoters, cCRE CTCF sites, TSS, chrom sizes; ~17 MB / ~11 MB), fetched via
`piaso.data.fetch_genome("hg38"|"mm10")`. These support ATAC/epigenomic analyses and are not needed
by the RNA workflows in this hub.

## License

PIASO-data has **no LICENSE file**; the README states genome files derive from public UCSC/ENCODE
annotations and tutorial datasets are **redistributed under CC BY 4.0 with attribution to original
sources**. Cite the original dataset paper (the `reference` column above) when using a fixture.
