# Ecosystem-wide gotchas

A tight reference of traps that span the PIASO ecosystem. Component-specific detail lives in
each `components/*.md`.

## Install / import

- **PIASO requires `matplotlib < 3.9` — pin it or `import piaso` fails.** PIASO 1.1.0 imports
  `piaso/plotting/color.py` at package import time, and that module calls
  `matplotlib.cm.get_cmap(...)` — an API **removed in matplotlib 3.9**. Because
  `piaso/__init__.py` loads `pl` → `color`, this runs on every `import piaso`, so under
  matplotlib ≥ 3.9 the import raises before you can do anything. PIASO's `pyproject.toml`
  only declares `matplotlib>=3.5.2` (no upper cap), so nothing enforces this for you. Always
  install with `pip install piaso-tools "matplotlib<3.9"` (the working env uses 3.8.4).
  Upstream fix would be `plt.get_cmap` / `matplotlib.colormaps`.

- **Emergene pins `annoy < 1.17.0`; PIASO `stitchSpace` (BBKNN) segfaults with newer annoy.**
  `pip install emergene` caps `annoy<1.17.0` for this reason. PIASO's `stitchSpace` uses
  BBKNN and warns hard that `annoy >= 1.17` causes a **BBKNN segfault** — pin
  `annoy==1.16.3` when running it. (`stitchSpace` also rejects `use_rep` embeddings
  containing NaN/Inf, which would otherwise segfault BBKNN.)

## Data / layer contracts

- **INFOG (`piaso.tl.infog`) needs RAW counts.** Give it raw UMI counts (from `adata.X` or a
  raw-counts layer); it rejects negatives. Keep a copy of the raw counts (e.g.
  `adata.layers['counts'] = adata.X.copy()`) before other steps overwrite `.X`. INFOG's
  output layer defaults to `infog`.

- **`piaso.tl.score` defaults to `layer='infog'` (NORMALIZED values), not raw counts.** It
  expects INFOG-normalized values and **errors if the `infog` layer is absent**. This is the
  opposite input from `infog` itself — run INFOG first, then score off its output layer.
  Genes not in `adata.var_names` are silently dropped. Single-set vs multi-set input changes
  the behavior drastically (writes to `adata` + returns None vs returns a tuple).

- **COSG expects NORMALIZED / log values in `.X` (or a layer).** Point `.X` (or `layer=`) at
  normalized values before calling `cosg.cosg(...)` — e.g. `adata.X = adata.layers['infog']`.
  Feeding raw counts gives wrong markers.

- **COSG is a dependency, not a re-export.** `piaso-tools` and `laris` both hard-depend on
  `cosg` (auto-installed) and call `cosg.cosg(...)` internally, but there is **no
  `piaso.cosg`** — to call COSG yourself you still `import cosg`.

## LARIS (spatial ligand–receptor)

- **Filter the LR database to genes actually present in the data.** `laris.datasets.lrDatabase`
  returns the full CellChatDB (mouse 3105 / human 2951 pairs). Keep only pairs whose ligand
  **and** receptor are both in `adata.var_names` before running, e.g.
  `lrdb[lrdb['ligand'].isin(present) & lrdb['receptor'].isin(present)]`.

- **The `groupby` column must be a categorical dtype.** Cast it explicitly, e.g.
  `adata.obs['CellTypes'] = pd.Categorical(adata.obs['Leiden'].astype(str))`, before calling
  `laris.tl.runLARIS`.

## PIASOmarkerDB

- **PIASOmarkerDB is a REMOTE REST API — it needs internet.** The client
  (`piaso.tl.queryPIASOmarkerDB` / `getMarkers` / `analyzeMarkers` / the `PIASOmarkerDB`
  class) issues HTTP calls to `https://piaso.org/piasomarkerdb/api/v1/`; it is **not
  wheel-bundled data**. It is Python-only, requires `requests`, and caches under
  `~/.piaso/markers`. Any offline/bundled use (e.g. an MCP tool) must proxy the live API or
  obtain a raw snapshot + redistribution license from the maintainers first.
