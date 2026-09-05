# piaso-mcp

<!-- Ownership marker for the official MCP Registry (verifies this PyPI package is ours). -->
mcp-name: io.github.genecell/piaso-mcp

A **local stdio MCP server** for the PIASO single-cell omics ecosystem — PIASO, COSG / COSGR,
cytome (Python + R), LARIS, Emergene, cytorete, PIASO-data. Serves **knowledge and public data
only — never compute, never user data, never your expression matrices.**

> **Generated artifact.** The doc tools read a snapshot of `canonical/` that `build.py` writes
> into `piaso_mcp/data/`. Do not hand-edit that snapshot — edit `canonical/` and re-run `build.py`.

## Tools

| Tool | What it does | Source |
|---|---|---|
| `search_docs(query)` | Search the ecosystem knowledge pack — works with **zero packages installed** | bundled `canonical/` |
| `get_api(function)` | Data-object contract (reads/writes, defaults) for a function; points at the live generated API reference | bundled `components/*.md` |
| `compare_implementations(function)` | COSG (Python) vs COSGR (R); `"cytome"` for the Python vs R cytome packages | bundled |
| `resolve_install(components, language)` | Exact install line(s) per language, incl. extras and the R routes (r-universe, conda-forge); never a matplotlib pin | bundled `meta.yaml` |
| `list_tutorials(topic, component)` | The **piaso.org tutorial index** (executed tutorials): title, URL, what it covers | bundled `meta.yaml` |
| `version_matrix()` | Versions this snapshot was tested against, roles, citations | bundled `meta.yaml` |
| `check_versions()` | Current PyPI release of each package vs the tested version — disclose drift | **PyPI JSON (network)** |
| `list_datasets()` / `get_dataset(name)` | The **PIASO-data registry** (id, format, size, cells, counts layer, md5, how to load) | **live `datasets.json` from GitHub**, cached 24 h; offline fallback to bundled `data.md` |
| `query_marker_db(...)` / `get_markers(...)` | Cell-type marker lookup | **live PIASOmarkerDB REST API** at `piaso.org` |
| `list_studies()` | Studies in PIASOmarkerDB (36) | live API |

The marker tools give **R / non-Python users their first programmatic path to PIASOmarkerDB** (the
`piaso.tl` client is Python-only). The registry, version and marker tools need network access;
everything else works offline.

## Run

```bash
uvx piaso-mcp            # zero-install
# or from source:
pip install -e mcp/ && piaso-mcp
python -m piaso_mcp      # equivalent
```

## Configure your agent

```jsonc
// Claude Code / Cursor / Windsurf: "mcpServers"; VS Code: "servers"; Zed: "context_servers"; Codex: TOML [mcp_servers.piaso]
{ "mcpServers": { "piaso": { "command": "uvx", "args": ["piaso-mcp"] } } }
```

Per-client details (uv prerequisite, PATH, Codex TOML, Zed shape) are in the repository README.

## Versions

`piaso-mcp 0.1.0` — snapshot tested against piaso-tools 1.2.3 · cosg 1.2.0 · cytome 0.3.1 ·
laris 0.13.0 · emergene 1.0.2 · cytorete 0.1.1 · COSGR 1.0.0 · cytome (R) 0.1.0 (2026-09-04).
`0.0.x` served the July 2026 snapshot (five packages, scanpy-based workflows).

## License

BSD-3-Clause.
