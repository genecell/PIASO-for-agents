# piaso-mcp

A **local stdio MCP server** for the PIASO ecosystem. Serves **knowledge and public data
only — never compute, never user data, never your expression matrices.**

> **Generated artifact.** The doc tools read a snapshot of `canonical/` that `build.py`
> writes into `piaso_mcp/data/`. Do not hand-edit that snapshot — edit `canonical/` and
> re-run `build.py`.

## Tools

| Tool | What it does | Source |
|---|---|---|
| `search_docs(query)` | Search the ecosystem knowledge pack — works with **zero packages installed** | bundled `canonical/` |
| `get_api(function)` | Exact current signature + data-object contract | bundled `api-map.md` |
| `compare_implementations(function)` | COSG (Python) vs COSGR (R) divergences | bundled `components/cosg.md` |
| `resolve_install(components, language)` | The exact install line for a chain, per language | bundled `meta.yaml` |
| `list_datasets()` | PIASO-data fixtures (URLs + loading code) | bundled `data.md` |
| `query_marker_db(...)` / `get_markers(...)` | Cell-type marker lookup | **proxies the live PIASOmarkerDB REST API** at `piaso.org` |
| `list_studies()` | Studies in PIASOmarkerDB | proxies the live API |

The marker tools give **R / non-Python users their first programmatic path to
PIASOmarkerDB** (the `piaso.tools` query API is Python-only). They require network access
to `piaso.org`; everything else works offline.

## Run

```bash
uvx piaso-mcp            # zero-install, once published
# or from source:
pip install -e mcp/ && piaso-mcp
python -m piaso_mcp      # equivalent
```

## Configure your agent

```jsonc
// Claude Code / Cursor / Windsurf: "mcpServers"; VS Code: "servers"; Zed: "context_servers"
{ "mcpServers": { "piaso": { "command": "uvx", "args": ["piaso-mcp"] } } }
```

## License

BSD-3-Clause.
