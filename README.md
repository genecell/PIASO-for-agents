# PIASO-for-agents

**Make the [PIASO](https://piaso.org) single-cell omics ecosystem first-class for any
coding agent — Claude Code, Cursor, Copilot, Codex, Windsurf, Cline, Aider — from one
canonical, agent-neutral knowledge pack.**

Maintained by the PIASO maintainers (Gord Fishell Lab, Harvard Medical School / Broad
Institute). Every agent-specific format (Claude skill, Cursor rules, `AGENTS.md`,
`llms.txt`, MCP server) is a **generated artifact** built from `canonical/` — never a
hand-maintained copy. A CI drift check (`python build.py --check`) fails the build if any
`dist/` artifact is out of sync with `canonical/`, and the code-block test suite re-runs on
every component release, so the guidance cannot silently rot.

## The ecosystem

Independently-installable packages under [github.com/genecell](https://github.com/genecell):

| Component | Package | Language | Role |
|---|---|---|---|
| [PIASO](https://github.com/genecell/PIASO) | `piaso-tools` | Python (+Rust) | Umbrella single-cell toolkit — see the capability table below |
| [COSG](https://github.com/genecell/COSG) | `cosg` | Python | Fast, specific marker-gene identification |
| [COSGR](https://github.com/genecell/COSGR) | `COSG` | R | COSG for Seurat / SingleCellExperiment |
| [LARIS](https://github.com/genecell/LARIS) | `laris` | Python | Ligand–receptor interaction in spatial transcriptomics |
| [Emergene](https://github.com/genecell/Emergene) | `emergene` | Python | Individual-cell differential expression across conditions |
| [PIASO-data](https://github.com/genecell/PIASO-data) | — | data | Genome references + tutorial datasets |

Each component is **independently installable** — you can `pip install cosg` (or `laris`, or
`emergene`) on its own, so a COSG-only user is a first-class citizen. Note the dependency
direction, though: installing `piaso-tools` (and `laris`) also pulls in `cosg`, so a PIASO
user always has COSG available. The hub's unique value is documenting how the components
**compose**, and the cross-component choices no single repo can make (e.g. SCALAR vs LARIS for
ligand–receptor: spatial data → LARIS, dissociated single-cell → SCALAR).

### Inside `piaso-tools`

The `piaso` package is itself a toolkit. Full references live in
[`canonical/components/piaso.md`](canonical/components/piaso.md). Grouped by what is a
PIASO-introduced method vs. a convenience wrapper around a standard step:

**Methods introduced by PIASO**

| Capability | Entry point | What it does |
|---|---|---|
| INFOG normalization | `piaso.tl.infog` | Information-content normalization of raw UMI counts + HVG selection |
| GDR (marker-gene-guided DR) | `piaso.tl.runGDR` / `runGDRParallel` | Embedding whose axes are per-cluster COSG-marker scores; also does batch integration |
| Gene-set scoring | `piaso.tl.score` | Optimized expression-matched-control gene-set enrichment scoring — Rust-accelerated |
| Cell-type prediction | `piaso.tl.predictCellTypeByMarker` / `predictCellTypeByGDR` | Marker-based and reference-based annotation |
| SCALAR (single-cell LR) | `piaso.tl.runSCALAR` | Cell-type-resolved ligand–receptor inference for dissociated scRNA-seq |
| Marker-guided integration | `piaso.tl.stitchSpace` | Batch correction of an embedding via COSG-marker graph pruning |
| PIASOmarkerDB | `piaso.tl.queryPIASOmarkerDB` / `getMarkers` / `analyzeMarkers` | Client for the curated PIASO marker-gene database (live API) |

**Utilities & standard building blocks**

| Capability | Entry point | What it does |
|---|---|---|
| SVD embedding | `piaso.tl.runSVDLazy` / `runSVD` | Convenience wrapper around truncated SVD with INFOG-aware HVG (SVD itself is a standard method) |
| Local sub-clustering | `piaso.tl.leiden_local` | Re-cluster selected groups locally |
| Preprocessing (`piaso.pp`) | `piaso.pp.table` / `getCrossCategories` / `rotateSpatialCoordinates` | Table/cross-tab helpers and spatial-coordinate rotation |
| Plotting (`piaso.pl`) | `piaso.pl.plot_embeddings_split` / `plot_features_violin` / `plotConfusionMatrix` / LR plots | Embedding, violin, confusion-matrix, and ligand–receptor plots |

## Install (per agent)

Users work in **their own** analysis repos, so drop the right snippet into your setup. All of
these are generated from `canonical/` and live under [`dist/`](dist/).

**Claude Code** — add this repo as a plugin marketplace and install the `piaso` skill:
```bash
claude plugin marketplace add genecell/PIASO-for-agents
claude plugin install piaso@PIASO-for-agents
```

**Claude.ai (web app)** — upload the generated skill as a Skill (Pro/Max/Team/Enterprise, with
code execution enabled). Download the [`dist/claude/skills/piaso/`](dist/claude/skills/piaso)
folder, zip it, then in claude.ai go to **Settings → Capabilities → Skills → Create skill** and
upload the zip:
```bash
# from a clone of this repo:
cd dist/claude/skills && zip -r piaso-skill.zip piaso    # -> upload piaso-skill.zip in claude.ai
```
The local MCP server below is stdio-only, so it does **not** work in the web app — use the Skill
upload (or the `llms.txt` URL) on claude.ai; use MCP in Claude Code / Cursor / Codex.

**Cursor** — download the rule into your project's `.cursor/rules/`:
```bash
curl -L https://raw.githubusercontent.com/genecell/PIASO-for-agents/master/dist/cursor/.cursor/rules/piaso.mdc \
  -o .cursor/rules/piaso.mdc
```

**GitHub Copilot** — copy the instructions file into your repo:
```bash
curl -L https://raw.githubusercontent.com/genecell/PIASO-for-agents/master/dist/copilot/.github/copilot-instructions.md \
  -o .github/copilot-instructions.md
```

**OpenAI Codex** — add the `AGENTS.md` pointer below to your project's `AGENTS.md` (Codex's
primary instructions file), and/or register the MCP server (see the **MCP server** section
below — Codex is covered there).

**AGENTS.md (Aider / Zed / Codex / any AGENTS.md-aware agent)** — append the hub pointer to
your project's `AGENTS.md` (or copy [`dist/agents/AGENTS.md`](dist/agents/AGENTS.md)):
> This project uses the PIASO single-cell omics ecosystem. Agent-neutral, tested docs for
> every component (Python + R), plus the cross-component decision rules, live at
> https://github.com/genecell/PIASO-for-agents

**llms.txt (any model with web access)** — point the tool at:
```
https://raw.githubusercontent.com/genecell/PIASO-for-agents/master/dist/llms/llms.txt
```
(and `llms-full.txt` alongside it). These can also be served from `https://piaso.org/llms.txt`.

## MCP server

`piaso-mcp` serves the PIASO ecosystem docs **plus** the live PIASOmarkerDB — no
Python packages required. Tools: `search_docs`, `get_api`, `compare_implementations`,
`resolve_install`, `list_datasets`, and the live DB proxies `query_marker_db`,
`get_markers`, `list_studies`. It is a **local stdio** server (not a hosted remote endpoint),
so it works in Claude Code / Cursor / VS Code / Windsurf / Zed / Codex / Cline, but **not** in
the claude.ai web app — use the Skill upload there.

### Prerequisite (all clients): `uv`

The server runs via `uvx`, which ships with **`uv`**. This is the one thing "no packages
needed" doesn't cover — install it once:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh   # macOS / Linux
# or:  pipx install uv   |   pip install --user uv   |   brew install uv   |   winget install astral-sh.uv
```

Then confirm it's reachable: `uvx --version`. **If that says "command not found"**, uv's bin
dir isn't on your PATH — either add it, or replace `"uvx"` in the configs below with the
**absolute path** from `which uvx` (Windows: `where uvx`). First launch downloads the package
(~30 s); later launches are cached.

The MCP **config key and file location differ per client** — pick your agent below.

### Claude Code — key `mcpServers`

Easiest is the CLI (no hand-editing, and it handles the PATH issue in one line):

```bash
claude mcp add piaso --scope user -- uvx piaso-mcp
# uvx not on PATH? use its absolute path:
claude mcp add piaso --scope user -- "$(which uvx)" piaso-mcp

claude mcp get piaso        # verify → Status: ✔ Connected
```

Or edit `~/.claude.json` (user) / project `.mcp.json`:

```jsonc
{ "mcpServers": { "piaso": { "command": "uvx", "args": ["piaso-mcp"] } } }
```

### Cursor — key `mcpServers`

File: `~/.cursor/mcp.json` (global) or `.cursor/mcp.json` (per project). Same shape as Claude Code:

```jsonc
{ "mcpServers": { "piaso": { "command": "uvx", "args": ["piaso-mcp"] } } }
```

Enable it under **Settings → MCP**.

### Windsurf — key `mcpServers`

File: `~/.codeium/windsurf/mcp_config.json` (open via **Settings → Cascade → MCP Servers → Manage → raw config**):

```jsonc
{ "mcpServers": { "piaso": { "command": "uvx", "args": ["piaso-mcp"] } } }
```

### VS Code (GitHub Copilot, Agent mode) — key `servers` (note: not `mcpServers`)

Workspace file `.vscode/mcp.json`, or user `settings.json` under `"mcp"`. VS Code also wants a `type`:

```jsonc
// .vscode/mcp.json
{ "servers": { "piaso": { "type": "stdio", "command": "uvx", "args": ["piaso-mcp"] } } }
```

Or one-shot from the terminal:

```bash
code --add-mcp '{"name":"piaso","command":"uvx","args":["piaso-mcp"]}'
```

### Zed — key `context_servers` (different shape)

File: `~/.config/zed/settings.json`. Zed nests under `context_servers` and marks custom servers with `"source": "custom"`:

```jsonc
{
  "context_servers": {
    "piaso": { "source": "custom", "command": "uvx", "args": ["piaso-mcp"], "env": {} }
  }
}
```

### Codex (OpenAI Codex CLI) — TOML, table `[mcp_servers.<name>]` (not JSON!)

Codex is the odd one out: its config is **TOML**, in `~/.codex/config.toml`. Add a table:

```toml
[mcp_servers.piaso]
command = "uvx"
args = ["piaso-mcp"]
# uvx not on PATH? give the absolute path from `which uvx`:
# command = "/home/you/.local/bin/uvx"
```

Or use the CLI (handles the file for you):

```bash
codex mcp add piaso -- uvx piaso-mcp
codex mcp list        # verify it's registered
```

### Cline / Continue (VS Code extensions) — key `mcpServers`

Cline: **MCP Servers → Configure** (writes `cline_mcp_settings.json`). Continue: `~/.continue/config` (`mcpServers`). Both use the standard shape:

```jsonc
{ "mcpServers": { "piaso": { "command": "uvx", "args": ["piaso-mcp"] } } }
```

---

**After configuring, restart the client** — MCP tools are loaded at startup, so a running
session won't see the server until it's relaunched. If it doesn't connect, 99% of the time
it's the `uv`/PATH prerequisite above.

## Repository layout

```
canonical/       # the ONLY hand-written content (agent-neutral markdown + meta.yaml)
build.py         # canonical/ -> all targets (pure text transforms); --check is the CI drift guard
dist/            # ALL GENERATED — never hand-edited (claude/ agents/ cursor/ copilot/ llms/ mcp/)
mcp/             # piaso-mcp source (local stdio server; serves knowledge + public data only)
tests/           # executes every canonical code block (Python + R) against PIASO-data fixtures
.claude-plugin/  # marketplace + plugin manifest (repo root, for `claude plugin marketplace add`)
.github/         # sync-check + test CI (re-runs on component releases + nightly)
```

## Citation

Cite each component by its own paper — see [`canonical/meta.yaml`](canonical/meta.yaml).
PIASO: Wu, S.J., Dai, M. *et al.* *Nature* (2026), DOI `10.1038/s41586-025-09996-8`.

## License

BSD-3-Clause. See [`LICENSE`](LICENSE).
