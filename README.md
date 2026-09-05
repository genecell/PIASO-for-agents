# PIASO-for-agents

**Make the [PIASO](https://piaso.org) single-cell omics ecosystem first-class for any
coding agent — Claude Code, Cursor, Copilot, Codex, Windsurf, Cline, Aider — from one
canonical, agent-neutral knowledge pack.**

Maintained by **[The Fishell Laboratory](https://fishelllab.hms.harvard.edu)** (Harvard
Medical School / Broad Institute). Every agent-specific format (Claude skill, Cursor rules, `AGENTS.md`,
`llms.txt`, MCP server) is a **generated artifact** built from `canonical/` — never a
hand-maintained copy. A CI drift check (`python build.py --check`) fails the build if any
`dist/` artifact is out of sync with `canonical/`, and the code-block test suite runs every
canonical block against the **pinned component versions** on every push, nightly, and on
component releases, so the guidance cannot silently rot.

**Hub v0.2.0 · piaso-mcp 0.1.0 — tested against piaso-tools 1.2.3 · cosg 1.2.0 · cytome 0.3.1 ·
laris 0.13.0 · emergene 1.0.2 · cytorete 0.1.1 · COSGR 1.0.0 · cytome (R) 0.1.0 (2026-09-04).**

## The ecosystem

Independently-installable packages under [github.com/genecell](https://github.com/genecell), in
four layers. Dependencies run one way (`cytorete → piaso-tools → cosg + cytome`; `laris → cosg`),
and `pip install piaso-tools` already brings COSG and cytome.

| Layer | Component | Package | Language | Role |
|---|---|---|---|---|
| Analysis | [PIASO](https://github.com/genecell/PIASO) | `piaso-tools` | Python + Rust | Self-contained pipeline — reading 10x data, QC, doublets, **INFOG**, SVD / **GDR**, Leiden / UMAP, **PIASOscore**, annotation, **SCALAR**, PIASOmarkerDB client, plotting, `piaso.data`. **No scanpy required.** |
| Storage | [cytome](https://github.com/genecell/cytome) | `cytome` | Python | Single-file SQLite `.cytome`: matrices, SQL-queryable cell/gene tables, embeddings, graphs, fragments, tissue images, provenance — what every component streams from |
| | [cytome (R)](https://github.com/genecell/cytome-r) | `cytome` (r-universe) | R | Read / write / stream the same file into Seurat or SingleCellExperiment, no Python |
| Methods | [COSG](https://github.com/genecell/COSG) | `cosg` | Python | Marker genes by cosine specificity — analytic p-values, GPU, batch-aware, streams from cytome |
| | [COSGR](https://github.com/genecell/COSGR) | `COSG` (r-universe / conda-forge) | R | COSG for Seurat / SingleCellExperiment |
| | [LARIS](https://github.com/genecell/LARIS) | `laris` | Python | Ligand–receptor interaction in **spatial** transcriptomics; exact p-values; cross-condition comparison |
| | [Emergene](https://github.com/genecell/Emergene) | `emergene` | Python | Individual-cell differential expression across conditions |
| | [cytorete](https://github.com/genecell/cytorete) | `cytorete` | Python | Cell-type-resolved gene regulatory networks (regulons) on the PIASO stack |
| Data | [PIASO-data](https://github.com/genecell/PIASO-data) | — | data | Tutorial datasets (Zenodo, incl. five `.cytome` atlases) + genome references; registry read by `piaso.data` |

Each component is **independently installable** — a COSG-only, cytome-only or LARIS-only user is
a first-class citizen, and every `canonical/components/*.md` assumes nothing else is installed.
The hub's unique value is documenting how the components **compose**, and the cross-component
choices no single repo can make: **SCALAR vs LARIS** (dissociated vs spatial ligand–receptor —
same CellChatDB either way), **AnnData vs `.cytome`** (in memory vs streamed — same function
calls), **COSG vs cytorete** (marker genes vs the TFs that drive them), **Python vs R** (COSG →
COSGR, cytome → cytome (R); everything else via a `.cytome` handoff), and **which annotation
route** (marker sets, reference projection, joint embedding, or a gene list against PIASOmarkerDB).

### Inside `piaso-tools`

Full reference: [`canonical/components/piaso.md`](canonical/components/piaso.md). Every function
takes `data=` as an AnnData, an open cytome Dataset or a `.cytome` path.

**Methods introduced by PIASO**

| Capability | Entry point | What it does |
|---|---|---|
| INFOG normalization | `piaso.tl.infog` | Information-content normalization of raw UMI counts + informative-gene selection |
| GDR (marker-gene-guided DR) | `piaso.tl.runGDR` / `runGDRParallel` / `projectGDR` | Embedding whose axes are per-group COSG-marker scores; integrates batches by identity; frozen reference spaces |
| Gene-set scoring (PIASOscore) | `piaso.tl.score` | Expression-matched-control scoring with per-cell p-values; whole pathway databases in one Rust matmul |
| Cell-type prediction | `piaso.tl.predictCellTypeByMarker` / `predictCellTypeByGDR` | Marker-set and reference-based annotation |
| SCALAR (single-cell LR) | `piaso.tl.specificity_matrix` + `runSCALAR` | Cell-type-resolved ligand–receptor inference for dissociated data, CellChatDB via `piaso.data.load_lr_database` |
| Marker-guided integration | `piaso.tl.stitchSpace` | Batch correction of an embedding via COSG-marker graph pruning |
| PIASOmarkerDB | `piaso.tl.getMarkers` / `analyzeMarkers` | Client for the curated marker database (36 studies, live API) |
| Motif scanning | `piaso.pp.scan_motifs` + `piaso.data` motif/genome loaders | The Rust PWM engine cytorete builds on |

**Pipeline building blocks (scanpy-free)**

| Capability | Entry point |
|---|---|
| Read 10x / Cell Ranger | `piaso.pp.read_10x_h5`, `read_10x`, `importCellRanger` (→ cytome) |
| QC, doublets, filtering | `piaso.pp.calculateCellMetrics`, `scrublet`, `filter_cells`, `calculateGroupMetrics` |
| Embedding, graph, clusters, UMAP | `piaso.tl.runSVD`, `neighbors`, `leiden`, `umap`, `leiden_local`, `runHarmony` |
| Datasets, genomes, motif DBs, CellChatDB | `piaso.data.load_dataset`, `fetch_genome`, `fetch_2bit`, `fetch_jaspar`, `load_lr_database` |
| Plotting | `piaso.pl.embedding`, `dotplot`, `violin`, `scatter`, `sankey`, `stackedBarplot`, `plot_embeddings_split` (+ tissue-image overlays on cytomes), `piaso.settings.set_figure_params` |

## What an agent gets

- `canonical/overview.md` — the router: task → component table and the seven decision rules.
- `canonical/components/` — self-sufficient references for PIASO, COSG (+ COSGR), cytome (+ R),
  LARIS, Emergene, cytorete, with executed code blocks and the data-object contract of every call.
- `canonical/workflows/` — end-to-end scRNA-seq (scanpy-free), streaming on a `.cytome`,
  marker-based annotation + reference projection, PIASOmarkerDB annotation, ligand–receptor
  (SCALAR and LARIS), spatial transcriptomics, gene regulatory networks.
- `canonical/gotchas.md` (layer contracts, deprecated names, the `as_dict` tuple, species-cased
  prefixes), `canonical/data.md` (registry, fixtures), and the **piaso.org tutorial index**
  (generated into every target) so the agent can point the user at the executed tutorial for
  their platform.

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
https://piaso.org/llms.txt          # and https://piaso.org/llms-full.txt
```
These are the hub's `dist/llms/piaso.org/` files (absolute links); the relative-link versions are
at `dist/llms/`.

## MCP server

`piaso-mcp` serves the PIASO ecosystem docs, the **piaso.org tutorial index**, the **PIASO-data
registry** and the **live PIASOmarkerDB** — no Python packages required. Tools: `search_docs`,
`get_api`, `compare_implementations`, `resolve_install`, `list_tutorials`, `version_matrix`,
`check_versions` (PyPI vs tested versions), `list_datasets` / `get_dataset` (live registry), and
the live DB proxies `query_marker_db`, `get_markers`, `list_studies`. It is a **local stdio**
server (not a hosted remote endpoint), so it works in Claude Code / Cursor / VS Code / Windsurf /
Zed / Codex / Cline, but **not** in the claude.ai web app — use the Skill upload there.

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
canonical/       # the ONLY hand-written content (agent-neutral markdown + meta.yaml, incl. the tutorial index)
build.py         # canonical/ -> all targets (pure text transforms); --check is the CI drift guard
dist/            # ALL GENERATED — never hand-edited (claude/ agents/ cursor/ copilot/ llms/ mcp/)
mcp/             # piaso-mcp source (local stdio server; serves knowledge + public data only)
tests/           # executes every canonical code block (Python + R) on the fixtures; heavy spatial/regulon runs nightly
.claude-plugin/  # marketplace + plugin manifest (repo root, for `claude plugin marketplace add`)
.github/         # sync-check + test CI (re-runs on component releases + nightly) + PyPI / MCP-registry publish
```

Related tooling (independent projects, listed on piaso.org's *Agents and project tooling* page):
[stato](https://stato.hiniki.com) — structured expertise management for long computational
projects; [PlanDrop](https://plandrop.ai) — plan-review-execute for Claude Code on remote machines.

## Changelog

See [`CHANGELOG.md`](CHANGELOG.md) — hub content releases and `piaso-mcp` versions, with the component versions each was tested against.

## Citation

Cite each component by its own paper — see [`canonical/meta.yaml`](canonical/meta.yaml).
PIASO: Wu, S.J., Dai, M. *et al.* *Nature* (2026), DOI `10.1038/s41586-025-09996-8`. COSG /
COSGR: Dai M, Pei X, Wang X-J, *Briefings in Bioinformatics* 23(2):bbab579 (2022). LARIS: Dai M,
Török T, Sun D, et al., bioRxiv (2025), DOI `10.1101/2025.11.26.690796`. cytome and cytorete have
no paper yet — cite the repositories.

## Maintainers

Developed and maintained by **[The Fishell Laboratory](https://fishelllab.hms.harvard.edu)**
(Harvard Medical School / Broad Institute).
Contact: Min Dai — dai@broadinstitute.org.

## License

BSD-3-Clause. See [`LICENSE`](LICENSE).
