#!/usr/bin/env python3
"""build.py — generate every agent target from canonical/ into dist/.

Pure text transformation. Reads canonical/meta.yaml + canonical/*.md and writes:
  dist/claude/skills/piaso/{SKILL.md, references/*.md, LICENSE.txt}
  dist/agents/AGENTS.md                    (hub root)
  dist/agents/components/<repo>/AGENTS.md  (fan-out, one per component repo)
  dist/cursor/.cursor/rules/piaso.mdc
  dist/copilot/.github/copilot-instructions.md
  dist/llms/{llms.txt, llms-full.txt}      (relative links — for the hub)
  dist/llms/piaso.org/{llms.txt, llms-full.txt}  (absolute raw-GitHub links — drop-in for piaso.org/web/public/)
  dist/mcp/                                (copied from mcp/ source + the bundled canonical snapshot)

Also generates one derived canonical file, tutorials.md (the piaso.org tutorial index from
meta.yaml), so every target and the MCP snapshot carry it.

No cleverness: flat-file targets (Copilot, one Cursor rule) get overview + a workflow digest
capped to the format's size budget (see docs/recon/format-matrix.md); multi-file targets
(Claude skill) map canonical files 1:1 into references/.

Run:  python build.py            # build all
      python build.py --check    # build to a temp dir and diff against dist/ (CI sync-check)
"""
from __future__ import annotations
import argparse, shutil, sys, tempfile, filecmp, re
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parent
CANON = ROOT / "canonical"
META = yaml.safe_load((CANON / "meta.yaml").read_text())

# Claude's description field is the trigger and is hard-capped at 1024 chars.
DESC_CAP = 1024
COPILOT_CAP = 20000  # soft; must hold the whole overview (router + decision rules) + the workflow list
RAW = "https://raw.githubusercontent.com/genecell/PIASO-for-agents/master/canonical/"


def read(p: Path) -> str:
    return p.read_text() if p.exists() else ""


def canon_file(rel: str) -> str:
    return read(CANON / rel)


def workflow_files() -> list[Path]:
    d = CANON / "workflows"
    return sorted(x for x in d.glob("*.md") if x.name.lower() != "readme.md")


def component_files() -> list[Path]:
    return sorted((CANON / "components").glob("*.md"))


def code_components() -> list[dict]:
    return [c for c in META["components"] if c["id"] != "piaso-data"]


# ---------------------------------------------------------------- derived text
def tested_against() -> str:
    parts = [f"{c.get('pypi') or c['name']} {c['version_last_tested']}" for c in META["components"]
             if c.get("version_last_tested")]
    return f"Tested against: {' · '.join(parts)} ({META['hub']['tested_on']})."


def tutorials_md() -> str:
    """The piaso.org tutorial index as markdown (generated from meta.yaml)."""
    base = META["hub"]["tutorials"].rstrip("/") + "/"
    out = ["# piaso.org tutorials — index", "",
           f"Executed, human-reviewed tutorials at {META['hub']['tutorials']} (generated API reference: "
           f"{META['hub']['api_reference']}). Each runs against the real dataset it names; the numbers and "
           "figures are what the code produced. Grouped by topic; `components` says which packages it uses.", ""]
    by_topic: dict[str, list[dict]] = {}
    for t in META.get("tutorials", []):
        by_topic.setdefault(t["topic"], []).append(t)
    for topic, items in by_topic.items():
        out.append(f"## {topic}")
        out.append("")
        for t in items:
            comps = ", ".join(t.get("components", [])) or "—"
            out.append(f"- [{t['title']}]({base}{t['slug']}/) — {t['covers']} *({comps})*")
        out.append("")
    return "\n".join(out)


# ---------------------------------------------------------------- trigger text
def build_description() -> str:
    """One-paragraph trigger from meta.triggers, hard-capped at 1024 chars."""
    t = META["triggers"]
    names = ", ".join(t["component_names"])
    tasks = "; ".join(t["task_phrasings"])
    fmts = ", ".join(t["file_formats"])
    desc = (
        "PIASO single-cell omics ecosystem (scRNA-seq, spatial; Python + R). "
        f"Use for: {tasks}. "
        f"Components (each a standalone trigger): {names}. "
        f"Fires on {fmts}, and on any component name even without 'PIASO'. "
        "Includes SCALAR-vs-LARIS, AnnData-vs-cytome, COSG-vs-cytorete and Python-vs-R choices."
    )
    if len(desc) > DESC_CAP:
        raise SystemExit(f"skill description is {len(desc)} chars > {DESC_CAP}; trim meta.yaml triggers")
    return desc


def digest(max_chars: int) -> str:
    """overview.md + a compressed workflow list, capped — for flat-file targets."""
    parts = [canon_file("overview.md").strip()]
    wf = ["\n## Workflows (see the hub for full code)"]
    for f in workflow_files():
        first = next((ln for ln in read(f).splitlines() if ln.startswith("# ")), f.stem)
        wf.append(f"- **{first.lstrip('# ').strip()}** (`workflows/{f.name}`)")
    wf.append(f"\nTutorial index: {META['hub']['tutorials']} · {tested_against()}")
    parts.append("\n".join(wf))
    out = "\n\n".join(parts)
    return out[:max_chars].rstrip() + ("\n" if len(out) <= max_chars else "\n\n*(truncated — see the hub)*\n")


# ---------------------------------------------------------------- targets
def build_claude(dist: Path) -> None:
    skill = dist / "claude" / "skills" / "piaso"
    refs = skill / "references"
    refs.mkdir(parents=True, exist_ok=True)
    front = f"---\nname: piaso\ndescription: {build_description()}\nlicense: BSD-3-Clause\n---\n\n"
    body = canon_file("overview.md").strip() + "\n\n## References\n"
    body += "\nComponent references (each self-sufficient):\n"
    for f in component_files():
        body += f"- `references/components/{f.name}`\n"
    body += "\nWorkflow references (cross-component tasks):\n"
    for f in workflow_files():
        body += f"- `references/workflows/{f.name}`\n"
    body += ("\nAlso: `references/gotchas.md` (layer contracts, deprecated names), `references/data.md` "
             "(fixtures, registry), `references/tutorials.md` (the piaso.org tutorial index — point the "
             "user at the executed tutorial for their platform).\n\n" + tested_against() + "\n")
    (skill / "SKILL.md").write_text(front + body)
    # references map 1:1
    (refs / "components").mkdir(exist_ok=True)
    (refs / "workflows").mkdir(exist_ok=True)
    for f in component_files():
        shutil.copy(f, refs / "components" / f.name)
    for f in workflow_files():
        shutil.copy(f, refs / "workflows" / f.name)
    for extra in ("gotchas.md", "data.md"):
        if (CANON / extra).exists():
            shutil.copy(CANON / extra, refs / extra)
    (refs / "tutorials.md").write_text(tutorials_md())
    shutil.copy(ROOT / "LICENSE", skill / "LICENSE.txt")


def _agents_body(scope: str) -> str:
    return (
        f"# AGENTS.md — {scope}\n\n"
        "This repository is part of the **PIASO single-cell omics ecosystem**. Full, "
        "cross-component, agent-neutral documentation (with executed code blocks for every "
        "component in Python and R) lives in the hub:\n"
        "**https://github.com/genecell/PIASO-for-agents** — tutorials: "
        f"{META['hub']['tutorials']} — API reference: {META['hub']['api_reference']}\n\n"
        "## Ecosystem at a glance\n"
        + "\n".join(
            f"- **{c['name']}** (`{c.get('pypi') or c.get('import', '')}`, {'/'.join(c['language'])}): "
            f"{c.get('role', '')} — install `{c.get('install', 'see hub')}`"
            for c in code_components()
        )
        + "\n\n## Cross-component decision rules\n"
        + "\n".join(f"- **{r['id']}**: {r['rule'].strip()}" for r in META["decision_rules"])
        + f"\n\n{tested_against()} For full API, workflows, and citations, read the hub.\n"
        + _maintainer_footer()
    )


def _maintainer_footer() -> str:
    h = META["hub"]
    lab, url, aff = h.get("lab", ""), h.get("lab_url", ""), h.get("affiliation", "")
    return f"\n---\nMaintained by **[{lab}]({url})** ({aff}).\n" if lab else ""


def build_agents(dist: Path) -> None:
    root = dist / "agents"
    root.mkdir(parents=True, exist_ok=True)
    (root / "AGENTS.md").write_text(_agents_body("PIASO ecosystem hub"))
    for c in META["components"]:
        repo = c["repo"].rstrip("/").split("/")[-1]
        d = root / "components" / repo
        d.mkdir(parents=True, exist_ok=True)
        (d / "AGENTS.md").write_text(_agents_body(f"genecell/{repo}"))


def build_cursor(dist: Path) -> None:
    d = dist / "cursor" / ".cursor" / "rules"
    d.mkdir(parents=True, exist_ok=True)
    front = (
        "---\n"
        f"description: {build_description()}\n"
        "globs: [\"*.h5ad\", \"*.rds\", \"*.cytome\", \"**/*.py\", \"**/*.R\"]\n"
        "alwaysApply: false\n"
        "---\n\n"
    )
    (d / "piaso.mdc").write_text(front + digest(20000))


def build_copilot(dist: Path) -> None:
    d = dist / "copilot" / ".github"
    d.mkdir(parents=True, exist_ok=True)
    (d / "copilot-instructions.md").write_text(
        "# PIASO ecosystem — Copilot instructions\n\n" + digest(COPILOT_CAP)
    )


def _llms_index(link) -> str:
    hub = META["hub"]
    idx = [f"# {hub['name']}", "", f"> {hub['description'].strip()}", "",
           f"{tested_against()} Tutorials: {hub['tutorials']} · API reference: {hub['api_reference']}", ""]
    idx.append("## Components")
    doc_stems = {f.stem for f in component_files()}
    for c in code_components():
        # components without their own file are documented inside their counterpart's file
        stem = c["id"] if c["id"] in doc_stems else c.get("counterpart", c["id"])
        note = "" if c["id"] in doc_stems else f" (documented with {stem})"
        idx.append(f"- [{c['name']}]({link(f'components/{stem}.md')}): {'/'.join(c['language'])}, "
                   f"{c.get('role', '')}; install `{c.get('install', 'see hub')}`{note}")
    idx.append("\n## Workflows")
    for f in workflow_files():
        first = next((ln for ln in read(f).splitlines() if ln.startswith("# ")), f.stem)
        idx.append(f"- [{first.lstrip('# ').strip()}]({link(f'workflows/{f.name}')})")
    idx.append("\n## Tutorials (piaso.org, executed)")
    base = hub["tutorials"].rstrip("/") + "/"
    for t in META.get("tutorials", []):
        idx.append(f"- [{t['title']}]({base}{t['slug']}/): {t['covers']}")
    idx.append("\n## Optional")
    idx.append(f"- [gotchas]({link('gotchas.md')})")
    idx.append(f"- [data fixtures and registry]({link('data.md')})")
    idx.append(f"- [tutorial index]({link('tutorials.md')})")
    for s in hub.get("see_also", []):
        idx.append(f"- [{s['name']}]({s['url']}) — {s['role']}")
    idx.append(_maintainer_footer())
    return "\n".join(idx) + "\n"


def _llms_full() -> str:
    full = [f"# {META['hub']['name']} — full\n\n{tested_against()}\n"]
    for name in ["overview.md"] + [f"components/{f.name}" for f in component_files()] \
            + [f"workflows/{f.name}" for f in workflow_files()] + ["gotchas.md", "data.md"]:
        p = CANON / name
        if p.exists():
            full.append(f"\n\n<!-- ===== {name} ===== -->\n\n" + p.read_text())
    full.append("\n\n<!-- ===== tutorials.md ===== -->\n\n" + tutorials_md())
    return "".join(full)


def build_llms(dist: Path) -> None:
    d = dist / "llms"
    d.mkdir(parents=True, exist_ok=True)
    (d / "llms.txt").write_text(_llms_index(lambda rel: rel))
    (d / "llms-full.txt").write_text(_llms_full())
    # the piaso.org copy: every canonical link absolute (raw GitHub), same content otherwise
    site = d / "piaso.org"
    site.mkdir(exist_ok=True)
    (site / "llms.txt").write_text(_llms_index(lambda rel: RAW + rel))
    (site / "llms-full.txt").write_text(_llms_full())


def build_mcp(dist: Path) -> None:
    src = ROOT / "mcp"
    if not src.exists():
        return
    dest = dist / "mcp"
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(src, dest, ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.egg-info", "dist", "build"))
    # populate the server's bundled canonical snapshot (generated, never hand-edited)
    data = dest / "piaso_mcp" / "data"
    if data.exists():
        shutil.rmtree(data)
    (data / "components").mkdir(parents=True, exist_ok=True)
    (data / "workflows").mkdir(parents=True, exist_ok=True)
    shutil.copy(CANON / "meta.yaml", data / "meta.yaml")
    for name in ("overview.md", "gotchas.md", "data.md"):
        if (CANON / name).exists():
            shutil.copy(CANON / name, data / name)
    (data / "tutorials.md").write_text(tutorials_md())
    for f in component_files():
        shutil.copy(f, data / "components" / f.name)
    for f in workflow_files():
        shutil.copy(f, data / "workflows" / f.name)


TARGETS = {
    "claude": build_claude, "agents": build_agents, "cursor": build_cursor,
    "copilot": build_copilot, "llms": build_llms, "mcp": build_mcp,
}


def build_all(dist: Path) -> None:
    if dist.exists():
        shutil.rmtree(dist)
    dist.mkdir(parents=True)
    for name, fn in TARGETS.items():
        fn(dist)
    print("built:", ", ".join(TARGETS), f"| skill description {len(build_description())}/{DESC_CAP} chars")


def check() -> int:
    """Build to a temp dir and diff against dist/ — the CI sync-check."""
    with tempfile.TemporaryDirectory() as tmp:
        tmpdist = Path(tmp) / "dist"
        build_all(tmpdist)
        diff = _dircmp_diff(ROOT / "dist", tmpdist)
        if diff:
            print("DRIFT: dist/ is out of sync with canonical/. Run `python build.py`.")
            for d in diff:
                print("  ", d)
            return 1
    print("sync-check OK: dist/ matches canonical/")
    return 0


def _dircmp_diff(a: Path, b: Path) -> list[str]:
    out: list[str] = []

    def walk(x: Path, y: Path, rel: str = ""):
        cmp = filecmp.dircmp(x, y, ignore=["__pycache__", "dist", "build"])
        for n in cmp.left_only:
            out.append(f"only in dist/: {rel}{n}")
        for n in cmp.right_only:
            out.append(f"missing from dist/: {rel}{n}")
        for n in cmp.diff_files:
            out.append(f"differs: {rel}{n}")
        for sub in cmp.common_dirs:
            walk(x / sub, y / sub, rel + sub + "/")

    if not a.exists():
        return ["dist/ does not exist"]
    walk(a, b)
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="CI sync-check against dist/")
    args = ap.parse_args()
    if args.check:
        return check()
    build_all(ROOT / "dist")
    return 0


if __name__ == "__main__":
    sys.exit(main())
