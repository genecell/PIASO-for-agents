#!/usr/bin/env python3
"""build.py — generate every agent target from canonical/ into dist/.

Pure text transformation. Reads canonical/meta.yaml + canonical/*.md and writes:
  dist/claude/skills/piaso/{SKILL.md, references/*.md, LICENSE.txt}
  dist/agents/AGENTS.md                    (hub root)
  dist/agents/components/<repo>/AGENTS.md  (fan-out, one per component repo)
  dist/cursor/.cursor/rules/piaso.mdc
  dist/copilot/.github/copilot-instructions.md
  dist/llms/{llms.txt, llms-full.txt}
  dist/mcp/                                (copied from mcp/ source; server is generated-once)

No cleverness: flat-file targets (Copilot, one Cursor rule) get overview + a workflow digest
capped to the format's size budget (see docs/recon/format-matrix.md); multi-file targets
(Claude skill) map canonical files 1:1 into references/.

Run:  python build.py            # build all
      python build.py --check    # build to a temp dir and diff against dist/ (CI sync-check)
"""
from __future__ import annotations
import argparse, shutil, sys, tempfile, filecmp, os
from pathlib import Path
import yaml

ROOT = Path(__file__).resolve().parent
CANON = ROOT / "canonical"
META = yaml.safe_load((CANON / "meta.yaml").read_text())

# Claude's description field is the trigger and is hard-capped at 1024 chars.
DESC_CAP = 1024
COPILOT_CAP = 6000  # ~2 pages, soft; keep the flat file tight


def read(p: Path) -> str:
    return p.read_text() if p.exists() else ""


def canon_file(rel: str) -> str:
    return read(CANON / rel)


def workflow_files() -> list[Path]:
    d = CANON / "workflows"
    return sorted(x for x in d.glob("*.md") if x.name.lower() != "readme.md")


def component_files() -> list[Path]:
    return sorted((CANON / "components").glob("*.md"))


# ---------------------------------------------------------------- trigger text
def build_description() -> str:
    """One-paragraph trigger from meta.triggers, capped at 1024 chars."""
    t = META["triggers"]
    names = ", ".join(t["component_names"])
    tasks = "; ".join(t["task_phrasings"])
    fmts = ", ".join(t["file_formats"])
    desc = (
        "PIASO single-cell omics ecosystem for scRNA-seq / spatial analysis. "
        f"Use for: {tasks}. "
        f"Components (each a standalone trigger): {names}. "
        "Covers Python (scanpy/AnnData) and R (Seurat/SCE, via COSGR). "
        f"Fires on {fmts}, and on any of the component names even without 'PIASO'. "
        "Includes the single-cell (SCALAR) vs spatial (LARIS) ligand-receptor choice."
    )
    if len(desc) > DESC_CAP:
        desc = desc[: DESC_CAP - 1].rstrip() + "…"
    return desc


def digest(max_chars: int) -> str:
    """overview.md + a compressed workflow list, capped — for flat-file targets."""
    parts = [canon_file("overview.md").strip()]
    wf = ["\n## Workflows (see the hub for full code)"]
    for f in workflow_files():
        first = next((ln for ln in read(f).splitlines() if ln.startswith("# ")), f.stem)
        wf.append(f"- **{first.lstrip('# ').strip()}** (`workflows/{f.name}`)")
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
    body += "\nAlso: `references/gotchas.md`, `references/data.md`.\n"
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
    shutil.copy(ROOT / "LICENSE", skill / "LICENSE.txt")


def _agents_body(scope: str) -> str:
    return (
        f"# AGENTS.md — {scope}\n\n"
        "This repository is part of the **PIASO single-cell omics ecosystem**. Full, "
        "cross-component, agent-neutral documentation (with runnable, tested code blocks "
        "for every component in Python and R) lives in the hub:\n"
        "**https://github.com/genecell/PIASO-for-agents**\n\n"
        "## Ecosystem at a glance\n"
        + "\n".join(
            f"- **{c['id']}** (`{c.get('pypi') or c.get('import','')}`, {'/'.join(c['language'])}): "
            f"install `{c.get('install','see hub')}`"
            for c in META["components"] if c["id"] != "piaso-data"
        )
        + "\n\n## Cross-component decision rules\n"
        + "\n".join(f"- **{r['id']}**: {r['rule'].strip()}" for r in META["decision_rules"])
        + "\n\nFor full API, workflows, and citations, read the hub.\n"
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
        "globs: [\"*.h5ad\", \"*.rds\", \"**/*.py\", \"**/*.R\"]\n"
        "alwaysApply: false\n"
        "---\n\n"
    )
    (d / "piaso.mdc").write_text(front + digest(9000))


def build_copilot(dist: Path) -> None:
    d = dist / "copilot" / ".github"
    d.mkdir(parents=True, exist_ok=True)
    (d / "copilot-instructions.md").write_text(
        "# PIASO ecosystem — Copilot instructions\n\n" + digest(COPILOT_CAP)
    )


def build_llms(dist: Path) -> None:
    d = dist / "llms"
    d.mkdir(parents=True, exist_ok=True)
    hub = META["hub"]
    idx = [f"# {hub['name']}", "", f"> {hub['description'].strip()}", ""]
    idx.append("## Components")
    doc_stems = {f.stem for f in component_files()}
    for c in META["components"]:
        if c["id"] == "piaso-data":
            continue
        # resolve to an existing doc file; components without their own file
        # (e.g. cosgr) are documented inside their counterpart's file (cosg.md).
        stem = c["id"] if c["id"] in doc_stems else c.get("counterpart", c["id"])
        note = "" if c["id"] in doc_stems else f" (documented with {stem})"
        idx.append(f"- [{c['id']}](components/{stem}.md): {'/'.join(c['language'])}, "
                   f"install `{c.get('install','see hub')}`{note}")
    idx.append("\n## Workflows")
    for f in workflow_files():
        first = next((ln for ln in read(f).splitlines() if ln.startswith("# ")), f.stem)
        idx.append(f"- [{first.lstrip('# ').strip()}](workflows/{f.name})")
    idx.append("\n## Optional")
    idx.append("- [gotchas](gotchas.md)")
    idx.append("- [data fixtures](data.md)")
    idx.append(_maintainer_footer())
    (d / "llms.txt").write_text("\n".join(idx) + "\n")
    # llms-full.txt = everything concatenated
    full = [f"# {hub['name']} — full\n"]
    for name in ["overview.md"] + [f"components/{f.name}" for f in component_files()] \
            + [f"workflows/{f.name}" for f in workflow_files()] + ["gotchas.md", "data.md"]:
        p = CANON / name
        if p.exists():
            full.append(f"\n\n<!-- ===== {name} ===== -->\n\n" + p.read_text())
    (d / "llms-full.txt").write_text("".join(full))


def build_mcp(dist: Path) -> None:
    src = ROOT / "mcp"
    if not src.exists():
        return
    dest = dist / "mcp"
    if dest.exists():
        shutil.rmtree(dest)
    shutil.copytree(src, dest, ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.egg-info"))
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
    print("built:", ", ".join(TARGETS))


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
        cmp = filecmp.dircmp(x, y)
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
