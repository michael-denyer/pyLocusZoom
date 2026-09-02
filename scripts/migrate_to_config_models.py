"""Rewrite plot calls from flat keywords to the config models.

Rewrites every ``plot`` and ``plot_stacked`` call whose keywords include
``chrom`` so that each option is passed inside the model that declares it,
and every Manhattan and Miami call so that the column names and chromosome
order are passed as a ``GenomeWideConfig``. Files are edited in place and
the needed names are added to an existing ``from pylocuszoom import`` or as
a new import after the last one; run ``ruff check --fix`` and
``ruff format`` afterwards.

Usage:
    uv run --with libcst python scripts/migrate_to_config_models.py PATH...

A ``.ipynb`` path is rewritten cell by cell without re-executing it, and a
``.md`` path has each ``python`` fenced block rewritten in place.
"""

import re
import subprocess
import sys
from pathlib import Path

import libcst as cst

GROUPS = {
    "columns": ("ColumnConfig", ("pos_col", "p_col", "rs_col")),
    "display": (
        "DisplayConfig",
        ("snp_labels", "label_top_n", "show_recombination", "auto_genes", "figsize"),
    ),
    "ld": ("LDConfig", ("lead_pos", "ld_reference_file", "ld_col")),
    "panels": (
        "PanelInputs",
        (
            "genes_df",
            "exons_df",
            "recomb_df",
            "eqtl_df",
            "eqtl_gene",
            "eqtl_threshold",
            "finemapping_df",
            "finemapping_cs_col",
            "ld_heatmap_df",
            "ld_heatmap_snp_ids",
            "ld_heatmap_height",
            "ld_heatmap_metric",
        ),
    ),
}
GENOMEWIDE = ("chrom_col", "pos_col", "p_col", "custom_chrom_order")
GENOMEWIDE_METHODS = {
    "plot_manhattan",
    "plot_qq",
    "plot_manhattan_stacked",
    "plot_manhattan_qq",
    "plot_manhattan_qq_stacked",
    "plot_miami",
}


def _keywords(call: cst.Call) -> set:
    return {a.keyword.value for a in call.args if a.keyword is not None}


def _nested(name: str, model: str, args: list) -> cst.Arg:
    inner = [a.with_changes(comma=cst.MaybeSentinel.DEFAULT) for a in args]
    return cst.Arg(
        keyword=cst.Name(name),
        value=cst.Call(func=cst.Name(model), args=inner),
        equal=cst.AssignEqual(cst.SimpleWhitespace(""), cst.SimpleWhitespace("")),
    )


class Rewriter(cst.CSTTransformer):
    def __init__(self) -> None:
        self.needed: set = set()

    def leave_Call(self, original: cst.Call, updated: cst.Call) -> cst.Call:
        func = updated.func
        if not isinstance(func, cst.Attribute):
            return updated
        method = func.attr.value
        if method in ("plot", "plot_stacked") and "chrom" in _keywords(updated):
            return self._regroup(updated, GROUPS)
        if method in GENOMEWIDE_METHODS and _keywords(updated) & set(GENOMEWIDE):
            return self._regroup(updated, {"config": ("GenomeWideConfig", GENOMEWIDE)})
        return updated

    def _regroup(self, call: cst.Call, groups: dict) -> cst.Call:
        owner = {k: name for name, (_, keys) in groups.items() for k in keys}
        kept, moved = [], {name: [] for name in groups}
        for arg in call.args:
            key = arg.keyword.value if arg.keyword is not None else None
            if key in owner:
                moved[owner[key]].append(arg)
            else:
                kept.append(arg)
        if not any(moved.values()):
            return call
        for name, (model, _) in groups.items():
            if moved[name]:
                self.needed.add(model)
                kept.append(_nested(name, model, moved[name]))
        kept = [a.with_changes(comma=cst.MaybeSentinel.DEFAULT) for a in kept]
        return call.with_changes(args=kept)


def _add_import(module: cst.Module, names: set) -> cst.Module:
    class Importer(cst.CSTTransformer):
        done = False

        def leave_ImportFrom(self, original, updated):
            if self.done or updated.module is None:
                return updated
            if cst.Module([]).code_for_node(updated.module) != "pylocuszoom":
                return updated
            if isinstance(updated.names, cst.ImportStar):
                return updated
            present = {n.name.value for n in updated.names}
            extra = [cst.ImportAlias(cst.Name(n)) for n in sorted(names - present)]
            if not extra:
                self.done = True
                return updated
            aliases = [
                a.with_changes(comma=cst.MaybeSentinel.DEFAULT) for a in updated.names
            ]
            self.done = True
            return updated.with_changes(names=aliases + extra)

    importer = Importer()
    module = module.visit(importer)
    if importer.done:
        return module
    line = cst.parse_statement(f"from pylocuszoom import {', '.join(sorted(names))}\n")
    body = list(module.body)
    last = max(
        (
            i
            for i, s in enumerate(body)
            if isinstance(s, cst.SimpleStatementLine)
            and any(isinstance(x, (cst.Import, cst.ImportFrom)) for x in s.body)
        ),
        default=-1,
    )
    body.insert(last + 1, line)
    return module.with_changes(body=body)


def rewrite_source(source: str) -> tuple:
    module = cst.parse_module(source)
    rewriter = Rewriter()
    module = module.visit(rewriter)
    return module, rewriter.needed


def rewrite_py(path: Path) -> None:
    module, needed = rewrite_source(path.read_text())
    if needed:
        module = _add_import(module, needed)
        path.write_text(module.code)
        print(f"rewrote {path}: {', '.join(sorted(needed))}")


def rewrite_notebook(path: Path) -> None:
    import nbformat

    nb = nbformat.read(path, as_version=4)
    needed: set = set()
    for cell in nb.cells:
        if cell.cell_type != "code":
            continue
        try:
            module, cell_needed = rewrite_source(cell.source)
        except cst.ParserSyntaxError:
            continue
        if cell_needed:
            cell.source = module.code.rstrip("\n")
            needed |= cell_needed
    if not needed:
        return
    for cell in nb.cells:
        if cell.cell_type == "code" and "from pylocuszoom import" in cell.source:
            cell.source = _add_import(
                cst.parse_module(cell.source), needed
            ).code.rstrip("\n")
            break
    nbformat.write(nb, path)
    print(f"rewrote {path}: {', '.join(sorted(needed))}")


def _format(source: str) -> str:
    """Sort the imports and format a rewritten block the way ruff would a file."""
    sorted_imports = subprocess.run(
        [
            "ruff",
            "check",
            "--fix",
            "--select",
            "I",
            "--stdin-filename",
            "block.py",
            "-",
        ],
        input=source,
        capture_output=True,
        text=True,
        check=False,
    ).stdout
    return subprocess.run(
        ["ruff", "format", "--stdin-filename", "block.py", "-"],
        input=sorted_imports or source,
        capture_output=True,
        text=True,
        check=True,
    ).stdout


def rewrite_markdown(path: Path) -> None:
    text = path.read_text()
    needed: set = set()

    def rewrite_block(match: re.Match) -> str:
        source = match.group(1)
        try:
            module, block_needed = rewrite_source(source)
        except cst.ParserSyntaxError:
            return match.group(0)
        if not block_needed:
            return match.group(0)
        needed.update(block_needed)
        return "```python\n" + _format(_add_import(module, block_needed).code) + "```"

    rewritten = re.sub(r"```python\n(.*?)```", rewrite_block, text, flags=re.S)
    if needed:
        path.write_text(rewritten)
        print(f"rewrote {path}: {', '.join(sorted(needed))}")


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        path = Path(arg)
        if path.suffix == ".ipynb":
            rewrite_notebook(path)
        elif path.suffix == ".md":
            rewrite_markdown(path)
        else:
            rewrite_py(path)
