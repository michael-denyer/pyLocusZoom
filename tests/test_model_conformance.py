"""The Python agrees with the formal models in ``specs/`` where the two can be compared.

Two links are checked. Every ``file.py::symbol`` a model cites must exist, and
no model may cite a line number. And the Python functions must return what the
Lean models pin for the models' own concrete ``#guard`` vectors.

The vectors are an explicit table, not parsed out of the Lean: each model
encodes its values differently (doubled edges, ``none`` for an exception,
constructor syntax for errors), so a parser would be a second transcription.
Each row quotes its ``#guard`` verbatim instead, and
``test_vector_table_matches_the_models`` fails when a quoted guard is no longer
in the model or a guard on a transcribed function has no row.

``LiftWindow.lean`` pins no concrete vector as a ``#guard``; its properties are
tested on the Python by the hypothesis tests in ``tests/test_liftover.py``.
"""

import ast
import re
import warnings
from pathlib import Path
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest
import requests

from pylocuszoom import _http
from pylocuszoom.backends._coerce import split_pixels
from pylocuszoom.backends.composition import cell_edges, heatmap_highlight_rects
from pylocuszoom.backends.plotly_layout import _Panel, secondary_axis_key
from pylocuszoom.config import GenomeWideStyle
from pylocuszoom.exceptions import EnsemblAPIError
from pylocuszoom.gene_track import assign_gene_positions
from pylocuszoom.manhattan import GenomeLayout
from pylocuszoom.miami_plotter import MiamiPlotter
from tests.figure_probes import PROBES

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src" / "pylocuszoom"
SPEC_FILES = sorted(
    path
    for pattern in ("lean/*.lean", "tla/*.tla", "tla/*.matrix")
    for path in (ROOT / "specs").glob(pattern)
)

LINE_CITATION = re.compile(
    r"\.py`?:\d|\.py`? l\.|\bl\.\s?\d+|(?<![\w)\]]):\d+|\blines? \d+"
)
SYMBOL_CITATION = re.compile(r"([\w/]+\.py)::([A-Za-z_][\w.]*)")


def _source_file(name):
    """The file under ``src/pylocuszoom`` a citation's path names, or None."""
    name = name.removeprefix("src/pylocuszoom/")
    hits = [path for path in SRC.rglob("*.py") if path.as_posix().endswith("/" + name)]
    return min(hits, key=lambda path: len(path.parts), default=None)


def _defined_names(path):
    """Functions, classes, methods and assigned names of one module, dotted."""
    found = set()

    def walk(node, prefix):
        for child in ast.iter_child_nodes(node):
            if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef)):
                found.add(prefix + child.name)
            elif isinstance(child, ast.ClassDef):
                found.add(prefix + child.name)
                walk(child, prefix + child.name + ".")
            elif isinstance(child, ast.Assign):
                found.update(
                    prefix + t.id for t in child.targets if isinstance(t, ast.Name)
                )
            elif isinstance(child, ast.AnnAssign) and isinstance(
                child.target, ast.Name
            ):
                found.add(prefix + child.target.id)

    walk(ast.parse(path.read_text()), "")
    return found


@pytest.mark.parametrize("spec", SPEC_FILES, ids=lambda path: path.name)
def test_models_cite_python_by_symbol_that_exists(spec):
    """A renamed or deleted function fails here until its model follows."""
    problems = []
    for number, line in enumerate(spec.read_text().splitlines(), 1):
        if LINE_CITATION.search(line):
            problems.append(f"{spec.name}:{number} cites a line number: {line.strip()}")
        for file_name, symbol in SYMBOL_CITATION.findall(line):
            source = _source_file(file_name)
            if source is None or symbol.rstrip(".") not in _defined_names(source):
                problems.append(f"{spec.name}:{number} cites {file_name}::{symbol}")

    assert problems == []


def test_every_model_cites_a_symbol():
    """The citation check is not passing on models that cite nothing."""
    uncited = [
        path.name
        for path in SPEC_FILES
        if path.suffix != ".matrix" and not SYMBOL_CITATION.search(path.read_text())
    ]

    assert uncited == []


class Raises:
    """The model's ``none``: the Python raises this exception."""

    def __init__(self, exception):
        self.exception = exception


def _halved(rows):
    """Undo HeatmapCells' doubling of every edge."""
    return [tuple(value / 2 for value in row) for row in rows]


def _gene_rows(genes, start, end):
    frame = pd.DataFrame(genes, columns=["start", "end"])
    return assign_gene_positions(frame, start, end)


def _genome(chroms):
    """One frame holding the positions of chromosomes "1", "2", ... in order."""
    return pd.DataFrame(
        {
            "chr": [
                str(i + 1) for i, positions in enumerate(chroms) for _ in positions
            ],
            "pos": [pos for positions in chroms for pos in positions],
            "p_value": 0.01,
        }
    )


def _points_are_laid_out(gap, chroms, i, p, j, q):
    """GenomeLayout's ``SeparationOK`` and ``TotalOK`` for points (i, p), (j, q)."""
    order = [str(n + 1) for n in range(len(chroms))]
    layout = GenomeLayout.from_frames(
        [_genome(chroms)], chrom_col="chr", pos_col="pos", order=order, gap=gap
    )
    x_p = layout.offsets[order[i]] + p
    x_q = layout.offsets[order[j]] + q
    last = order[-1]
    return (
        x_p + gap < x_q
        and layout.offsets[order[0]] == 0
        and layout.total_length == layout.offsets[last] + max(chroms[-1]) + gap
    )


def _miami_highlight(chroms, gap, chrom, start, end):
    """The x-span ``plot_miami`` shades for one region, or None when it skips it."""
    frame = _genome(chroms)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)  # a skipped region warns
        figure = MiamiPlotter(species="canine").plot_miami(
            frame,
            frame,
            highlight_regions=[(str(chrom + 1), start, end)],
            style=GenomeWideStyle(chrom_gap=gap),
        )
    bands = PROBES["matplotlib"].region_highlights(figure, 0)
    return (bands[0].x0, bands[0].x1) if bands else None


def _retry_loop(max_retries, retry_delay, errors):
    """Run ``_with_retries`` over failing attempts: how it ended, calls, sleeps."""
    script = iter(errors)
    calls = []

    def attempt():
        calls.append(1)
        raise next(script)

    with patch.object(_http.time, "sleep") as sleep:
        try:
            _http._with_retries(
                attempt, what="GET", max_retries=max_retries, retry_delay=retry_delay
            )
        except ValueError:
            ended = "rejected"
        except requests.RequestException as error:
            ended = ("raised", error.attempts)
    return ended, len(calls), [call.args[0] for call in sleep.call_args_list]


def _request_json(max_retries, outcomes):
    """Run ``request_json`` over failing outcomes: the branch taken and the calls.

    An outcome is an exception to raise or an HTTP status to answer with.
    """
    responses = [
        outcome
        if isinstance(outcome, Exception)
        else MagicMock(ok=False, status_code=outcome, text="body")
        for outcome in outcomes
    ]
    with (
        patch.object(_http.time, "sleep"),
        patch.object(_http.requests, "get", side_effect=responses) as get,
    ):
        try:
            _http.request_json(
                "https://example.invalid",
                {},
                error_cls=EnsemblAPIError,
                service="X",
                max_retries=max_retries,
                retry_delay=1,
            )
        except ValueError:
            return "rejected"
        except EnsemblAPIError as error:
            counted = re.search(r"failed after (\d+) attempts", str(error))
            branch = ("count", int(counted[1])) if counted else "fatal"
            return branch, get.call_count


def _panel(row, col, n_cols):
    return _Panel(None, row, col, n_cols)


def _cross_hits(n_rows, n_cols):
    """(p, q): panel p's secondary axis has the name of panel q's primary axis."""
    cells = [(r, c) for r in range(1, n_rows + 1) for c in range(1, n_cols + 1)]
    return [
        (p, q)
        for p in cells
        for q in cells
        if secondary_axis_key(_panel(*p, n_cols).secondary_ref())
        == _panel(*q, n_cols).axis("yaxis")
    ]


CONNECTION = requests.ConnectionError("reset")
HG38 = [[10583, 248946422], [10019, 242183529], [9, 198235559]]
TWO_CHROMS = [[1000000, 100000000], [5000000, 60000000]]

# (model, the model's #guard, the Python on the same input, the pinned value)
VECTORS = [
    (
        "GeneRows",
        "assign_gene_positions [⟨1, 500⟩, ⟨600, 100⟩, ⟨200, 300⟩] 1 1001 = [0, 0, 0]",
        lambda: _gene_rows([(1, 500), (600, 100), (200, 300)], 1, 1001),
        [0, 0, 0],
    ),
    (
        "GeneRows",
        "assign_gene_positions [⟨1, 500⟩, ⟨600, 100⟩] 1 1001 = [0, 0]",
        lambda: _gene_rows([(1, 500), (600, 100)], 1, 1001),
        [0, 0],
    ),
    (
        "GenomeLayout",
        "allOK ⟨1000000, [[10583, 248946422], [10019, 242183529], [9, 198235559]],"
        " 0, 248946422, 1, 10019⟩",
        lambda: _points_are_laid_out(1000000, HG38, 0, 248946422, 1, 10019),
        True,
    ),
    (
        "GenomeLayout",
        "allOK ⟨1000000, [[10583, 248946422], [10019, 242183529], [9, 198235559]],"
        " 1, 242183529, 2, 9⟩",
        lambda: _points_are_laid_out(1000000, HG38, 1, 242183529, 2, 9),
        True,
    ),
    (
        "GenomeLayout",
        "(⟨1, [[1], [1]], 0, 1, 3⟩ : Region).span = (1, 1)",
        lambda: _miami_highlight([[1], [1]], 1, 0, 1, 3),
        (1, 1),
    ),
    (
        "GenomeLayout",
        "highlight (from_frames [[1], [1]] 1) 0 2 3 = none",
        lambda: _miami_highlight([[1], [1]], 1, 0, 2, 3),
        None,
    ),
    (
        "GenomeLayout",
        "highlight (from_frames [[1000000, 100000000], [5000000, 60000000]] 1000000)"
        " 0 150000000 160000000 = none",
        lambda: _miami_highlight(TWO_CHROMS, 1000000, 0, 150000000, 160000000),
        None,
    ),
    (
        "GenomeLayout",
        "highlight (from_frames [[1000000, 100000000], [5000000, 60000000]] 1000000)"
        " 0 90000000 160000000 = some (90000000, 100000000)",
        lambda: _miami_highlight(TWO_CHROMS, 1000000, 0, 90000000, 160000000),
        (90000000, 100000000),
    ),
    (
        "HeatmapCells",
        "cell_edges [] == none",
        lambda: cell_edges([]),
        Raises(IndexError),
    ),
    (
        "HeatmapCells",
        "cell_edges [5, 5, 7] == some [(10, 10), (10, 12), (12, 16)]",
        lambda: cell_edges([5, 5, 7]),
        _halved([(10, 10), (10, 12), (12, 16)]),
    ),
    (
        "HeatmapCells",
        "cell_edges [1, 5, 5, 9] == some [(-2, 6), (6, 10), (10, 14), (14, 22)]",
        lambda: cell_edges([1, 5, 5, 9]),
        _halved([(-2, 6), (6, 10), (10, 14), (14, 22)]),
    ),
    (
        "HeatmapCells",
        "cell_edges [3, 2, 1] == some [(7, 5), (5, 3), (3, 1)]",
        lambda: cell_edges([3, 2, 1]),
        _halved([(7, 5), (5, 3), (3, 1)]),
    ),
    (
        "HeatmapCells",
        "cell_edges [0, 10, 5] == some [(-10, 10), (10, 15), (15, 5)]",
        lambda: cell_edges([0, 10, 5]),
        _halved([(-10, 10), (10, 15), (15, 5)]),
    ),
    (
        "HeatmapCells",
        "heatmap_highlight_rects 0 [5, 5, 7] [0, 1, 2]"
        " == some [(10, -1, 0, 2), (10, 1, 0, 2), (10, 3, 0, 2)]",
        lambda: heatmap_highlight_rects(0, [5, 5, 7], [0, 1, 2]),
        _halved([(10, -1, 0, 2), (10, 1, 0, 2), (10, 3, 0, 2)]),
    ),
    (
        "HeatmapCells",
        "heatmap_highlight_rects 0 [1, 2, 3] [0, 1] == none",
        lambda: heatmap_highlight_rects(0, [1, 2, 3], [0, 1]),
        Raises(IndexError),
    ),
    (
        "PlotlyAxes",
        "crossHits 100 1 = [((1, 1), (100, 1))]",
        lambda: _cross_hits(100, 1),
        [((1, 1), (100, 1))],
    ),
    (
        "PlotlyAxes",
        "crossHits 50 2 = [((1, 1), (50, 2))]",
        lambda: _cross_hits(50, 2),
        [((1, 1), (50, 2))],
    ),
    ("PlotlyAxes", "(crossHits 99 1).isEmpty", lambda: _cross_hits(99, 1), []),
    ("PlotlyAxes", "(crossHits 33 3).isEmpty", lambda: _cross_hits(33, 3), []),
    (
        "PlotlyAxes",
        "secondary_ref 1 = some 100",
        lambda: _panel(1, 1, 1).secondary_ref(),
        "y100",
    ),
    (
        "PlotlyAxes",
        "subplot_idx 1 3 2 = subplot_idx 2 1 2",
        lambda: _panel(1, 3, 2).subplot_idx,
        _panel(2, 1, 2).subplot_idx,
    ),
    (
        "PlotlyAxes",
        "subplot_idx 2 0 2 = subplot_idx 1 2 2",
        lambda: _panel(2, 0, 2).subplot_idx,
        _panel(1, 2, 2).subplot_idx,
    ),
    (
        "PlotlyAxes",
        "axis (subplot_idx 1 0 2) = axis (subplot_idx 1 1 2)",
        lambda: _panel(1, 0, 2).axis("yaxis"),
        _panel(1, 1, 2).axis("yaxis"),
    ),
    (
        "PlotlyAxes",
        "split_pixels_even 800 0 = none",
        lambda: split_pixels(800, None, 0),
        Raises(ZeroDivisionError),
    ),
    (
        "PlotlyAxes",
        "split_pixels_even 800 3 = some [266, 266, 266]",
        lambda: split_pixels(800, None, 3),
        [266, 266, 266],
    ),
    (
        "PlotlyAxes",
        "split_pixels_ratio 800 [3, 1] = some [600, 200]",
        lambda: split_pixels(800, [3, 1], 2),
        [600, 200],
    ),
    (
        "PlotlyAxes",
        "split_pixels_ratio 800 [0, 0] = none",
        lambda: split_pixels(800, [0, 0], 2),
        Raises(ZeroDivisionError),
    ),
    (
        "PlotlyAxes",
        "split_pixels_ratio 800 [] = some []",
        lambda: split_pixels(800, [], 0),
        [],
    ),
    (
        "RetryLoop",
        "_with_retries 0 1 [.error .connection] = (.rejected, init)",
        lambda: _retry_loop(0, 1, [CONNECTION]),
        ("rejected", 0, []),
    ),
    (
        "RetryLoop",
        "_with_retries (-1) 1 [.error .connection] = (.rejected, init)",
        lambda: _retry_loop(-1, 1, [CONNECTION]),
        ("rejected", 0, []),
    ),
    (
        "RetryLoop",
        "_with_retries 3 (-1) [.error .connection] = (.rejected, init)",
        lambda: _retry_loop(3, -1, [CONNECTION]),
        ("rejected", 0, []),
    ),
    (
        "RetryLoop",
        "(request_json 0 1 [.error .connection]).1 = .rejected",
        lambda: _request_json(0, [CONNECTION]),
        "rejected",
    ),
    (
        "RetryLoop",
        "(request_json 3 1 [.error .connection, .error .connection,"
        " .error (.http (some 429))]).1 = .countBranch 3 3",
        lambda: _request_json(3, [CONNECTION, CONNECTION, 429]),
        (("count", 3), 3),
    ),
    (
        "RetryLoop",
        "(request_json 3 1 [.error (.http (some 429)), .error (.http (some 429)),"
        " .error .connection]).1 = .countBranch 3 3",
        lambda: _request_json(3, [429, 429, CONNECTION]),
        (("count", 3), 3),
    ),
    (
        "RetryLoop",
        "(request_json 3 1 [.error (.http (some 404))]).1"
        " = .fatalBranch 1 (.http (some 404))",
        lambda: _request_json(3, [404]),
        ("fatal", 1),
    ),
    (
        "RetryLoop",
        "(request_json 3 1 [.error .connection, .error (.http (some 404))]).1"
        " = .fatalBranch 2 (.http (some 404))",
        lambda: _request_json(3, [CONNECTION, 404]),
        ("fatal", 2),
    ),
    (
        "RetryLoop",
        "(runLoop 3 [.error .connection, .error .connection, .error .connection])"
        ".2.sleeps = [1, 2]",
        lambda: _retry_loop(3, 1, [CONNECTION] * 3),
        (("raised", 3), 3, [1, 2]),
    ),
]


@pytest.mark.parametrize(
    ("run", "pinned"),
    [
        pytest.param(run, pinned, id=f"{model}: {guard}")
        for model, guard, run, pinned in VECTORS
    ],
)
def test_python_returns_what_the_model_pins(run, pinned):
    """A Python change that leaves a modelled value fails on the model's vector."""
    if isinstance(pinned, Raises):
        with pytest.raises(pinned.exception):
            run()
    else:
        assert run() == pinned


def _guards(model):
    """Every ``#guard`` of one Lean model, on one line with single spaces."""
    text = (ROOT / "specs" / "lean" / f"{model}.lean").read_text()
    statements = re.findall(r"^#guard (.*(?:\n[ \t]+\S.*)*)", text, flags=re.MULTILINE)
    return [" ".join(statement.split()) for statement in statements]


def _head(guard):
    """The first name in a guard, such as ``cell_edges``."""
    return re.search(r"[A-Za-z_][\w.]*", guard)[0]


def test_vector_table_matches_the_models():
    """A row quotes a guard the model still has, and no vector guard lacks a row."""
    quoted = {(model, " ".join(guard.split())) for model, guard, _, _ in VECTORS}
    heads = {(model, _head(guard)) for model, guard in quoted}
    in_models = {
        (model, guard)
        for model in {model for model, _ in quoted}
        for guard in _guards(model)
        if (model, _head(guard)) in heads
    }

    assert quoted == in_models
