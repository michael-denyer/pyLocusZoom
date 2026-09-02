"""Guard the one-name-one-schema rule for pytest fixtures.

A fixture name that means two different DataFrame shapes forces every reader
to walk up to the nearest enclosing class before they can tell what a test
runs on. This test fails when a name regains a second schema.
"""

import ast
import collections
from pathlib import Path

TESTS_DIR = Path(__file__).parent


def _return_schema(fn: ast.FunctionDef) -> tuple[str, ...]:
    """Return the column names of every DataFrame a fixture builds.

    Reading the return expression alone misses a fixture that assembles frames
    into a list and returns the name, so every such fixture reported the empty
    schema and a second shape under the same name went unnoticed. Collecting
    the keys of each ``pd.DataFrame({...})`` literal in the body catches those.
    """
    columns: set[str] = set()
    for node in ast.walk(fn):
        if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)):
            continue
        if node.func.attr != "DataFrame" or not node.args:
            continue
        frame = node.args[0]
        if not isinstance(frame, ast.Dict):
            continue
        columns.update(key.value for key in frame.keys if isinstance(key, ast.Constant))
    return tuple(sorted(columns))


def _fixture_schemas() -> dict[str, set[tuple[str, ...]]]:
    schemas: dict[str, set[tuple[str, ...]]] = collections.defaultdict(set)
    for path in sorted(TESTS_DIR.glob("*.py")):
        tree = ast.parse(path.read_text())
        for node in ast.walk(tree):
            if not isinstance(node, ast.FunctionDef):
                continue
            if any("fixture" in ast.unparse(d) for d in node.decorator_list):
                schemas[node.name].add(_return_schema(node))
    return schemas


def test_every_fixture_name_has_one_schema():
    """No fixture name may describe two different data shapes."""
    offenders = {
        name: shapes for name, shapes in _fixture_schemas().items() if len(shapes) > 1
    }

    assert not offenders, "Fixture names bound to more than one schema: " + "; ".join(
        f"{name} -> {sorted(shapes)}" for name, shapes in offenders.items()
    )
