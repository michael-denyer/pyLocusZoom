"""Tests for package exports."""


def test_ensure_recomb_maps_exported():
    """Test that ensure_recomb_maps is exported from package."""
    from pylocuszoom import ensure_recomb_maps

    assert callable(ensure_recomb_maps)


def test_prepare_finemapping_for_plotting_exported():
    """Test that prepare_finemapping_for_plotting is exported from package."""
    from pylocuszoom import prepare_finemapping_for_plotting

    assert callable(prepare_finemapping_for_plotting)


def test_version_matches_pyproject():
    """Package __version__ should match the version in pyproject.toml."""
    import re
    from pathlib import Path

    import pylocuszoom

    pyproject = Path(__file__).parent.parent / "pyproject.toml"
    text = pyproject.read_text()
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
    assert match, "Could not find version in pyproject.toml"

    expected = match.group(1)
    assert pylocuszoom.__version__ == expected, (
        f"__version__={pylocuszoom.__version__!r} != pyproject.toml version={expected!r}"
    )
