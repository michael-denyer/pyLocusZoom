"""Tests for package exports."""


def test_ensure_recomb_maps_exported():
    """Test that ensure_recomb_maps is exported from package."""
    from pylocuszoom import ensure_recomb_maps

    assert callable(ensure_recomb_maps)


def test_plot_finemapping_exported():
    """Test that plot_finemapping is exported from package."""
    from pylocuszoom import plot_finemapping

    assert callable(plot_finemapping)


def test_version_matches_pyproject():
    """Package __version__ should match the version in pyproject.toml."""
    from pathlib import Path

    import tomllib

    import pylocuszoom

    pyproject = Path(__file__).parent.parent / "pyproject.toml"
    with open(pyproject, "rb") as f:
        meta = tomllib.load(f)

    expected = meta["project"]["version"]
    assert pylocuszoom.__version__ == expected, (
        f"__version__={pylocuszoom.__version__!r} != pyproject.toml version={expected!r}"
    )
