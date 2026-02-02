"""Tests for package exports."""


def test_ensure_recomb_maps_exported():
    """Test that ensure_recomb_maps is exported from package."""
    from pylocuszoom import ensure_recomb_maps

    assert callable(ensure_recomb_maps)


def test_plot_finemapping_exported():
    """Test that plot_finemapping is exported from package."""
    from pylocuszoom import plot_finemapping

    assert callable(plot_finemapping)
