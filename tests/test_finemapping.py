"""Tests for fine-mapping/SuSiE data handling."""

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from pylocuszoom import DisplayConfig, PanelInputs
from pylocuszoom.backends.hover import HoverConfig
from pylocuszoom.backends.matplotlib_backend import MatplotlibBackend
from pylocuszoom.config import RegionConfig
from pylocuszoom.finemapping import (
    FinemappingValidationError,
    filter_by_credible_set,
    filter_finemapping_by_region,
    get_credible_sets,
    get_top_pip_variants,
    prepare_finemapping_for_plotting,
    validate_finemapping_df,
)
from pylocuszoom.panels.finemapping import FinemappingPanel
from pylocuszoom.plotter import LocusZoomPlotter

DRAW_REGION = RegionConfig(chrom=1, start=1, end=1_000_000)


@pytest.fixture
def finemapping_df():
    """Create sample fine-mapping DataFrame."""
    return pd.DataFrame(
        {
            "pos": [1000, 2000, 3000, 4000, 5000],
            "pip": [0.95, 0.02, 0.8, 0.1, 0.01],
            "cs": [1, 0, 2, 2, 0],
            "rs": ["rs1", "rs2", "rs3", "rs4", "rs5"],
        }
    )


class TestValidateFinemappingDf:
    """Tests for validate_finemapping_df function."""

    def test_valid_df(self, finemapping_df):
        """Should not raise for valid DataFrame."""
        validate_finemapping_df(finemapping_df)

    def test_missing_pos_col(self, finemapping_df):
        """Should raise for missing position column."""
        df = finemapping_df.drop(columns=["pos"])
        with pytest.raises(FinemappingValidationError, match="Missing columns"):
            validate_finemapping_df(df)

    def test_missing_pip_col(self, finemapping_df):
        """Should raise for missing PIP column."""
        df = finemapping_df.drop(columns=["pip"])
        with pytest.raises(FinemappingValidationError, match="Missing columns"):
            validate_finemapping_df(df)

    def test_invalid_pip_values(self, finemapping_df):
        """Should raise for PIP values outside [0, 1]."""
        df = finemapping_df.copy()
        df.loc[0, "pip"] = 1.5
        with pytest.raises(FinemappingValidationError, match="values > 1"):
            validate_finemapping_df(df)

    def test_custom_column_names(self):
        """Should accept custom column names."""
        df = pd.DataFrame({"position": [1000, 2000], "probability": [0.5, 0.3]})
        validate_finemapping_df(df, pos_col="position", pip_col="probability")


class TestFilterFinemappingByRegion:
    """Tests for filter_finemapping_by_region function."""

    def test_filter_by_position(self, finemapping_df):
        """Should filter to region bounds."""
        result = filter_finemapping_by_region(
            finemapping_df, chrom=1, start=1500, end=4500
        )
        assert len(result) == 3
        assert set(result["pos"]) == {2000, 3000, 4000}

    def test_filter_with_chrom(self):
        """Should filter by chromosome when column exists."""
        df = pd.DataFrame(
            {
                "chr": ["1", "1", "2"],
                "pos": [1000, 2000, 1500],
                "pip": [0.5, 0.3, 0.8],
            }
        )
        result = filter_finemapping_by_region(df, chrom=1, start=0, end=3000)
        assert len(result) == 2


class TestGetCredibleSets:
    """Tests for get_credible_sets function."""

    def test_returns_unique_cs(self, finemapping_df):
        """Should return sorted unique credible set IDs."""
        result = get_credible_sets(finemapping_df)
        assert result == [1, 2]

    def test_excludes_zero(self):
        """Should exclude cs=0 (not in credible set)."""
        df = pd.DataFrame({"pos": [1, 2, 3], "pip": [0.5, 0.3, 0.2], "cs": [0, 0, 0]})
        result = get_credible_sets(df)
        assert result == []

    def test_no_cs_column(self):
        """Should return empty list if no cs column."""
        df = pd.DataFrame({"pos": [1000, 2000], "pip": [0.5, 0.3]})
        result = get_credible_sets(df)
        assert result == []


class TestFilterByCredibleSet:
    """Tests for filter_by_credible_set function."""

    def test_filter_to_cs(self, finemapping_df):
        """Should filter to specific credible set."""
        result = filter_by_credible_set(finemapping_df, cs_id=2)
        assert len(result) == 2
        assert set(result["pos"]) == {3000, 4000}

    def test_missing_cs_col(self):
        """Should raise if cs column missing."""
        df = pd.DataFrame({"pos": [1000], "pip": [0.5]})
        with pytest.raises(FinemappingValidationError):
            filter_by_credible_set(df, cs_id=1)


class TestGetTopPipVariants:
    """Tests for get_top_pip_variants function."""

    def test_returns_top_n(self, finemapping_df):
        """Should return top N by PIP."""
        result = get_top_pip_variants(finemapping_df, n=2)
        assert len(result) == 2
        assert list(result["pip"]) == [0.95, 0.8]

    def test_respects_threshold(self, finemapping_df):
        """Should filter by PIP threshold."""
        result = get_top_pip_variants(finemapping_df, n=10, pip_threshold=0.5)
        assert len(result) == 2


class TestPrepareFinemappingForPlotting:
    """Tests for prepare_finemapping_for_plotting function."""

    def test_sorts_by_position(self, finemapping_df):
        """Should sort by position."""
        # Shuffle first
        df = finemapping_df.sample(frac=1, random_state=42)
        result = prepare_finemapping_for_plotting(df)
        assert list(result["pos"]) == sorted(finemapping_df["pos"])

    def test_filters_by_region(self, finemapping_df):
        """Should filter by region when specified."""
        result = prepare_finemapping_for_plotting(
            finemapping_df, chrom=1, start=1500, end=3500
        )
        assert len(result) == 2


class TestDrawFinemapping:
    """Tests for the fine-mapping panel's drawing.

    Assertions query the rendered matplotlib axes directly per
    CLAUDE.md's observable-outputs rule.
    """

    @pytest.fixture
    def rendering_axes(self):
        """Provide a real MatplotlibBackend and matplotlib Axes.

        Yields (backend, ax) and closes the figure after the test.
        """
        backend = MatplotlibBackend()
        fig, ax = plt.subplots()
        yield backend, ax

    def test_pip_line_carries_the_input_values(self, rendering_axes):
        """PIP values are rendered as a single line on the axes."""
        backend, ax = rendering_axes
        df = pd.DataFrame({"pos": [1000, 2000, 3000], "pip": [0.1, 0.5, 0.2]})

        FinemappingPanel.from_frame(df, DRAW_REGION, None).draw(backend, ax)

        lines = ax.get_lines()
        assert len(lines) >= 1, "expected a PIP line on the axes"
        y = list(lines[0].get_ydata())
        assert y == [0.1, 0.5, 0.2]

    def test_each_credible_set_gets_its_own_collection(self, rendering_axes):
        """Each credible set contributes its own scatter collection."""
        backend, ax = rendering_axes
        df = pd.DataFrame(
            {
                "pos": [1000, 2000, 3000, 4000],
                "pip": [0.1, 0.5, 0.2, 0.05],
                "cs": [1, 1, 2, 0],
            }
        )

        panel = FinemappingPanel.from_frame(df, DRAW_REGION, "cs")
        assert panel.credible_sets == [1, 2], "cs=0 is not a credible set"
        panel.draw(backend, ax)

        assert len(ax.collections) >= 2, (
            f"expected >=2 scatter collections for CS 1 and CS 2, "
            f"got {len(ax.collections)}"
        )

    def test_pip_line_renders_without_a_credible_set_column(self, rendering_axes):
        """PIP line renders even when no credible-set column is provided."""
        backend, ax = rendering_axes
        df = pd.DataFrame({"pos": [1000, 2000, 3000], "pip": [0.1, 0.5, 0.2]})

        FinemappingPanel.from_frame(df, DRAW_REGION, None).draw(backend, ax)

        assert len(ax.get_lines()) >= 1

    def test_only_points_above_the_pip_threshold_scatter(self, rendering_axes):
        """The panel scatters the variants that clear PIP_SCATTER_THRESHOLD."""
        backend, ax = rendering_axes
        df = pd.DataFrame({"pos": [1000, 2000, 3000], "pip": [0.005, 0.5, 0.002]})

        FinemappingPanel.from_frame(df, DRAW_REGION, None).draw(backend, ax)

        assert len(ax.get_lines()) >= 1
        assert len(ax.collections) == 1
        offsets = ax.collections[0].get_offsets()
        assert len(offsets) == 1
        assert offsets[0][0] == 2000 and offsets[0][1] == pytest.approx(0.5)

    def test_empty_frame_draws_nothing(self, rendering_axes):
        """An empty panel draws no line and no points."""
        backend, ax = rendering_axes
        panel = FinemappingPanel(
            data=pd.DataFrame({"pos": [], "pip": []}),
            height=1.5,
            cs_col=None,
            credible_sets=[],
            hover=HoverConfig(pos_col="pos", extra_cols={"pip": "PIP"}),
        )

        panel.draw(backend, ax)

        assert len(ax.get_lines()) == 0
        assert len(ax.collections) == 0


class TestPlotterDelegation:
    """Tests for plotter delegation to specialized classes."""

    def test_finemapping_panel_renders_the_supplied_pips(self):
        """plot_stacked() draws the fine-mapping frame onto its own PIP panel.

        Asserts on the plotted points rather than on a call to the renderer, so
        a change to the internal dispatch path cannot break this test without
        changing what the reader sees.
        """
        plotter = LocusZoomPlotter(species=None, backend="matplotlib", log_level=None)

        gwas_df = pd.DataFrame({"pos": [1000, 2000], "p_value": [0.01, 0.001]})
        fm_df = pd.DataFrame({"pos": [1000, 2000], "pip": [0.5, 0.3], "cs": [1, 1]})

        fig = plotter.plot_stacked(
            [gwas_df],
            chrom=1,
            start=1,
            end=3000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(finemapping_df=fm_df),
        )

        pip_axes = [ax for ax in fig.get_axes() if ax.get_ylabel() == "PIP"]
        assert len(pip_axes) == 1, "fine-mapping data should add one PIP panel"
        plotted = pip_axes[0].collections[0].get_offsets().tolist()
        assert plotted == [[1000.0, 0.5], [2000.0, 0.3]]


class TestFinemappingManyCredibleSets:
    """Test fine-mapping plot with many credible sets cycles colors correctly."""

    def test_finemapping_12_credible_sets_cycle_the_palette(
        self, regional_plotter, small_regional_gwas_df
    ):
        """Draw all 12 credible sets in the palette's 10 colours, cycled."""
        from pylocuszoom.colors import CREDIBLE_SET_COLORS

        n_variants = 60
        positions = list(range(1000000, 1000000 + n_variants * 10000, 10000))
        credible_sets = [((i // 5) % 12) + 1 for i in range(n_variants)]

        finemapping_df = pd.DataFrame(
            {
                "pos": positions,
                "pip": [0.8 if i % 5 == 0 else 0.1 for i in range(n_variants)],
                "cs": credible_sets,
            }
        )

        fig = regional_plotter.plot_stacked(
            [small_regional_gwas_df],
            chrom=1,
            start=900000,
            end=1700000,
            display=DisplayConfig(show_recombination=False),
            panels=PanelInputs(finemapping_df=finemapping_df, finemapping_cs_col="cs"),
        )

        pip_ax = fig.get_axes()[1]
        drawn_points = sum(len(c.get_offsets()) for c in pip_ax.collections)
        drawn_colours = {
            tuple(face) for c in pip_ax.collections for face in c.get_facecolor()
        }
        assert drawn_points == n_variants
        assert len(drawn_colours) == len(CREDIBLE_SET_COLORS)

    def test_credible_set_color_cycling(self):
        """Verify get_credible_set_color cycles correctly for cs > 10."""
        from pylocuszoom.colors import CREDIBLE_SET_COLORS, get_credible_set_color

        # Test that colors cycle after 10
        for cs_id in range(1, 25):
            color = get_credible_set_color(cs_id)
            expected_idx = (cs_id - 1) % len(CREDIBLE_SET_COLORS)
            expected_color = CREDIBLE_SET_COLORS[expected_idx]
            assert color == expected_color, (
                f"CS {cs_id}: got {color}, expected {expected_color}"
            )
