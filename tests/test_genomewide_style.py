"""GenomeWideStyle: palette, points, fonts and chromosome axis of genome-wide plots."""

import colorcet as cc
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import pytest
from bokeh.models import Scatter
from matplotlib.collections import PathCollection
from pydantic import ValidationError as PydanticValidationError

from pylocuszoom import GenomeWideStyle, ManhattanPlotter, MiamiPlotter
from pylocuszoom.manhattan import CHROMOSOME_GAP

DEFAULT_PALETTE = list(cc.b_glasbey_bw_minc_20_maxl_70)


@pytest.fixture
def four_chrom_df():
    """Four chromosomes, ten SNPs each, with a strong hit on chromosome 2."""
    rng = np.random.default_rng(7)
    return pd.DataFrame(
        {
            "chr": np.repeat([1, 2, 3, 4], 10),
            "pos": np.tile(np.arange(1, 11) * 1_000_000, 4),
            "p_value": np.concatenate(
                [rng.uniform(1e-3, 1, 10), [1e-9], rng.uniform(1e-3, 1, 29)]
            ),
        }
    )


@pytest.fixture
def plotter():
    return ManhattanPlotter(species="human")


def _manhattan_collections(ax):
    """Scatter collections of a Manhattan panel, one per chromosome, in order."""
    return [c for c in ax.collections if isinstance(c, PathCollection)]


def _hex(rgba):
    return mcolors.to_hex(rgba[:3])


def _xtick_texts(ax):
    return [tick.get_text() for tick in ax.get_xticklabels()]


class TestModel:
    def test_defaults(self):
        style = GenomeWideStyle()

        assert style.palette is None
        assert style.point_size is None
        assert style.point_alpha is None
        assert style.title_fontsize == 14
        assert style.panel_title_fontsize is None
        assert style.axis_label_fontsize is None
        assert style.tick_label_fontsize is None
        assert style.tick_step == 1
        assert style.tick_rotation is None
        assert style.chrom_gap == CHROMOSOME_GAP == 1_000_000
        assert style.line_style == "--"
        assert style.line_width == 1.0
        assert style.title_fontweight == "bold"
        assert style.point_edge_width is None
        assert style.y_headroom == 0.1
        assert style.manhattan_qq_width_ratio == 2.5

    def test_is_frozen(self):
        with pytest.raises(PydanticValidationError):
            GenomeWideStyle().tick_step = 2

    @pytest.mark.parametrize(
        "field, value",
        [
            ("tick_step", 0),
            ("point_size", 0),
            ("point_alpha", 0),
            ("point_alpha", 1.5),
            ("title_fontsize", 0),
            ("panel_title_fontsize", -1),
            ("axis_label_fontsize", 0),
            ("tick_label_fontsize", 0),
            ("chrom_gap", -1),
            ("palette", ()),
            ("palette", ("#d60000", "not-a-colour")),
            ("line_style", "dashed"),
            ("line_width", 0),
            ("title_fontweight", "heavy"),
            ("point_edge_width", -0.1),
            ("y_headroom", -0.1),
            ("manhattan_qq_width_ratio", 0),
        ],
    )
    def test_rejects_invalid_values(self, field, value):
        with pytest.raises(PydanticValidationError):
            GenomeWideStyle(**{field: value})

    def test_palette_is_normalised_to_hex(self):
        assert GenomeWideStyle(palette=["red", (0, 0, 1)]).palette == (
            "#ff0000",
            "#0000ff",
        )

    def test_is_exported(self):
        import pylocuszoom

        assert "GenomeWideStyle" in pylocuszoom.__all__


class TestDefaultsAreUnchanged:
    """With no style, and with GenomeWideStyle(), figures render as before."""

    def test_manhattan_matplotlib_defaults(self, plotter, four_chrom_df):
        fig = plotter.plot_manhattan(four_chrom_df)

        (ax,) = fig.get_axes()
        collections = _manhattan_collections(ax)
        assert [_hex(c.get_facecolors()[0]) for c in collections] == (
            DEFAULT_PALETTE[:4]
        )
        assert {float(s) for c in collections for s in c.get_sizes()} == {10.0}
        assert {c.get_alpha() for c in collections} == {None}
        assert _xtick_texts(ax) == ["1", "2", "3", "4"]
        assert {t.get_fontsize() for t in ax.get_xticklabels()} == {8}
        assert {t.get_rotation() for t in ax.get_xticklabels()} == {0}
        assert ax.title.get_fontsize() == 14
        assert ax.xaxis.label.get_fontsize() == 12
        assert ax.yaxis.label.get_fontsize() == 12
        # The second chromosome starts one default gap past the first's end.
        second = collections[1].get_offsets()[:, 0].min()
        assert second == 10_000_000 + CHROMOSOME_GAP + 1_000_000

    def test_manhattan_qq_matplotlib_defaults(self, plotter, four_chrom_df):
        fig = plotter.plot_manhattan_qq(four_chrom_df, title="Study")

        manhattan_ax, qq_ax = fig.get_axes()
        assert manhattan_ax.title.get_fontsize() == 12
        assert qq_ax.title.get_fontsize() == 12
        assert qq_ax.xaxis.label.get_fontsize() == 12
        assert fig._suptitle.get_fontsize() == 14

    @pytest.mark.parametrize(
        "method",
        [
            "plot_manhattan",
            "plot_qq",
            "plot_manhattan_qq",
            "plot_manhattan_stacked",
            "plot_manhattan_qq_stacked",
        ],
    )
    @pytest.mark.parametrize("backend", ["plotly", "bokeh"])
    def test_explicit_default_style_matches_no_style(
        self, method, backend, four_chrom_df
    ):
        plotter = ManhattanPlotter(species="human", backend=backend)
        frames = (
            [four_chrom_df, four_chrom_df]
            if method.endswith("stacked")
            else four_chrom_df
        )

        plain = getattr(plotter, method)(frames, title="T")
        styled = getattr(plotter, method)(frames, title="T", style=GenomeWideStyle())

        if backend == "plotly":
            assert styled.to_json() == plain.to_json()
        else:
            assert _bokeh_signature(styled) == _bokeh_signature(plain)


def _bokeh_plots(layout):
    """Every bokeh Plot in a layout, depth first in layout order."""
    from bokeh.models import Plot

    if isinstance(layout, Plot):
        return [layout]
    return [plot for child in layout.children for plot in _bokeh_plots(child)]


def _bokeh_scatter_glyphs(plot):
    """The drawn Scatter glyph of each renderer, without selection variants."""
    from bokeh.models import GlyphRenderer

    return [
        renderer.glyph
        for renderer in plot.renderers
        if isinstance(renderer, GlyphRenderer) and isinstance(renderer.glyph, Scatter)
    ]


def _bokeh_signature(layout):
    """Every styling property this feature can touch, for a whole bokeh layout."""
    return [
        (
            p.title.text,
            p.title.text_font_size,
            [a.major_label_text_font_size for a in p.xaxis + p.yaxis],
            [a.axis_label_text_font_size for a in p.xaxis + p.yaxis],
            [a.major_label_orientation for a in p.xaxis],
            [(g.fill_alpha, g.line_alpha) for g in _bokeh_scatter_glyphs(p)],
        )
        for p in _bokeh_plots(layout)
    ]


class TestPalette:
    def test_colours_cycle_over_chromosomes(self, plotter, four_chrom_df):
        style = GenomeWideStyle(palette=["#111111", "#222222", "#333333"])

        fig = plotter.plot_manhattan(four_chrom_df, style=style)

        collections = _manhattan_collections(fig.get_axes()[0])
        assert [_hex(c.get_facecolors()[0]) for c in collections] == [
            "#111111",
            "#222222",
            "#333333",
            "#111111",
        ]

    def test_categorical_plot_takes_the_palette(self, plotter):
        df = pd.DataFrame({"cat": ["a", "b", "c"], "p_value": [0.1, 0.01, 0.001]})

        fig = plotter.plot_manhattan(
            df, category_col="cat", style=GenomeWideStyle(palette=["#abcdef"])
        )

        collections = _manhattan_collections(fig.get_axes()[0])
        assert {_hex(c.get_facecolors()[0]) for c in collections} == {"#abcdef"}

    def test_plotly_markers_take_the_palette(self, four_chrom_df):
        plotter = ManhattanPlotter(species="human", backend="plotly")

        fig = plotter.plot_manhattan(
            four_chrom_df, style=GenomeWideStyle(palette=["#111111", "#222222"])
        )

        assert [trace.marker.color for trace in fig.data] == [
            "#111111",
            "#222222",
            "#111111",
            "#222222",
        ]


class TestPoints:
    def test_point_size_sets_manhattan_and_qq_markers(self, plotter, four_chrom_df):
        fig = plotter.plot_manhattan_qq(
            four_chrom_df, style=GenomeWideStyle(point_size=75)
        )

        sizes = {
            float(s)
            for ax in fig.get_axes()
            for c in _manhattan_collections(ax)
            for s in c.get_sizes()
        }
        assert sizes == {75.0}

    def test_point_alpha_sets_matplotlib_alpha(self, plotter, four_chrom_df):
        fig = plotter.plot_manhattan_qq(
            four_chrom_df, style=GenomeWideStyle(point_alpha=0.4)
        )

        alphas = {
            c.get_alpha() for ax in fig.get_axes() for c in _manhattan_collections(ax)
        }
        assert alphas == {0.4}

    def test_point_alpha_sets_plotly_opacity(self, four_chrom_df):
        plotter = ManhattanPlotter(species="human", backend="plotly")

        fig = plotter.plot_manhattan(
            four_chrom_df, style=GenomeWideStyle(point_alpha=0.4)
        )

        assert {trace.marker.opacity for trace in fig.data} == {0.4}

    def test_point_alpha_sets_bokeh_glyph_alpha(self, four_chrom_df):
        plotter = ManhattanPlotter(species="human", backend="bokeh")

        fig = plotter.plot_manhattan(
            four_chrom_df, style=GenomeWideStyle(point_alpha=0.4)
        )

        glyphs = [g for p in _bokeh_plots(fig) for g in _bokeh_scatter_glyphs(p)]
        assert glyphs
        assert {(g.fill_alpha, g.line_alpha) for g in glyphs} == {(0.4, 0.4)}


class TestFonts:
    STYLE = GenomeWideStyle(
        title_fontsize=30,
        panel_title_fontsize=26,
        axis_label_fontsize=20,
        tick_label_fontsize=18,
    )

    def test_manhattan_qq_matplotlib_fonts(self, plotter, four_chrom_df):
        fig = plotter.plot_manhattan_qq(four_chrom_df, title="Study", style=self.STYLE)

        assert fig._suptitle.get_fontsize() == 30
        assert fig.get_axes()[1].title.get_fontsize() == 26
        for ax in fig.get_axes():
            assert ax.xaxis.label.get_fontsize() == 20
            assert ax.yaxis.label.get_fontsize() == 20
            assert {t.get_fontsize() for t in ax.get_xticklabels()} == {18}
            assert {t.get_fontsize() for t in ax.get_yticklabels()} == {18}

    def test_single_plot_title_is_a_panel_title(self, plotter, four_chrom_df):
        fig = plotter.plot_qq(
            four_chrom_df, style=GenomeWideStyle(panel_title_fontsize=9)
        )

        assert fig.get_axes()[0].title.get_fontsize() == 9

    def test_stacked_title_takes_title_fontsize(self, plotter, four_chrom_df):
        fig = plotter.plot_manhattan_stacked(
            [four_chrom_df, four_chrom_df],
            title="Stack",
            style=GenomeWideStyle(title_fontsize=22),
        )

        assert fig.get_axes()[0].title.get_fontsize() == 22

    def test_plotly_fonts(self, four_chrom_df):
        plotter = ManhattanPlotter(species="human", backend="plotly")

        fig = plotter.plot_manhattan_qq(four_chrom_df, title="Study", style=self.STYLE)

        layout = fig.layout
        assert layout.title.font.size == 30
        assert {a.font.size for a in layout.annotations} == {26}
        for axis in (layout.xaxis, layout.yaxis, layout.xaxis2, layout.yaxis2):
            assert axis.title.font.size == 20
            assert axis.tickfont.size == 18

    def test_bokeh_fonts(self, four_chrom_df):
        plotter = ManhattanPlotter(species="human", backend="bokeh")

        fig = plotter.plot_manhattan_qq(four_chrom_df, title="Study", style=self.STYLE)

        plots = _bokeh_plots(fig)
        assert plots[0].title.text == "Study"
        assert plots[0].title.text_font_size == "30pt"
        assert plots[1].title.text_font_size == "26pt"
        for plot in plots:
            for axis in plot.xaxis + plot.yaxis:
                assert axis.axis_label_text_font_size == "20pt"
                assert axis.major_label_text_font_size == "18pt"


class TestLines:
    STYLE = GenomeWideStyle(line_style="-", line_width=1.5)

    def test_threshold_lines_and_qq_diagonal_take_the_style(
        self, plotter, four_chrom_df
    ):
        fig = plotter.plot_manhattan_qq(
            four_chrom_df, suggestive_threshold=1e-5, style=self.STYLE
        )

        manhattan_ax, qq_ax = fig.get_axes()
        lines = manhattan_ax.get_lines() + qq_ax.get_lines()
        assert len(lines) == 3
        assert {(line.get_linestyle(), line.get_linewidth()) for line in lines} == {
            ("-", 1.5)
        }

    def test_miami_threshold_lines_take_the_style(self, four_chrom_df):
        fig = MiamiPlotter(species="human").plot_miami(
            four_chrom_df, four_chrom_df, style=self.STYLE
        )

        lines = [line for ax in fig.get_axes() for line in ax.get_lines()]
        assert lines
        assert {(line.get_linestyle(), line.get_linewidth()) for line in lines} == {
            ("-", 1.5)
        }

    def test_plotly_threshold_line_takes_the_style(self, four_chrom_df):
        plotter = ManhattanPlotter(species="human", backend="plotly")

        fig = plotter.plot_manhattan(four_chrom_df, style=self.STYLE)

        (shape,) = fig.layout.shapes
        assert (shape.line.dash, shape.line.width) == ("solid", 1.5)


class TestLayout:
    @pytest.mark.parametrize("headroom", [0.1, 0.3])
    def test_headroom_is_left_above_the_significance_line(
        self, plotter, four_chrom_df, headroom
    ):
        fig = plotter.plot_manhattan(
            four_chrom_df,
            significance_threshold=1e-12,
            style=GenomeWideStyle(y_headroom=headroom),
        )

        assert fig.axes[0].get_ylim()[1] == pytest.approx(12 * (1 + headroom))

    def test_width_ratio_sets_manhattan_to_qq_width(self, plotter, four_chrom_df):
        fig = plotter.plot_manhattan_qq(
            four_chrom_df, style=GenomeWideStyle(manhattan_qq_width_ratio=2)
        )

        manhattan, qq = fig.axes[:2]
        ratio = manhattan.get_position().width / qq.get_position().width
        assert ratio == pytest.approx(2)


class TestPointEdges:
    def test_zero_removes_manhattan_and_qq_outlines(self, plotter, four_chrom_df):
        fig = plotter.plot_manhattan_qq(
            four_chrom_df, style=GenomeWideStyle(point_edge_width=0)
        )

        widths = {
            float(w)
            for ax in fig.get_axes()
            for c in _manhattan_collections(ax)
            for w in c.get_linewidths()
        }
        assert widths == {0.0}

    def test_unset_keeps_each_method_width(self, plotter, four_chrom_df):
        fig = plotter.plot_manhattan_qq(four_chrom_df)

        manhattan_ax, qq_ax = fig.get_axes()
        assert {
            float(w)
            for c in _manhattan_collections(manhattan_ax)
            for w in c.get_linewidths()
        } == {0.1}
        assert {
            float(w) for c in _manhattan_collections(qq_ax) for w in c.get_linewidths()
        } == {0.02}


class TestTitleWeight:
    STYLE = GenomeWideStyle(title_fontweight="normal")

    def test_default_titles_are_bold(self, plotter, four_chrom_df):
        fig = plotter.plot_manhattan_qq(four_chrom_df, title="Study")

        assert fig._suptitle.get_fontweight() == "bold"
        assert fig.get_axes()[1].title.get_fontweight() == "bold"

    def test_normal_weight_on_suptitle_and_panel_titles(self, plotter, four_chrom_df):
        fig = plotter.plot_manhattan_qq(four_chrom_df, title="Study", style=self.STYLE)

        assert fig._suptitle.get_fontweight() == "normal"
        assert fig.get_axes()[1].title.get_fontweight() == "normal"

    def test_normal_weight_on_stacked_title(self, plotter, four_chrom_df):
        fig = plotter.plot_manhattan_stacked(
            [four_chrom_df, four_chrom_df], title="Stack", style=self.STYLE
        )

        assert fig.get_axes()[0].title.get_fontweight() == "normal"

    def test_plotly_panel_titles_are_not_bold(self, four_chrom_df):
        plotter = ManhattanPlotter(species="human", backend="plotly")

        fig = plotter.plot_manhattan_qq(four_chrom_df, title="Study", style=self.STYLE)

        texts = [a.text for a in fig.layout.annotations]
        assert texts
        assert not any("<b>" in text for text in texts)

    def test_bokeh_titles_are_not_bold(self, four_chrom_df):
        plotter = ManhattanPlotter(species="human", backend="bokeh")

        fig = plotter.plot_manhattan_qq(four_chrom_df, title="Study", style=self.STYLE)

        titled = [p.title for p in _bokeh_plots(fig) if p.title.text]
        assert titled
        assert {t.text_font_style for t in titled} == {"normal"}


class TestChromosomeAxis:
    def test_tick_step_labels_every_nth_chromosome(self, plotter, four_chrom_df):
        fig = plotter.plot_manhattan(four_chrom_df, style=GenomeWideStyle(tick_step=2))

        assert _xtick_texts(fig.get_axes()[0]) == ["1", "3"]

    def test_tick_step_on_plotly(self, four_chrom_df):
        plotter = ManhattanPlotter(species="human", backend="plotly")

        fig = plotter.plot_manhattan(four_chrom_df, style=GenomeWideStyle(tick_step=3))

        assert list(fig.layout.xaxis.ticktext) == ["1", "4"]

    def test_tick_rotation(self, plotter, four_chrom_df):
        fig = plotter.plot_manhattan(
            four_chrom_df, style=GenomeWideStyle(tick_rotation=90)
        )

        assert {t.get_rotation() for t in fig.get_axes()[0].get_xticklabels()} == {90}

    def test_tick_rotation_on_plotly(self, four_chrom_df):
        plotter = ManhattanPlotter(species="human", backend="plotly")

        fig = plotter.plot_manhattan(
            four_chrom_df, style=GenomeWideStyle(tick_rotation=90)
        )

        assert fig.layout.xaxis.tickangle == -90

    def test_chrom_gap_moves_the_next_chromosome(self, plotter, four_chrom_df):
        fig = plotter.plot_manhattan(four_chrom_df, style=GenomeWideStyle(chrom_gap=0))

        collections = _manhattan_collections(fig.get_axes()[0])
        assert collections[1].get_offsets()[:, 0].min() == 10_000_000 + 1_000_000

    def test_miami_takes_the_style(self, four_chrom_df):
        fig = MiamiPlotter(species="human").plot_miami(
            four_chrom_df,
            four_chrom_df,
            style=GenomeWideStyle(tick_step=2, palette=["#123456"], point_size=5),
        )

        top, bottom = fig.get_axes()
        assert _xtick_texts(bottom) == ["1", "3"]
        colours = {
            _hex(c.get_facecolors()[0])
            for ax in (top, bottom)
            for c in _manhattan_collections(ax)
        }
        assert colours == {"#123456"}
        sizes = {
            float(s)
            for ax in (top, bottom)
            for c in _manhattan_collections(ax)
            for s in c.get_sizes()
        }
        assert sizes == {5.0}


def test_reproduces_lab_disco_manhattan_style(four_chrom_df):
    """Recreate lab_disco's ManhattanStyle defaults with GenomeWideStyle.

    lab_disco hands qqman ``cmap=cc.m_glasbey_dark, cmap_var=42``; qqman keeps
    every ``int(256 / 42)``-th colour of the map, which is ``glasbey_dark[::6]``.
    """
    lab_style = GenomeWideStyle(
        palette=cc.glasbey_dark[::6],
        point_size=75,
        title_fontsize=30,
        panel_title_fontsize=26,
        axis_label_fontsize=20,
        tick_label_fontsize=18,
        tick_step=2,
        tick_rotation=90,
    )

    fig = ManhattanPlotter(species="human").plot_manhattan_qq(
        four_chrom_df,
        title="Lab study",
        suggestive_threshold=1e-5,
        style=lab_style,
    )

    manhattan_ax, qq_ax = fig.get_axes()
    collections = _manhattan_collections(manhattan_ax)
    assert [_hex(c.get_facecolors()[0]) for c in collections] == [
        "#d60000",
        cc.glasbey_dark[6],
        cc.glasbey_dark[12],
        cc.glasbey_dark[18],
    ]
    assert _xtick_texts(manhattan_ax) == ["1", "3"]
    assert {t.get_rotation() for t in manhattan_ax.get_xticklabels()} == {90}
    assert {t.get_fontsize() for t in manhattan_ax.get_xticklabels()} == {18}
    assert {t.get_fontsize() for t in qq_ax.get_yticklabels()} == {18}
    assert manhattan_ax.xaxis.label.get_fontsize() == 20
    assert qq_ax.title.get_fontsize() == 26
    assert fig._suptitle.get_fontsize() == 30
    assert {float(s) for c in collections for s in c.get_sizes()} == {75.0}
