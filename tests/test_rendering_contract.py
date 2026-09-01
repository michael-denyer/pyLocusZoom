"""Contract tests for the semantic Manhattan/QQ renderer."""

from types import SimpleNamespace

import pandas as pd
import pytest

from pylocuszoom._miami_renderer import MiamiRenderer
from pylocuszoom._rendering import ManhattanQQRenderer
from pylocuszoom.backends import get_backend
from pylocuszoom.manhattan import prepare_manhattan_data
from pylocuszoom.qq import prepare_qq_data


class RecordingBackend:
    """Small primitive adapter used to test the renderer's semantic seam."""

    supports_hover = False

    def __init__(self):
        self.calls = []

    def _record(self, name, *args, **kwargs):
        self.calls.append((name, args, kwargs))
        return SimpleNamespace()

    def create_figure(self, **kwargs):
        self._record("create_figure", **kwargs)
        return SimpleNamespace(), [SimpleNamespace() for _ in range(kwargs["n_panels"])]

    def create_figure_grid(self, **kwargs):
        self._record("create_figure_grid", **kwargs)
        n_axes = kwargs["n_rows"] * kwargs["n_cols"]
        return SimpleNamespace(), [SimpleNamespace() for _ in range(n_axes)]

    def finalize_layout(self, *args, **kwargs):
        self._record("finalize_layout", *args, **kwargs)

    def scatter(self, *args, **kwargs):
        return self._record("scatter", *args, **kwargs)

    def fill_between(self, *args, **kwargs):
        return self._record("fill_between", *args, **kwargs)

    def line(self, *args, **kwargs):
        return self._record("line", *args, **kwargs)

    def axhline(self, *args, **kwargs):
        return self._record("axhline", *args, **kwargs)

    def set_xlim(self, *args, **kwargs):
        self._record("set_xlim", *args, **kwargs)

    def set_ylim(self, *args, **kwargs):
        self._record("set_ylim", *args, **kwargs)

    def set_xticks(self, *args, **kwargs):
        self._record("set_xticks", *args, **kwargs)

    def set_xlabel(self, *args, **kwargs):
        self._record("set_xlabel", *args, **kwargs)

    def set_ylabel(self, *args, **kwargs):
        self._record("set_ylabel", *args, **kwargs)

    def set_title(self, *args, **kwargs):
        self._record("set_title", *args, **kwargs)

    def set_suptitle(self, *args, **kwargs):
        self._record("set_suptitle", *args, **kwargs)

    def add_panel_label(self, *args, **kwargs):
        self._record("add_panel_label", *args, **kwargs)

    def set_yticks(self, *args, **kwargs):
        self._record("set_yticks", *args, **kwargs)

    def axvline(self, *args, **kwargs):
        return self._record("axvline", *args, **kwargs)

    def add_text(self, *args, **kwargs):
        return self._record("add_text", *args, **kwargs)

    def add_legend(self, *args, **kwargs):
        return self._record("add_legend", *args, **kwargs)


class FullCapabilityBackend(RecordingBackend):
    """RecordingBackend that also opts into every optional capability.

    RecordingBackend deliberately implements no capability protocol, so it
    doubles as the backend that declines them. Renderers whose whole figure is
    a heatmap or a forest plot need one that accepts.
    """

    def add_heatmap(self, *args, **kwargs):
        return self._record("add_heatmap", *args, **kwargs)

    def add_colorbar(self, *args, **kwargs):
        return self._record("add_colorbar", *args, **kwargs)

    def highlight_heatmap_snp(self, *args, **kwargs):
        self._record("highlight_heatmap_snp", *args, **kwargs)

    def hbar(self, *args, **kwargs):
        return self._record("hbar", *args, **kwargs)

    def errorbar_h(self, *args, **kwargs):
        return self._record("errorbar_h", *args, **kwargs)


@pytest.fixture
def prepared_data():
    df = pd.DataFrame(
        {
            "chrom": [1, 1, 2, 2],
            "pos": [1_000_000, 2_000_000, 1_000_000, 2_000_000],
            "p": [1e-8, 0.01, 1e-6, 0.5],
        }
    )
    return (
        prepare_manhattan_data(df, species="human"),
        prepare_qq_data(df),
    )


def test_renderer_owns_manhattan_panel_policy(prepared_data):
    manhattan_df, _ = prepared_data
    backend = RecordingBackend()

    ManhattanQQRenderer(backend).render_manhattan(
        manhattan_df,
        figsize=(12, 5),
        significance_threshold=5e-8,
        title="Contract Manhattan",
    )

    names = [name for name, _, _ in backend.calls]
    assert names[0] == "create_figure"
    assert "scatter" in names
    assert "axhline" in names
    assert names[-1] == "finalize_layout"
    assert "set_title" in names
    assert "set_xticks" in names


def test_renderer_owns_qq_panel_policy(prepared_data):
    _, qq_df = prepared_data
    backend = RecordingBackend()

    ManhattanQQRenderer(backend).render_qq(
        qq_df,
        figsize=(6, 6),
        show_confidence_band=True,
        show_lambda=True,
        title=None,
    )

    names = [name for name, _, _ in backend.calls]
    assert names[0] == "create_figure"
    assert "fill_between" in names
    assert "line" in names
    assert "scatter" in names
    assert names[-1] == "finalize_layout"


@pytest.mark.parametrize("backend_name", ["matplotlib", "plotly", "bokeh"])
def test_same_prepared_intent_renders_on_each_builtin_backend(
    backend_name, prepared_data
):
    """All built-in adapters accept the same prepared rendering intent."""
    manhattan_df, qq_df = prepared_data
    backend = get_backend(backend_name)
    renderer = ManhattanQQRenderer(backend)
    figures = [
        renderer.render_manhattan(
            manhattan_df,
            figsize=(12, 5),
            significance_threshold=5e-8,
            title="Contract Manhattan",
        ),
        renderer.render_qq(
            qq_df,
            figsize=(6, 6),
            show_confidence_band=True,
            show_lambda=True,
            title=None,
        ),
    ]
    assert all(figure is not None for figure in figures)


def test_miami_region_highlight_is_optional_for_legacy_backends(prepared_data):
    """A pre-existing custom backend need not implement the new capability."""
    manhattan_df, _ = prepared_data
    backend = RecordingBackend()

    figure = MiamiRenderer(backend).render(
        manhattan_df,
        manhattan_df,
        pos_col="pos",
        p_col="p",
        rs_col=None,
        top_threshold=None,
        bottom_threshold=None,
        top_label=None,
        bottom_label=None,
        top_snp_annotations=None,
        bottom_snp_annotations=None,
        highlight_regions=[("1", 1_000_000, 2_000_000)],
        highlight_color="yellow",
        highlight_alpha=0.3,
        figsize=(12, 8),
        title=None,
    )

    assert figure is not None


def test_regional_renderer_skips_heatmap_panel_without_the_capability():
    """A backend that declines SupportsHeatmap loses the panel, not the figure.

    The heatmap is one panel among several in a regional plot, so an adapter
    without it degrades the way the SNP-label and recombination gates do,
    rather than failing the whole render.
    """
    import numpy as np

    from pylocuszoom._regional import RegionalPlotComposer

    backend = RecordingBackend()
    ld_matrix = pd.DataFrame(np.eye(3), index=list("abc"), columns=list("abc"))

    composer = RegionalPlotComposer(backend, genomewide_line=5e-8)

    composer.render_heatmap_panel(
        ax=SimpleNamespace(),
        fig=SimpleNamespace(),
        ld_matrix=ld_matrix,
        x_positions=[1, 2, 3],
        snp_ids=["a", "b", "c"],
        metric="r2",
        lead_snp_id="a",
        start=0,
        end=10,
    )

    assert backend.calls == []


def test_ld_heatmap_renderer_owns_panel_policy():
    """LDHeatmapRenderer drives the heatmap, its scale, ticks, and layout."""
    import numpy as np

    from pylocuszoom._ld_heatmap_renderer import LDHeatmapRenderer

    backend = FullCapabilityBackend()

    LDHeatmapRenderer(backend).render(
        np.eye(3),
        ["rs1", "rs2", "rs3"],
        lead_idx=0,
        highlight_indices=[2],
        metric="dprime",
        figsize=(8.0, 8.0),
        title="Contract Heatmap",
        show_colorbar=True,
    )

    names = [name for name, _, _ in backend.calls]
    assert names[0] == "create_figure"
    assert names[-1] == "finalize_layout"
    assert names.count("add_heatmap") == 1
    # One highlight for the lead, one for each extra SNP.
    assert names.count("highlight_heatmap_snp") == 2
    assert "set_xticks" in names and "set_yticks" in names
    assert "set_title" in names

    colorbar = next(
        kwargs for name, _, kwargs in backend.calls if name == "add_colorbar"
    )
    assert colorbar["label"] == "D'"


def test_ld_heatmap_renderer_skips_the_colorbar_when_not_asked():
    import numpy as np

    from pylocuszoom._ld_heatmap_renderer import LDHeatmapRenderer

    backend = FullCapabilityBackend()

    LDHeatmapRenderer(backend).render(
        np.eye(2),
        ["rs1", "rs2"],
        lead_idx=None,
        highlight_indices=[],
        metric="r2",
        figsize=(8.0, 8.0),
        title=None,
        show_colorbar=False,
    )

    names = [name for name, _, _ in backend.calls]
    assert "add_heatmap" in names
    assert "add_colorbar" not in names
    assert "highlight_heatmap_snp" not in names
    assert "set_title" not in names


def test_stats_renderer_owns_phewas_panel_policy():
    """PheWAS grouping, the significance line, and axis policy live above the seam."""
    from pylocuszoom._stats_renderer import StatsRenderer

    backend = FullCapabilityBackend()
    df = pd.DataFrame(
        {
            "phenotype": ["a", "b", "c"],
            "p_value": [1e-9, 1e-4, 0.3],
            "neglog10p": [9.0, 4.0, 0.52],
            "category": ["x", "x", "y"],
        }
    )

    StatsRenderer(backend).render_phewas(
        df,
        variant_id="rs1",
        phenotype_col="phenotype",
        p_col="p_value",
        category_col="category",
        effect_col=None,
        significance_threshold=5e-8,
        figsize=(10.0, 8.0),
    )

    names = [name for name, _, _ in backend.calls]
    assert names[0] == "create_figure"
    assert names[-1] == "finalize_layout"
    # One scatter per category, not one per phenotype.
    assert names.count("scatter") == 2
    assert names.count("axvline") == 1
    assert "set_yticks" in names


def test_stats_renderer_draws_no_significance_line_for_none():
    from pylocuszoom._stats_renderer import StatsRenderer

    backend = FullCapabilityBackend()
    df = pd.DataFrame(
        {
            "phenotype": ["a"],
            "p_value": [1e-9],
            "neglog10p": [9.0],
            "category": ["x"],
        }
    )

    StatsRenderer(backend).render_phewas(
        df,
        variant_id="rs1",
        phenotype_col="phenotype",
        p_col="p_value",
        category_col="category",
        effect_col=None,
        significance_threshold=None,
        figsize=(10.0, 8.0),
    )

    assert "axvline" not in [name for name, _, _ in backend.calls]


def test_stats_renderer_owns_forest_panel_policy():
    """Forest error bars, markers, null line, and x-padding live above the seam."""
    from pylocuszoom._stats_renderer import StatsRenderer

    backend = FullCapabilityBackend()
    df = pd.DataFrame(
        {
            "study": ["A", "B"],
            "beta": [0.1, 0.4],
            "ci_lower": [0.0, 0.2],
            "ci_upper": [0.2, 0.6],
        }
    )

    StatsRenderer(backend).render_forest(
        df,
        variant_id="rs1",
        study_col="study",
        effect_col="beta",
        ci_lower_col="ci_lower",
        ci_upper_col="ci_upper",
        weight_col=None,
        null_value=0.0,
        effect_label="Beta",
        figsize=(8.0, 6.0),
    )

    names = [name for name, _, _ in backend.calls]
    assert names[0] == "create_figure"
    assert names[-1] == "finalize_layout"
    # One batched errorbar and one batched scatter, not one call per study.
    assert names.count("errorbar_h") == 1
    assert names.count("scatter") == 1
    assert names.count("axvline") == 1

    xlim = next(args for name, args, _ in backend.calls if name == "set_xlim")
    # 10% padding either side of the CI span, which spans 0.0 to 0.6.
    assert xlim[1:] == pytest.approx((-0.06, 0.66))


def test_coloc_renderer_owns_panel_policy():
    """Colocalisation axes, thresholds, legend, and layout live above the seam."""
    from pylocuszoom._coloc_renderer import ColocRenderer

    backend = FullCapabilityBackend()
    merged = pd.DataFrame(
        {
            "neglog10_gwas": [8.0, 3.0],
            "neglog10_eqtl": [6.0, 2.0],
            "color": ["#FF0000", "#0000FF"],
            "rs": ["rs1", "rs2"],
        }
    )

    ColocRenderer(backend).render(
        merged,
        lead_idx=0,
        merged_rs_col="rs",
        ld_col_merged=None,
        gwas_threshold=5e-8,
        eqtl_threshold=1e-5,
        show_correlation=True,
        color_by_effect=False,
        h4_posterior=0.92,
        title="Contract Coloc",
        figsize=(8.0, 8.0),
    )

    names = [name for name, _, _ in backend.calls]
    assert names[0] == "create_figure"
    assert names[-1] == "finalize_layout"
    assert "scatter" in names
    # One threshold line per axis.
    assert names.count("axhline") == 1
    assert names.count("axvline") == 1
    assert "set_xlabel" in names and "set_ylabel" in names
    assert "set_title" in names
    # Correlation and H4 posterior are both annotations.
    assert names.count("add_text") >= 2
