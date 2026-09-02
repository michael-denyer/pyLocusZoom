"""Each regional panel type builds itself and renders through the composer."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from pylocuszoom._data import prepare_pvalue_data
from pylocuszoom._regional import (
    AssociationPanel,
    EqtlPanel,
    FinemappingPanel,
    GenePanel,
    HeatmapPanel,
    RegionalFigurePlan,
    RegionalPlotComposer,
)
from pylocuszoom.config import ColumnConfig, DisplayConfig, RegionConfig
from tests.test_rendering_contract import FullCapabilityBackend, RecordingBackend

REGION = RegionConfig(chrom=1, start=1_000_000, end=2_000_000)


def _gwas():
    return pd.DataFrame(
        {
            "rs": ["rs1", "rs2", "rs3"],
            "ps": [1_100_000, 1_500_000, 1_900_000],
            "p_wald": [1e-8, 1e-5, 1e-3],
        }
    )


def _association(**overrides):
    fields = dict(
        data=prepare_pvalue_data(_gwas(), "p_wald"),
        height=4.0,
        columns=ColumnConfig(),
        display=DisplayConfig(snp_labels=False),
        ld_col=None,
        lead_pos=1_500_000,
        recomb_df=None,
    )
    fields.update(overrides)
    return AssociationPanel(**fields)


def _render(panel, backend=None):
    backend = backend or RecordingBackend()
    plan = RegionalFigurePlan(
        chrom=REGION.chrom,
        start=REGION.start,
        end=REGION.end,
        panels=[panel],
        figsize=(8.0, panel.height),
    )
    RegionalPlotComposer(backend, genomewide_line=7.3).render(plan)
    return backend


def _names(backend):
    return [name for name, _, _ in backend.calls]


def test_association_panel_draws_points_line_and_lead():
    names = _names(_render(_association()))

    assert names[0] == "create_figure"
    assert names.count("scatter") == 2
    assert "axhline" in names
    assert names[-1] == "finalize_layout"


def test_association_panel_label_and_ld_legend():
    data = prepare_pvalue_data(_gwas(), "p_wald").assign(R2=[1.0, 0.5, 0.1])
    names = _names(
        _render(
            _association(data=data, ld_col="R2", panel_label="A", add_ld_legend=True)
        )
    )

    assert "add_panel_label" in names
    assert "add_legend" in names


def test_finemapping_panel_from_frame_filters_sorts_and_renders():
    fm = pd.DataFrame(
        {"pos": [1_900_000, 500, 1_100_000], "pip": [0.2, 0.9, 0.8], "cs": [1, 1, 1]}
    )

    panel = FinemappingPanel.from_frame(fm, REGION, cs_col="cs")

    assert list(panel.data["pos"]) == [1_100_000, 1_900_000]
    assert panel.height == 1.5
    backend = _render(panel)
    assert "add_legend" in _names(backend)
    assert [a[1] for n, a, _ in backend.calls if n == "set_ylabel"] == ["PIP"]


def test_eqtl_panel_renders_threshold_line():
    eqtl = pd.DataFrame({"pos": [1_200_000, 5], "p_value": [1e-6, 1e-3]})
    panel = EqtlPanel.from_frame(eqtl, REGION, gene=None, threshold=1e-5)

    assert list(panel.data["pos"]) == [1_200_000]

    names = _names(_render(panel))

    assert names.count("scatter") == 1
    assert "axhline" in names


@pytest.mark.parametrize(
    ("starts", "ends", "expected_height"),
    [
        ([1_100_000], [1_200_000], 1.0),
        ([1_100_000, 1_150_000, 1_160_000], [1_400_000, 1_450_000, 1_460_000], 2.0),
    ],
)
def test_gene_panel_height_grows_with_stacked_rows(starts, ends, expected_height):
    genes = pd.DataFrame(
        {
            "chr": ["1"] * len(starts),
            "start": starts,
            "end": ends,
            "gene_name": [f"G{i}" for i in range(len(starts))],
        }
    )

    panel = GenePanel.from_genes(genes, REGION, exons_df=None)

    assert panel.height == expected_height
    assert len(panel.data) == len(starts)


def test_gene_panel_from_genes_drops_other_chromosomes_and_renders():
    genes = pd.DataFrame(
        {
            "chr": ["1", "2"],
            "start": [1_100_000, 1_100_000],
            "end": [1_200_000, 1_200_000],
            "gene_name": ["G1", "G2"],
        }
    )

    panel = GenePanel.from_genes(genes, REGION, exons_df=None)

    assert list(panel.data["gene_name"]) == ["G1"]
    assert "set_xlim" in _names(_render(panel))


def _ld_matrix(ids):
    return pd.DataFrame(np.eye(len(ids)), index=ids, columns=ids)


def test_heatmap_panel_from_matrix_keeps_region_snps_and_lead():
    ids = ["rs1", "rs2", "rs3", "rs9"]

    panel = HeatmapPanel.from_matrix(
        _ld_matrix(ids),
        ids,
        source=_association(),
        region=REGION,
        height=1.0,
        metric="r2",
    )

    assert panel is not None
    assert panel.snp_ids == ["rs1", "rs2", "rs3"]
    assert panel.x_positions == [1_100_000, 1_500_000, 1_900_000]
    assert panel.lead_snp_id == "rs2"
    assert panel.matrix.shape == (3, 3)


def test_heatmap_panel_from_matrix_is_none_without_overlap():
    ids = ["rs8", "rs9"]

    panel = HeatmapPanel.from_matrix(
        _ld_matrix(ids),
        ids,
        source=_association(),
        region=REGION,
        height=1.0,
        metric="r2",
    )

    assert panel is None


def test_heatmap_panel_renders_on_a_capable_backend():
    ids = ["rs1", "rs2", "rs3"]
    panel = HeatmapPanel.from_matrix(
        _ld_matrix(ids),
        ids,
        source=_association(),
        region=REGION,
        height=1.0,
        metric="r2",
    )

    names = _names(_render(panel, FullCapabilityBackend()))

    assert names.count("add_heatmap") == 1
    assert names.count("highlight_heatmap_snp") == 1


def test_render_panel_rejects_unknown_panel_types():
    composer = RegionalPlotComposer(RecordingBackend(), genomewide_line=7.3)
    plan = RegionalFigurePlan(chrom=1, start=1, end=2, panels=[], figsize=(1.0, 1.0))

    with pytest.raises(TypeError, match="No renderer for str"):
        composer.render_panel("not a panel", SimpleNamespace(), SimpleNamespace(), plan)
