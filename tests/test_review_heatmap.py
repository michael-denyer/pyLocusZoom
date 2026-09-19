"""Native heatmap cells use the same coordinate boundaries as SNP outlines."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from pylocuszoom import DisplayConfig, LocusZoomPlotter, PanelInputs
from pylocuszoom.backends import get_backend


@pytest.mark.parametrize("backend_name", ["matplotlib", "bokeh"])
@pytest.mark.parametrize(
    "coordinates, boundaries",
    [
        ([100, 200, 1000], [(50, 150), (150, 600), (600, 1400)]),
        ([0, 1, 2], [(-0.5, 0.5), (0.5, 1.5), (1.5, 2.5)]),
        ([100], [(99.5, 100.5)]),
    ],
)
def test_native_cell_boundaries(backend_name, coordinates, boundaries):
    backend = get_backend(backend_name)
    fig, axes = backend.create_figure([1], (5, 4))
    count = len(coordinates)
    result = backend.add_heatmap(
        axes[0], np.eye(count), coordinates, coordinates, ["white", "red"]
    )
    if backend_name == "matplotlib":
        if axes[0].images:
            xmin, xmax, _, _ = result.get_extent()
            edges = np.linspace(xmin, xmax, count + 1)
        else:
            edges = result.get_coordinates()[0, :, 0]
        actual = list(zip(edges[:-1], edges[1:]))
        plt.close(fig)
    else:
        data = axes[0].renderers[0].data_source.data
        actual = [
            (x - w / 2, x + w / 2) for x, w in zip(data["x"][:count], data["w"][:count])
        ]
    np.testing.assert_allclose(actual, boundaries)


def test_regional_heatmap_sorts_coordinates_and_matrix_together():
    frame = pd.DataFrame(
        {"pos": [100, 200, 1000], "p_value": [0.1, 0.01, 0.001], "rs": ["a", "b", "c"]}
    )
    matrix = pd.DataFrame([[0.9, 0.1, 0.2], [0.1, 0.8, 0.3], [0.2, 0.3, 0.7]])
    fig = LocusZoomPlotter(species=None, backend="bokeh", log_level=None).plot(
        frame,
        chrom=1,
        start=1,
        end=1500,
        display=DisplayConfig(show_recombination=False, snp_labels=False),
        panels=PanelInputs(ld_heatmap_df=matrix, ld_heatmap_snp_ids=["c", "a", "b"]),
    )
    cells = fig.children[-1].renderers[0].data_source.data
    bounds = [(x - w / 2, x + w / 2) for x, w in zip(cells["x"], cells["w"])]
    assert bounds == [
        (50, 150),
        (50, 150),
        (150, 600),
        (50, 150),
        (150, 600),
        (600, 1400),
    ]
    assert cells["value"] == [0.8, 0.3, 0.7, 0.1, 0.2, 0.9]


def test_regional_heatmap_rejects_duplicate_genomic_coordinates():
    frame = pd.DataFrame({"pos": [150, 150], "p_value": [0.1, 0.01], "rs": ["a", "b"]})
    with pytest.raises(ValueError, match="distinct genomic positions"):
        LocusZoomPlotter(species=None, log_level=None).plot(
            frame,
            chrom=1,
            start=100,
            end=200,
            display=DisplayConfig(show_recombination=False, snp_labels=False),
            panels=PanelInputs(
                ld_heatmap_df=pd.DataFrame(np.eye(2)), ld_heatmap_snp_ids=["a", "b"]
            ),
        )


def test_heatmap_highlights_selected_variant_at_duplicate_source_position():
    frame = pd.DataFrame(
        {
            "pos": [150, 150, 250],
            "p_value": [0.1, 1e-8, 0.001],
            "rs": ["weak", "strong", "other"],
        }
    )
    fig = LocusZoomPlotter(species=None, log_level=None).plot(
        frame,
        chrom=1,
        start=100,
        end=300,
        display=DisplayConfig(show_recombination=False, snp_labels=False),
        panels=PanelInputs(
            ld_heatmap_df=pd.DataFrame(np.eye(2)),
            ld_heatmap_snp_ids=["other", "strong"],
        ),
    )
    outlines = fig.axes[1].patches
    assert [patch.get_xy() for patch in outlines] == [(100, -0.5), (100, 0.5)]
    assert [patch.get_width() for patch in outlines] == [100, 100]


@pytest.mark.parametrize("with_genes", [False, True])
@pytest.mark.parametrize("with_recombination", [False, True])
def test_regional_colorbar_preserves_genomic_display_alignment(
    sample_genes_df, sample_recomb_df, with_genes, with_recombination
):
    frame = pd.DataFrame(
        {
            "pos": [1_100_000, 1_200_000, 1_900_000],
            "p_value": [0.1, 0.01, 0.001],
            "rs": ["a", "b", "c"],
        }
    )
    fig = LocusZoomPlotter(species=None, log_level=None).plot(
        frame,
        chrom=1,
        start=1_000_000,
        end=2_000_000,
        display=DisplayConfig(show_recombination=with_recombination, snp_labels=False),
        panels=PanelInputs(
            ld_heatmap_df=pd.DataFrame(np.eye(3)),
            ld_heatmap_snp_ids=["a", "b", "c"],
            genes_df=sample_genes_df if with_genes else None,
            recomb_df=sample_recomb_df if with_recombination else None,
        ),
    )
    fig.canvas.draw()
    association = fig.axes[0]
    for panel in association.get_shared_x_axes().get_siblings(association):
        np.testing.assert_allclose(
            panel.transData.transform([[1_100_000, 0], [1_900_000, 0]])[:, 0],
            association.transData.transform([[1_100_000, 0], [1_900_000, 0]])[:, 0],
        )
    colorbar = next(axis for axis in fig.axes if axis.get_ylabel() == "R²")
    assert colorbar.get_position().x0 > association.get_position().x1
    assert colorbar.get_tightbbox(fig.canvas.get_renderer()).x1 <= fig.bbox.x1
