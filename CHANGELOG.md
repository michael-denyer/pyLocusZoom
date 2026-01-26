# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.0] - 2026-01-26

### Added
- Bioconda recipe for conda installation
- `adjustText` moved to default dependencies (was optional)
- **Interactive plotly backend** - use `backend="plotly"` for hover tooltips and pan/zoom
- **Interactive bokeh backend** - use `backend="bokeh"` for dashboard-ready plots

### Changed
- `plot()` and `plot_stacked()` now use backend protocol for all rendering (scatter, line, axes, layout)
- **Gene track now works with all backends** (plotly, bokeh, matplotlib)
- **Recombination overlay now works with all backends** - secondary y-axis with rate line and fill
- **LD legend now works with all backends** - r² color scale (lead SNP highlighted in plot, not legend)
- SNP labels remain matplotlib-only (interactive backends use hover tooltips instead)
- Default `genomewide_threshold` changed from 5e-7 to 5e-8 (standard GWAS significance)
- Gene track strand colors: forward strand now goldenrod (#DAA520), reverse strand light blue (#6BB3FF)
- Gene track directional arrows: black for forward, dark grey for reverse
- Added panel spacing (hspace=0.1) between stacked/fine-mapping panels for visual separation
- Tightened gene track internal spacing for more compact layout

### Fixed
- Bokeh backend `x_range=None` error when creating figures with shared x-axis
- Bokeh backend `legend_label=None` error in scatter plots
- Bokeh backend LD legend not rendering (empty scatter plots don't create legend glyphs)
- Bokeh backend deprecated `FuncTickFormatter` replaced with `CustomJSTickFormatter`
- Bokeh backend deprecated `circle()` method replaced with `scatter(marker=...)`
- Bokeh backend `FIXED_SIZING_MODE` validation warning in column layouts

## [0.2.0] - 2026-01-26

### Added
- Fine-mapping/SuSiE visualization with credible set coloring
- Example plots in `examples/` directory
- Plot generation script for documentation

### Fixed
- Ruff linting and formatting errors
- Bokeh security vulnerability (bumped to >= 3.8.2)
- `plot()` KeyError when `rs_col` column missing with `ld_reference_file` provided
- `plot_stacked()` now validates eQTL DataFrame columns before use
- `plot_stacked()` now validates list lengths for `lead_positions`, `panel_labels`, and `ld_reference_files`
- `calculate_ld()` docstring now documents `ValidationError` for missing PLINK files

### Changed
- Minimum Python version bumped to 3.10 (required by bokeh 3.8.2)
- Renamed species terminology: "dog" → "canine", "cat" → "feline"
- Clarified interactive backend status in README (coming soon)

## [0.1.0] - 2026-01-26

### Added
- Initial release of pyLocusZoom
- Regional association plots with LD coloring
- Gene and exon track visualization
- Recombination rate overlay (canine only)
- Automatic SNP labeling with adjustText
- Species support: Canine (CanFam3.1/CanFam4), Feline (FelCat9), custom
- CanFam4 coordinate liftover via pyliftover
- Stacked plots for multi-GWAS comparison
- eQTL overlay panel support
- PySpark DataFrame support
- Backend infrastructure for matplotlib, plotly, bokeh (matplotlib only active)
- Logging via loguru
- Comprehensive test suite

### Dependencies
- matplotlib >= 3.5.0
- pandas >= 1.4.0
- numpy >= 1.21.0
- loguru >= 0.7.0
- pyliftover >= 0.4
- plotly >= 5.0.0
- bokeh >= 3.8.2
- kaleido >= 0.2.0

[Unreleased]: https://github.com/michael-denyer/pyLocusZoom/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/michael-denyer/pyLocusZoom/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/michael-denyer/pyLocusZoom/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/michael-denyer/pyLocusZoom/releases/tag/v0.1.0
