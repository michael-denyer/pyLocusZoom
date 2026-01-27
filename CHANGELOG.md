# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `plot_phewas()` method for phenome-wide association study plots
- `plot_forest()` method for forest plots (meta-analysis visualization)
- PheWAS category color palette with 12 distinct colors
- Forest plot and PheWAS validation utilities
- Backend methods: `axvline()`, `hbar()`, `errorbar_h()` for new plot types
- Example plots for PheWAS and forest plots

### Changed
- Bumped minimum Plotly version to 5.15.0 (required for multiple legends feature)
- eQTL loaders now output `effect_size` column instead of `effect` for plotter compatibility

### Fixed
- SAIGE loader now prefers SPA-adjusted p-values (`p.value.NA`) over raw p-values when both present
- BED loader now handles BED12 format and files with more than 6 columns
- eQTL panel in `plot_stacked()` now filters by chromosome in addition to position
- Validation errors for non-numeric p-values or positions now show clear "must be numeric" message instead of runtime errors

## [0.5.0] - 2026-01-27

### Added
- Hover tooltips for fine-mapping scatter plots (Plotly/Bokeh backends)
- Hover tooltips for eQTL scatter plots (Plotly/Bokeh backends)
- Interactive HTML example plots for eQTL and fine-mapping (Plotly/Bokeh)
- Comprehensive marker and hover data tests for interactive backends

### Changed
- Plotly/Bokeh backends now hide grid lines for cleaner LocusZoom appearance
- Plotly/Bokeh backends now show black axis lines (matching matplotlib style)
- Plotly/Bokeh gene track panels now hide y-axis (ticks, labels, line, grid)
- Plotly/Bokeh backends now hide minor ticks and zero lines

## [0.4.0] - 2026-01-26

### Added
- **File format loaders** for common GWAS, eQTL, and fine-mapping formats:
  - GWAS: `load_gwas`, `load_plink_assoc`, `load_regenie`, `load_bolt_lmm`, `load_gemma`, `load_saige`, `load_gwas_catalog`
  - eQTL: `load_gtex_eqtl`, `load_eqtl_catalogue`, `load_matrixeqtl`
  - Fine-mapping: `load_susie`, `load_finemap`, `load_caviar`, `load_polyfun`
  - Gene annotations: `load_gtf`, `load_bed`, `load_ensembl_genes`
- Pydantic validation for file loaders with detailed error messages
- `py.typed` marker for PEP 561 type checking support
- Pre-commit configuration for automated linting
- GitHub issue templates for bug reports and feature requests
- Codecov badge in README

### Changed
- eQTL and fine-mapping legends now route through backend protocol (works with all backends)
- Simplified backend code with reduced duplication
- Backend protocol class diagram added to ARCHITECTURE.md

### Fixed
- Additional robustness improvements for edge cases

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

[Unreleased]: https://github.com/michael-denyer/pyLocusZoom/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/michael-denyer/pyLocusZoom/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/michael-denyer/pyLocusZoom/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/michael-denyer/pyLocusZoom/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/michael-denyer/pyLocusZoom/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/michael-denyer/pyLocusZoom/releases/tag/v0.1.0
