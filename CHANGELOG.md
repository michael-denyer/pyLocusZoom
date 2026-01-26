# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Example plots in `examples/` directory
- Plot generation script for documentation

### Fixed
- Ruff linting and formatting errors
- Bokeh security vulnerability (bumped to >= 3.8.2)

### Changed
- Minimum Python version bumped to 3.10 (required by bokeh 3.8.2)
- Clarified interactive backend status in README (coming soon)

## [0.1.0] - 2026-01-26

### Added
- Initial release of pyLocusZoom
- Regional association plots with LD coloring
- Gene and exon track visualization
- Recombination rate overlay (dog only)
- Automatic SNP labeling with adjustText
- Species support: Dog (CanFam3.1/CanFam4), Cat (FelCat9), custom
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

[Unreleased]: https://github.com/michael-denyer/pyLocusZoom/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/michael-denyer/pyLocusZoom/releases/tag/v0.1.0
