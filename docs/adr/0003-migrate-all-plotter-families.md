# ADR 0003: Migrate Plotter Families to Semantic Renderers

- Status: accepted
- Date: 2026-07-20

## Context

The first rendering seam pilot covered Manhattan and QQ figures, but the
remaining plotter families still embedded panel composition and backend
dispatch in their public classes. Miami also needed backend-specific region
highlighting logic.

## Decision

Route every standalone plotter family through a semantic renderer:

- Manhattan/QQ: `ManhattanQQRenderer`
- Miami, PheWAS/forest, colocalization, and LD heatmap use focused private
  modules: `_miami_renderer.py`, `_stats_renderer.py`,
  `_coloc_renderer.py`, and `_ld_heatmap_renderer.py`.
- Regional association panels and regional heatmap layout:
  `RegionalPlotComposer`

Add the optional `SupportsRegionHighlight` capability protocol so Miami does
not branch on backend names. Existing registered backends that do not implement
the capability keep their previous behavior.

## Consequences

Plotter classes retain validation and domain preparation, while renderers own
panel composition, axes, labels, legends, and layout. Existing primitive
backends remain compatible, and the semantic renderer classes provide focused
contract surfaces for future backend implementations.
