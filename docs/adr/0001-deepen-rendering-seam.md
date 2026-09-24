# Deepen the rendering seam

- Status: superseded. [ADR-0004](0004-complete-rendering-seam-and-capability-protocols.md)
  dropped the legacy primitive contract for custom backends in 2.0, with no
  shim, and [ADR-0007](0007-one-figure-plan.md) deleted the semantic
  renderers: every family is now a `FigurePlan` of panels that draw themselves.

Keep one external rendering seam that accepts prepared, backend-neutral figure intent. Semantic renderers now cover Manhattan/QQ, Miami, PheWAS/forest, colocalization, LD heatmaps, and regional association panels. They own panel composition, axes, legends, layout, and capability negotiation, while Matplotlib, Plotly, and Bokeh adapters handle native rendering mechanics. Legacy custom backends remain supported through the existing primitive contract so the documented extension point remains usable during migration.

## Considered options

- Keep primitive drawing calls in plotters: rejected because backend mechanics and capability branches remain spread across callers.
- Add a second public semantic seam above the current protocol: rejected because it adds indirection and splits the extension contract.
- Break custom backends immediately: rejected because custom registration is a documented extension point.
