# ADR 0005: Move Heatmap and Bar-Chart Drawing to Capability Protocols

- Status: accepted
- Date: 2026-07-25
- Target: 2.1 (additive change to the `PlotBackend` extension contract)

## Context

ADR-0004 settled how optional backend features are negotiated: an
`@runtime_checkable` protocol per capability, checked with `isinstance`, so a
custom backend declares support by implementing the methods and declines by
omitting them. It applied that rule to SNP labels, the secondary axis, and
region highlights, but left five drawing methods on the required `PlotBackend`
protocol that are not general primitives:

- `add_heatmap`, `add_colorbar`, and `highlight_heatmap_snp` are called only by
  `LDHeatmapRenderer` and `RegionalPlotComposer.render_heatmap_panel`, both of
  which exist to draw an LD matrix.
- `hbar` and `errorbar_h` are called only by `StatsRenderer.render_forest`.
  `hbar` has no in-tree caller at all.

A backend that never draws an LD heatmap or a forest plot still had to supply
all five to satisfy the protocol. Plotly's `add_colorbar` was a `pass` stub for
exactly this reason, and that stub silently discarded the caller's
`show_colorbar` and `label` arguments.

Two further pieces of the heatmap contract were duplicated rather than owned.
The bounds check and the row-then-column cell walk inside `highlight_heatmap_snp`
were written out separately in all three adapters, and the hover
column-name-to-number-format heuristic in `hover.py` had no production caller
because Plotly and Bokeh each re-derived it inline in `scatter`.

## Decision

Extract `SupportsHeatmap` (`add_heatmap`, `add_colorbar`,
`highlight_heatmap_snp`) and `SupportsBarCharts` (`hbar`, `errorbar_h`) as
`@runtime_checkable` protocols, and remove those methods from `PlotBackend`.
Method signatures are unchanged, so the three bundled backends satisfy both
protocols without modification.

Each caller negotiates the capability according to what it can still produce
without it:

- `RegionalPlotComposer.render_heatmap_panel` skips the panel with a debug log.
  The heatmap is one panel among several, so the rest of the regional figure
  still renders. This matches the existing SNP-label and recombination gates.
- `LDHeatmapRenderer` raises `TypeError` from its constructor, and
  `StatsRenderer.render_forest` raises `TypeError` before drawing. There the
  capability is the entire figure, so a silent skip would return a blank panel.
  The message names the backend class and the missing protocol.

Alongside, geometry and formatting that three adapters were each deriving move
to a single owner above the seam: `composition.heatmap_highlight_cells` returns
the cells a SNP highlight covers, and `hover.plotly_hovertemplate` /
`hover.bokeh_tooltips` build each backend's tooltip spec.

`hbar` is retained despite having no in-tree caller. It is a working, tested
public primitive, and 2.0 has shipped; removing it would break any downstream
caller for the sake of about ninety lines.

## Considered options

- Split each adapter's heatmap methods into a mixin module so the Plotly and
  Bokeh files drop under 1000 lines: rejected. The methods still have to resolve
  on the backend instance for the protocol to be satisfied, so a mixin buys a
  line count at the cost of one more hop between a reader's question and the
  code that answers it. The defect the review found was a fat required protocol,
  not file length.
- Fold `add_colorbar` into `add_heatmap` as a `show_colorbar` argument: rejected.
  It removes a protocol method, but changes a primitive's signature for every
  custom backend after 2.0 shipped. Making Plotly's `add_colorbar` enable the
  trace's own scale fixes the same two bugs inside the existing contract.
- Leave the five methods on `PlotBackend` and document them as optional in
  prose: rejected, because it is the boolean-flag failure mode ADR-0004 already
  rejected. A method a backend must define and may leave empty is not optional.

## Consequences

`PlotBackend` drops from 33 required methods to 28, and a custom backend that
implements none of the five optional capabilities still renders every regional,
Manhattan, Miami, colocalisation, and PheWAS plot. Adding a capability protocol
is now a repeatable move rather than a one-off, and the ADR-0004 rule covers
every optional method on the seam.

Two Plotly defects the `pass` stub concealed are fixed as part of the change:
`show_colorbar=False` is honoured, and a D' heatmap's scale is titled `D'`
rather than `R²`. This is a visible behaviour change for Plotly LD heatmaps and
is recorded in the changelog.
