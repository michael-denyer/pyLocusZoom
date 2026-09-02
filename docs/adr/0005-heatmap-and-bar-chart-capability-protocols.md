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
- `require_heatmap_backend` raises `TypeError`, which `LDHeatmapPlotter` calls
  when it is constructed, and
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

## Addendum: a value type is not a mixin

The rejection of a mixin above rules out one shape of decomposition, not all of
them. A mixin keeps the extracted code on the backend instance, so a reader
tracing a protocol method still resolves through the class to find it, and the
method count on the seam is unchanged. That hop is the cost the rejection names.

Extracting a value type and pure functions has no such hop.
`backends/plotly_layout.py` holds `_Panel`, the panel handle itself, which
owns Plotly's subplot axis naming, plus `configure_legend`,
`panel_y`, `x_range`, and `secondary_axis_key`. None of them touch the backend
instance, so a reader following `set_xlim` reads one call to a named function
rather than a method that might be overridden further up an inheritance chain.
`backends/_coerce.py` is the same move for the coercions out of matplotlib's
vocabulary that Plotly and Bokeh each need.

The rule this leaves: extract along the axis of what the code *is*, not to hit a
line count. Arithmetic over a figure, coercion between vocabularies, and
geometry a caller could compute belong outside the adapter. Anything that draws,
and anything the protocol names, stays on the backend class.

## Addendum: geometry a caller can compute does not belong on the protocol

`SupportsHeatmap` shipped with three methods. `add_heatmap` took `x_coords` and
`y_coords`; `highlight_heatmap_snp` took a SNP index and a SNP count. One grid,
two calls, and only one of them knew where the cells were. Each adapter closed
the gap its own way, so the lead-SNP outline drew at matrix indices on an axis
of base pairs whenever the regional composer supplied genomic coordinates, and
Plotly's shapes carried no axis reference and bound to the first subplot.

The method is removed. `composition.heatmap_highlight_rects(snp_idx, x_coords,
y_coords)` returns outline rectangles in data coordinates, and both renderers
draw them through the `add_rectangle` primitive every backend already
implements. `add_rectangle` gains `facecolor=None` for an unfilled outline.

This is the addendum rule above applied to a protocol method rather than to a
file: geometry a caller can compute is not drawing, so it moves out of the
adapter. What stays on `SupportsHeatmap` is the pair of calls only a rendering
library can make. Batching is the cost. Bokeh drew all outline cells in one
`rect()` call and now draws one renderer per cell, at most `2n - 1` renderers
for `n` SNPs on the panel.

## Addendum: `hbar` is removed and the protocol is renamed

The decision above kept `hbar` on the grounds that it was a working, tested
public primitive. That reasoning weighed a hypothetical downstream caller
against ninety lines and did not count the cost to anyone writing a backend:
a required member of an opt-in protocol has to be implemented correctly for a
figure this library never draws. `StatsRenderer.render_forest` draws with
`errorbar_h`, `scatter` and `axvline`; nothing in the tree ever called `hbar`.

`hbar` leaves the protocol and all three adapters, and `SupportsBarCharts`
becomes `SupportsErrorBars`, which now declares the one method it has. A
third-party backend that implements `hbar` keeps working; the method is simply
no longer part of any contract.

## Addendum: four of the five protocols are folded back in

The split above was correct in principle and one level too fine in practice.
All three shipped backends implemented `SupportsHeatmap`, `SupportsErrorBars`,
`SupportsSecondaryAxis` and `SupportsRegionHighlight`, so every `isinstance`
gate on them guarded a branch no backend could reach. Worse, the unreachable
branches had grown three different policies for the same question: the regional
heatmap panel logged and skipped, the Miami highlight returned silently, and the
forest plot and the standalone heatmap each raised a hand-written `TypeError`.

The four are folded back into `PlotBackend`. `require_heatmap_backend`, the two
bespoke `TypeError` messages and the `_SecondaryAxisBackend` intersection shim
in `composition.py` go with them, and `render_recombination_overlay` takes a
plain `PlotBackend`.

`SupportsSNPLabels` stays optional, and it is the one that earns it: it needs
adjustText, which has no plotly or bokeh equivalent, so the declining branch is
reachable by two of the three backends the library ships. The rule ADR-0004 set
still holds; the test it now has to pass is whether any shipped backend actually
declines.
