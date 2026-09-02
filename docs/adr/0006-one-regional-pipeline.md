# ADR 0006: One Regional Pipeline, Panels That Build Themselves

- Status: accepted
- Date: 2026-09-02
- Supersedes: the second consequence of ADR-0002 ("future regional modes
  can submit prepared panels without duplicating axis indices or layout
  branches"), which the code did not deliver

## Context

ADR-0002 introduced typed panel dataclasses and `RegionalPlotComposer` so
that single and stacked regional plots would share one composition path.
By 2.1 they did not. `plot()` and `plot_stacked()` were two copies of the
same pipeline, 194 and 289 lines, and the copies had drifted: eQTL and
fine-mapping panels existed only on `plot_stacked()`, the constructor's
`auto_genes` flag was honoured only by `plot()`, and one of the two heatmap
paths read the raw frame where the other read the prepared one. Every
regional feature was a two-place edit.

The composer's `render` was a five-arm `isinstance` chain whose job was to
unpack each frozen panel back into keyword arguments for a `render_*_panel`
method with a matching signature, so each panel type had two parameter
lists to keep in step. Panel construction (gene-row layout, heatmap
coordinate mapping, fine-mapping and eQTL preparation) lived in the
plotter, not in the panel types, and the eQTL block was a private copy of
the exported `prepare_eqtl_for_plotting` with two different rules.

`StackedPlotConfig` copied `PlotConfig`'s four nested fields instead of
extending it, and both carried a `from_kwargs` that restated every default
the nested models already declare.

## Decision

- `plot()` and `plot_stacked()` build their config, resolve their per-panel
  lists, and call one private `LocusZoomPlotter._render_regional`. The only
  mode difference is height policy, expressed as two arguments: the
  association panel height and the floor on the figure height. `plot()` is
  the single-element case. Both methods therefore accept every optional
  panel: `plot()` gains the eQTL and fine-mapping arguments, and
  `plot_stacked()` gains `auto_genes`.
- Each panel type owns its construction: `FinemappingPanel.from_frame`,
  `EqtlPanel.from_frame`, `GenePanel.from_genes`, and
  `HeatmapPanel.from_matrix` validate and prepare raw caller input through
  the exported helpers, so the composer trusts a prepared frame.
  `AssociationPanel` carries `ColumnConfig` and `DisplayConfig` rather than
  six loose fields copied out of them.
- `RegionalPlotComposer.render_panel` is a `functools.singledispatchmethod`
  taking the panel object. One registered method per panel type replaces
  the `isinstance` chain and the keyword-unpack signatures. The heatmap
  capability gate from ADR-0005 lives inside the `HeatmapPanel` arm.
- `StackedPlotConfig` extends `PlotConfig`. One generic `from_kwargs` routes
  each keyword to the nested model that declares it; defaults live in the
  `Field()` declarations and the public signatures only.
- The eQTL panel goes through `prepare_eqtl_for_plotting`. A zero eQTL
  p-value is dropped with the existing range warning, and `eqtl_gene` on a
  frame with no `gene` column raises `EQTLValidationError`. A warning that
  then draws an unfiltered panel is a silent failure dressed as a plot.
- `RegionConfig.start` is `>= 1`, matching `lead_pos`; coordinates are
  1-based throughout.

## Considered options

- A `RegionalRequest` dataclass wrapping the twenty arguments
  `_render_regional` takes: rejected for now. It would have two builders and
  one reader, and wrapping a keyword list in a frozen dataclass later is
  mechanical. It can be added when a third caller appears.
- A `render(composer, ax, plan)` method on each panel dataclass instead of
  dispatch in the composer: rejected. It is five pass-throughs that couple
  value types to the composer, and a reader looking for "how is an eQTL
  panel drawn" lands on a method that forwards to the composer anyway.
- A `dict[type, method]` dispatch table: an equal alternative to
  `singledispatchmethod`. The stdlib form was chosen because it has no
  table to maintain and a grep for `EqtlPanel` lands on its renderer.
- Keeping the inline eQTL rules and changing the exported helper to match:
  rejected. ADR-0002 named eQTL's strict `(0, 1]` domain as the reason
  `allow_zero` exists, and `tests/test_eqtl.py` pins both helper behaviours.

## Consequences

`plotter.py` drops from 754 lines to 596 and `plot_stacked()` from radon F
(42) to C (11); `config.py` drops from 436 lines to 314. A regional feature
is a one-place edit, and the answer to "what does an eQTL panel need" is
the `EqtlPanel` class.

Two public tightenings, both recorded in the changelog: `start=0` is
rejected, and `eqtl_gene` naming an absent column raises instead of
warning. Rendered output for every existing call is unchanged; the example
plots regenerate byte-identical after id normalisation.

Height policy stays dual (sum for `plot()`, floored at `figsize` for
`plot_stacked()`). The answer to "why sum here and clamp there" is existing
output, and unifying it would be a visible change for its own sake.

`tests/test_regional_plan.py` renders every panel type through the
contract tests' recording backend, so the panel model is tested directly
rather than only through matplotlib.
