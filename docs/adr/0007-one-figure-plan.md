# ADR 0007: One Figure Plan for Every Family

- Status: accepted
- Date: 2026-09-02
- Supersedes: the third decision of ADR-0006 (`RegionalPlotComposer.render_panel`
  as a `singledispatchmethod`) and the option it rejected ("a `render` method
  on each panel dataclass")

## Context

The same idea, a figure is an ordered list of prepared panels drawn onto
backend axes, was implemented three ways under three names. Regional plots
had `RegionalFigurePlan`, five panel values and `RegionalPlotComposer`, whose
five dispatch arms each forwarded to a `draw_*` function in the module the
panels lived in, threading a `fig` no arm read. Manhattan and QQ had panel
specs and free draw functions but no figure value, so the figure shape was
hard-coded in six `ManhattanQQRenderer` methods. PheWAS and forest had a
renderer class whose two methods shared nothing but the backend field, and
Miami, colocalisation and the standalone LD heatmap each had a frozen request
value and a `render_*` function that opened with `create_figure` and closed
with `finalize_layout`.

That put `create_figure` at ten call sites, `create_figure_grid` at two and
`finalize_layout` at fifteen, across seven modules, with the suptitle branch
written out three times. Figure-level policy had no owner, so a change to
layout fractions or a backend needing a different figure setup was a
seven-module edit, and a new plot family meant reading three existing ones
and picking a template.

`ManhattanQQRenderer.render_manhattan` and `render_categorical` were the same
eleven lines around a different spec; `_plotter_utils.add_significance_line`,
the PheWAS renderer and the eQTL panel each drew their own dashed threshold
line; and `LocusZoomPlotter` threaded its threshold through the composer while
every other family resolved it through `UNSET` into an `Optional[float]`.

## Decision

- `_figure.py` holds one frozen `FigurePlan` and one `render_figure`. The plan
  is the panels in row-major order plus the figure-level policy the families
  vary: `figsize`, `height_ratios`, `n_cols` and `width_ratios`, `sharex`, an
  `xlabel` for the bottom axis, `mb_xaxis` for a genomic x axis read in
  megabases, `highlights` spanning every panel, a `suptitle`, and the `top`
  and `hspace` layout fractions. `render_figure`
  creates the figure or grid, calls `draw(backend, ax)` on each panel, applies
  the figure-level policy, and finalizes. It is the only caller of
  `create_figure`, `create_figure_grid`, `set_suptitle` and `finalize_layout`
  outside `backends/`.
- Every panel value has `draw(backend, ax)` and nothing else reads the plan.
  The five regional panels carry their region and draw the way their
  `draw_*` functions did; the Manhattan and QQ specs forward to
  `render_manhattan_panel` and `render_qq_panel`.
- A plotter validates, prepares, builds a plan, and returns
  `render_figure(self._backend, plan)`. Plan builders issue no backend calls.
- The renderer classes go. `RegionalPlotComposer` is deleted, its threshold
  moving onto `AssociationPanel` as the same `Optional[float]` every other
  family draws from. `ManhattanQQRenderer`'s six methods become plan builders,
  with `render_manhattan` and `render_categorical` one builder over
  `manhattan_spec` or `categorical_spec`. `StatsRenderer`'s two methods become
  `PhewasPanel` and `ForestPanel`. `ColocRequest` and `LDHeatmapRequest` are
  `ColocPanel` and `LDHeatmapPanel`, panels that draw themselves. Miami stays
  a request value, `MiamiRequest`, because it describes a two-panel figure:
  `miami_plan` turns it into two `MiamiPanel`s and the plan's highlights.
- One `add_significance_line` draws every threshold line, horizontal or
  vertical, with the colour and opacity the site had.

## Considered options

- `draw(backend, ax, plan)`, as the reviews sketched: rejected. Once each
  regional panel carries its region, no panel reads the plan, and a parameter
  no implementation uses is a lie in the signature.
- A callable on the plan for steps between drawing and finalizing (the
  regional x-axis label and megabase formatting, the Miami highlights):
  rejected. A frozen value carrying a closure is dispatch by another name.
  The label became the plan's `xlabel` and the highlights its `highlights`,
  both generic.
- Each regional panel formatting its own x axis in megabases inside `draw`,
  beside the `set_xlim` it already issues, instead of a plan-level flag:
  tried and reverted. Plotly writes layout keys in insertion order, so
  formatting the bottom axis before labelling it reorders four regional
  exports, and labelling it from the panel needs an `x_label` on all five
  panel types. `mb_xaxis` is a figure-level property in the same sense as
  `sharex`, and `render_figure` labels before it formats, with the reason in
  a comment at that line.
- Keeping `RegionalFigurePlan` beside `FigurePlan`: rejected. Two names for
  one thing is the plurality this record exists to end.
- Deleting `MiamiRequest` and building the plan inline in `plot_miami`:
  rejected. The contract test wants a value it can construct without a
  plotter, and a fifteen-field request is the value.

## Consequences

`create_figure`, `create_figure_grid`, `set_suptitle` and `finalize_layout`
have one caller above the seam. `_regional.py` and `_rendering.py` are
deleted; `_stats_renderer.py`, `_coloc_renderer.py` and
`_ld_heatmap_renderer.py` hold panels and `_miami_renderer.py` holds the
request, the panel and the plan builder. The answer to "how do I add a plot
type" is one sentence: write a panel value with `draw`, and a plotter method
that returns a `FigurePlan`.

Rendered output for every existing call is unchanged; the example plots
regenerate identically after id normalisation at every commit of the
migration. `test_rendering_contract.py` keeps every assertion it made on the
primitive call sequence, re-pointed at the new names, and
`tests/test_figure.py` pins what `render_figure` sends to the backend for a
column plan, a grid plan, and the figure-level policy.

`plot()` and `plot_stacked()` still take no threshold argument; the panel
now has the field a later change will feed.
