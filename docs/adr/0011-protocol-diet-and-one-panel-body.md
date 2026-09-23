# ADR 0011: Protocol Diet, Typed Hover Roles and One Body per Panel

- Status: accepted
- Date: 2026-09-23
- Target: 5.0 (breaking change to the `PlotBackend` extension contract)
- Supersedes: the ADR-0005 option it rejected ("Fold `add_colorbar` into
  `add_heatmap`") and with it `add_colorbar` and the `Mappable` handle; the
  ADR-0007 decision that "the Manhattan and QQ specs forward to
  `render_manhattan_panel` and `render_qq_panel`", with the `manhattan_spec`
  and `categorical_spec` builders; and the ADR-0004 decision that
  "`supports_hover` remains a boolean property"

## Context

The protocol carried parameters no caller ever passed, and the backends paid
for each one. Every `add_legend` call passed `loc="upper right"`, so the
protocol default `"upper left"` was dead and plotly and bokeh each kept a
matplotlib-vocabulary location table with a silent fallback. No caller passed
`add_colorbar(orientation=)`, `add_text(rotation=)` or
`add_snp_labels(lead_pos=, region_span=)`, and `adjust` was always `True`.
`supports_hover` had one caller, the Manhattan panel, which skipped building a
hover frame on matplotlib; the regional panels built theirs unconditionally,
so one question had two policies. 4.1 and 4.2 added two shims,
`scatter_alpha` and `title_weight`, that left a keyword out so a backend
written before it existed kept working. ADR-0004 had already ruled that no
compatibility shim is added to the protocol.

Two primitives needed a handshake. `add_heatmap` returned an opaque
`Mappable` that only `add_colorbar` consumed, the one return value crossing
the seam, and plotly handed back `fig.data[-1]` because `add_trace` stores a
copy. ADR-0005 kept the two calls separate because folding them changed a
primitive's signature after 2.0 shipped. The protocol has broken three times
since, so that argument no longer holds. Both callers also repeated one draw
sequence: the heatmap, the colour bar labelled from the metric, and the
outline rectangles.

Hover formats were guessed from display names. `HoverDataBuilder` knew which
column was the SNP id, the position, the p-value and the r², then emitted a
plain frame; `hover.py` sniffed names for `"pos"` and `"ld"`, and plotly
assumed column 0 was the SNP id. A fine-mapping or eQTL frame has no SNP id,
so plotly showed the position bold, unlabelled and unformatted where bokeh
showed it labelled and comma-grouped.

The Manhattan panel was declared twice: `manhattan_spec` restated thirteen of
`ManhattanPanelSpec`'s fields with their defaults, and a test compared the
two lists. `ManhattanPanelSpec.draw` and `QQPanelSpec.draw` forwarded to free
functions nothing else called, where every other panel's `draw` is its body.

## Decision

- Delete the never-passed parameters: `add_legend(loc)`, drawn in the upper
  right always; `add_text(rotation)`; `add_snp_labels(lead_pos, region_span,
  adjust)`, whose selection the panel already does above the seam. Delete
  `supports_hover`, and the `scatter_alpha` and `title_weight` shims:
  `alpha=None` and `fontweight` are always passed.
- `add_heatmap` takes `colorbar_label: Optional[str]` and returns `None`;
  `None` draws no colour scale. `add_colorbar` and `Mappable` are deleted, so
  every drawing primitive returns `None`. `composition.draw_ld_heatmap` owns
  the heatmap, colour bar and outline sequence both LD heatmap panels draw.
- `scatter(hover_data=)` takes a `HoverData`: the display-named frame plus
  one `HoverRole` per column, assigned by `HoverDataBuilder` from the
  `HoverConfig` field the column came from. `plotly_hovertemplate` and
  `bokeh_tooltips` format by role, and both label every field.
- `ManhattanPanelSpec` takes the `PreparedManhattan` value, which carries the
  x and group columns its preparation created. The builders and the drift
  test go, and `ManhattanPanelSpec.draw` and `QQPanelSpec.draw` hold the
  drawing code themselves.

## Considered options

- Keep each dead parameter as documented vocabulary: rejected. A custom
  backend must implement and test a keyword the library never sends, and a
  silent fallback for an unknown `loc` hides a typo in the one caller that
  could send it.
- Keep `supports_hover` as a performance flag: rejected. The hover frame it
  saves costs about 5 ms over a million rows against hundreds of
  milliseconds for the scatter itself, and matplotlib ignores `hover_data`
  whether it is built or not.
- Fold the colour bar as a `show_colorbar: bool` beside a separate `label`:
  rejected. A label without a bar is not a state any caller needs, and one
  optional argument cannot express it.
- Keep the plain hover DataFrame and read roles from `DataFrame.attrs`:
  rejected. `attrs` does not survive every pandas operation, so the role
  would be lost without an error.

## Consequences

Custom backends break once, in one release, and the migration note lists
each signature change. The protocol loses `supports_hover` and
`add_colorbar`. The plotly and bokeh legend location tables and their
fallbacks are gone, and plotly's heatmap sets its colour bar when it adds the
trace, so no trace is looked up after the fact.

Plotly and bokeh render the same hover fields with the same formats, so
plotly fine-mapping and eQTL exports gain a labelled, formatted position, and
the SNP line of an association tooltip gains its label. Where a panel draws
does not change.
