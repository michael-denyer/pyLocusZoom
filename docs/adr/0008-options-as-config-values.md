# ADR 0008: Plot Methods Take Their Options as Config Values

- Status: accepted
- Date: 2026-09-02
- Supersedes: the fourth decision of ADR-0006 ("One generic `from_kwargs`
  routes each keyword to the nested model that declares it")

## Context

`plot()` and `plot_stacked()` declared 26 and 29 keyword parameters, 24 of
them identical. Each was restated in the `PlotConfig.from_kwargs(...)` call
that immediately followed the signature and again in the docstring, and
twelve of them were already fields of `PanelInputs`. Adding one option cost
six edits in `plotter.py` against one in `config.py`, and `plot_stacked`'s
docstring had given up documenting the shared arguments twice. `from_kwargs`
undid the flattening the signature had just performed by inspecting
`model_fields`, and in doing so constructed each nested model and then let
the parent revalidate it, so every validator ran twice.

The Manhattan family had the opposite problem. Its five methods and
`plot_miami` restated `chrom_col`, `pos_col`, `p_col` and `custom_chrom_order`
by hand, imported no config model and validated nothing at the boundary, so
`ARCHITECTURE.md`'s claim that inputs are validated at the API boundary was
true for one family.

Two shapes were considered for the regional methods.

- Pass the models as values: `plot(df, *, chrom, start, end, columns=
  ColumnConfig(), display=DisplayConfig(), ld=LDConfig(), panels=
  PanelInputs())`. Eight parameters, each option declared once, editor
  completion and static typing on every field, and a model built once in a
  notebook reused across calls. The common call gains one import line and
  one constructor per concern it uses.
- Keep the flat call and accept `**kwargs` routed by `from_kwargs`, typed
  through an `Unpack[TypedDict]` derived from the models. The user's call
  does not change, but pydantic cannot emit the `TypedDict`, so it would be
  a second declaration to keep in step by hand or generated code, static
  checkers would see nothing without it, and the double validation and the
  untyped router stay.

## Decision

- `plot()` and `plot_stacked()` take the region (`chrom`, `start`, `end`)
  and the four models as values. `plot_stacked()` keeps the per-panel lists
  (`lead_positions`, `panel_labels`, `ld_reference_files`) as flat
  parameters: they are parallel to `gwas_dfs`, not options of one panel.
  Both compose a `PlotConfig` or `StackedPlotConfig`, where the cross-model
  rules live; callers never build the composite. `from_kwargs` is deleted.
- The four models are exported from `pylocuszoom` and stay frozen.
- `auto_genes` is a `DisplayConfig` field, `None` inheriting the plotter's
  constructor setting, so `plot()` honours it too. `label_top_n` is
  `Optional[int]`, `None` meaning the method's own default (5 on `plot()`, 3
  per panel on `plot_stacked()`), which keeps the two defaults the old
  signatures carried without a second model.
- Both regional methods take `significance_threshold` through the `UNSET`
  sentinel, matching every other family: omit it to inherit the plotter's
  `genomewide_threshold`, pass `None` to draw no line.
- The genome-wide families take one `GenomeWideConfig` (`chrom_col`,
  `pos_col`, `p_col`, `custom_chrom_order`) and route every frame through
  `manhattan.prepare_genomewide_frames`, which validates each against those
  names before layout. `figsize`, `title`, the threshold and the per-method
  options stay per call: their defaults differ per figure shape.
  `rs_col` stays a `plot_miami` parameter because only the Miami hover reads
  it. `plot_qq` reads `config.p_col` only.
- Every option after the frame is keyword-only across both families.
- No compatibility shim, following ADR-0004 and the 3.0.0 precedent.
  `scripts/migrate_to_config_models.py` rewrites `.py`, `.ipynb` and `.md`
  callers in place with libcst.

## Consequences

- One declaration per option. `grep -c ld_heatmap_metric plotter.py` is 1, the
  line that reads it, from 6.
- `_render_regional` no longer reads a resolved `auto_genes` argument; it
  resolves the display value against the plotter itself.
- A caller migrating from 3.x writes one import line and wraps the options
  in the model that owns them; the migration script does it mechanically.
- Step 6 of the remediation (canonical column vocabulary) will touch
  `ColumnConfig` and `GenomeWideConfig` defaults, not the plot signatures.
