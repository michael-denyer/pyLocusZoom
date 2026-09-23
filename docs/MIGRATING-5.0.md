# Migrating to 5.0

5.0 makes three kinds of breaking change: inputs that 4.x tolerated now raise,
several functions take different arguments or are gone, and the `PlotBackend`
protocol that custom backends implement changes shape. This page lists what to
change. The [CHANGELOG](../CHANGELOG.md) has the reason for each change.

Most call-site changes are mechanical, and the migration script rewrites them
in `.py`, `.ipynb` and `.md` files:

```bash
uv run --with libcst python scripts/migrate_to_config_models.py PATH...
```

It moves flat `PanelInputs` fields into `EqtlInput`, `FinemappingInput` and
`LDHeatmapInput`, and flat `plot_coloc` options into `ColocConfig`. Run
`ruff check --fix` and `ruff format` on the files afterwards.

## Plot calls

- **`plot_coloc` takes a `ColocConfig`.** The column, lead, colouring,
  annotation and figure-size options move into `config=ColocConfig(...)`, which
  `pylocuszoom` now exports. `gwas_threshold`, `eqtl_threshold` and `title`
  stay per call, and `ColocConfig` no longer has threshold fields.

  ```python
  # 4.x
  plotter.plot_coloc(gwas_df, eqtl_df, ld_col="ld_r2", lead_snp="rs1", gwas_threshold=1e-6)

  # 5.0
  from pylocuszoom import ColocConfig

  plotter.plot_coloc(
      gwas_df,
      eqtl_df,
      config=ColocConfig(ld_col="ld_r2", lead_snp="rs1"),
      gwas_threshold=1e-6,
  )
  ```

- **Options are keyword-only on every plotter.** `plot_phewas(df, variant_id,
  ...)`, `plot_forest(df, variant_id, ...)` and `plot_ld_heatmap(ld_matrix,
  snp_ids, ...)` take every later argument by name, as the other plotters
  already did.
- **`PanelInputs` nests one model per optional panel.** `eqtl_df`, `eqtl_gene`
  and `eqtl_threshold` become `eqtl=EqtlInput(data=, gene=, threshold=)`;
  `finemapping_df` and `finemapping_cs_col` become
  `finemapping=FinemappingInput(data=, cs_col=)`; the four `ld_heatmap_*`
  fields become `ld_heatmap=LDHeatmapInput(matrix=, snp_ids=, height=,
  metric=)`. The migration script does this.
- **`LocusZoomPlotter(log_level=...)` is removed.** Logging is off until you
  call `enable_logging()`: replace `log_level="DEBUG"` with
  `enable_logging("DEBUG")`.
- **The GWAS loaders lose `pos_col`, `p_col` and `rs_col`,** and
  `load_gwas(**kwargs)` is removed. The loaders emit `chr`, `pos`, `p_value`
  and `rs`; rename afterwards if you need other names.
- **The LD functions need `species`.** `calculate_ld`,
  `calculate_pairwise_ld` and the two command builders defaulted to canine.
  Pass `species="canine"` to keep the old command; `species=None` means
  PLINK's default (human) chromosome set.

## Inputs that now raise

Every input error is a `pylocuszoom.ValidationError`, which subclasses
`ValueError`. That includes a config model's rejected value (it was pydantic's
`ValidationError`) and a keyword a config model does not declare.

- A column you name must exist. `LDConfig(ld_col=)`, a non-default
  `ColumnConfig(rs_col=)`, `plot_phewas(effect_col=)`, a non-default
  `category_col`, `plot_forest(weight_col=)` and a non-default
  `FinemappingInput(cs_col=)` used to switch the feature off silently.
- A regional frame needs a chromosome column, `chr` by default. Name another
  with `ColumnConfig(chrom_col=...)`, or pass `chrom_col=None` for a frame that
  holds only the region's chromosome. eQTL and fine-mapping frames follow the
  same rule through `EqtlInput(chrom_col=)` and `FinemappingInput(chrom_col=)`;
  `load_susie` output has no chromosome column.
- LD from a reference fileset needs a SNP id column, and every stacked panel
  that computes it needs a lead (`lead_positions=[...]` or
  `LDConfig(lead_pos=...)`).
- A lead position must lie inside the region, `start` to `end` inclusive. Omit
  it to auto-detect the strongest in-region SNP.
- `exons_df` is validated, the LD heatmap metric must be `"r2"` or `"dprime"`,
  and `plot_manhattan_qq_stacked` checks the length of `panel_labels`.
- Frames in the pre-4.0 `ps`/`p_wald` names are no longer read as `pos` and
  `p_value`: name them with `ColumnConfig` or `GenomeWideConfig`.

## Removed and moved names

| 4.x | 5.0 |
|-----|-----|
| `validate_gwas_df` and the other five `validate_*_df` wrappers | `validation.check(df, schemas.gwas_plot_spec(...))` and the other `*_plot_spec` builders |
| `RecombStatus`, `RecombResult`, `recomb_for_region` | `get_recombination_rate_for_region`, which raises a `PyLocusZoomError` subclass |
| `FileNotFoundError` from `load_recombination_map` | `RecombinationMapNotFound`, a `FileNotFoundError` subclass |
| bare `OSError` for an unwritable map cache | `DataDownloadError` |
| `FileNotFoundError` when PLINK is missing | `PlinkError` |
| `utils.ASSEMBLY_SYNONYMS`, `reference_genes.ucsc_genome_for_build` and the other build tables | `genome_build.resolve_build(...)` and its `GenomeBuild` fields |
| `liftover_region(..., species=)` | `liftover_region(..., build=)` |
| `recombination.download_liftover_chain`, `liftover_recombination_map` and the chain-path helpers | `LiftoverConfig(chain_path=...).resolve()` and `liftover_region` |
| `manhattan.prepare_manhattan_frames(..., chrom_col=, pos_col=, p_col=)` | `prepare_genomewide_frames(dfs, GenomeWideConfig(...), species=...)` |
| `ensembl.ENSEMBL_*` and `ucsc.UCSC_*` retry constants | nothing: they restated the shared retry defaults (30 s timeout, 3 attempts, 1 s delay) |
| `pylocuszoom.backends.Mappable` | nothing: `add_heatmap` returns `None` |

## Custom backends

A backend registered with `@register_backend` implements the `PlotBackend`
protocol in `pylocuszoom/backends/base.py`. 5.0 changes these members. Each
change removes something no built-in caller used, or folds two calls into one
([ADR-0011](adr/0011-protocol-diet-and-one-panel-body.md)).

| Member | 4.x | 5.0 |
|--------|-----|-----|
| `supports_hover` | property returning a bool | removed; hover data may reach any backend, and a static backend ignores it |
| `add_legend` | `(ax, entries, loc="upper left", title=None)` | `(ax, entries, title=None)`; draw the legend in the panel's upper-right corner, the only place any caller asked for |
| `add_text` | `(ax, x, y, text, fontsize=10, ha="center", va="bottom", rotation=0, color="black")` | `rotation` removed |
| `add_heatmap` | `(ax, data, x_coords, y_coords, cmap_colors, vmin=0.0, vmax=1.0) -> Mappable` | `(ax, data, x_coords, y_coords, cmap_colors, vmin=0.0, vmax=1.0, colorbar_label=None) -> None`; when `colorbar_label` is a string, draw a vertical colour scale beside the panel titled with it |
| `add_colorbar` | `(ax, mappable, label="R²", orientation="vertical")` | removed; see `add_heatmap` |
| `scatter` | `hover_data: Optional[pd.DataFrame]` | `hover_data: Optional[HoverData]`: `hover_data.frame` holds the columns under their display names and `hover_data.roles` one `HoverRole` per column. Format a column by its role (`hover.plotly_hovertemplate` and `hover.bokeh_tooltips` show how), not by its name |
| `scatter`, `set_title`, `set_suptitle` | `alpha` and `fontweight` left out when unset, so a backend written before 4.1 or 4.2 still drew | always passed: `alpha` may be `None` (the backend's default opacity) and `fontweight` is `"bold"` or `"normal"` |
| `SupportsSNPLabels.add_snp_labels` | `(ax, df, pos_col, neglog10p_col, rs_col, label_top_n, adjust=True, lead_pos=None, region_span=None) -> List[Any]` | `(ax, df, pos_col, neglog10p_col, rs_col, label_top_n) -> None`; the caller has already chosen the rows eligible for a label, so label the `label_top_n` strongest |

`HoverData` and `HoverRole` live in `pylocuszoom.backends.hover`. In the same
module, `HoverDataBuilder.build_dataframe` is now `HoverDataBuilder.build`
and returns a `HoverData`, and `plotly_hovertemplate` and `bokeh_tooltips`
take one.
