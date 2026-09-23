# ADR 0010: Named Means Required; Intake Raises Only PyLocusZoomError

- Status: accepted
- Date: 2026-09-23
- Amends: ADR-0002 (the p-value intake gains an explicit invalid-row policy),
  ADR-0008 (`PanelInputs` nests one model per optional panel) and ADR-0009
  (a regional frame's chromosome is a resolved column role, not a fallback)

## Context

`exceptions.py` promised that `except PyLocusZoomError` catches every library
error. Five of seven probed input errors bypassed it: plain `ValueError`s from
ad-hoc checks in `manhattan.py`, `qq.py`, `_data.py` and `plotter.py`, and
pydantic's own `ValidationError` from the config models, a class that shares
the name of `pylocuszoom.ValidationError` and none of its hierarchy.

Options the caller named were dropped when the frame did not carry them. A
misspelt `ld_col` drew every point grey; a misspelt PheWAS `effect_col`,
forest `weight_col` or `finemapping_cs_col` drew without the feature; a
missing `rs` column skipped LD computation with a log line, while the LD
heatmap raised for the same condition. `plot_manhattan_qq_stacked` accepted
one label for two frames. `exons_df` was never validated, so a misnamed
column surfaced as a `KeyError`. Any `ld_heatmap_metric` other than `"r2"`
labelled the colour bar D′. ADR-0006 already calls a warning followed by an
unfiltered panel "a silent failure dressed as a plot"; these were the same
failure without the warning.

The regional path had no chromosome role. `filter_by_region` filtered by
position alone whenever the frame lacked a `chr` column, so a whole-genome
frame whose chromosome column was named `chrom` plotted rows from every
chromosome at the requested positions. The genome-wide path compared raw
strings, so `chr1, chr2, chr10` laid out lexicographically.

Whether an invalid p-value was dropped or rejected was decided implicitly, by
whether a family's column contract carried `pvalue=`. A `Family` × `Tier`
registry and seven one-line `validate_*_df` wrappers sat between the callers
and `validation.check`, and three different exception types reported a
missing column. The pre-4.0 column aliases and the loader output-column
overrides were scheduled for removal in 5.0.

## Decision

- **Intake raises only `PyLocusZoomError`.** Every input error raises
  `pylocuszoom.ValidationError` or a subclass. It still subclasses
  `ValueError`, so `except ValueError` keeps working. The config models share
  one base whose constructor re-raises pydantic's error as
  `pylocuszoom.ValidationError`, naming each failing field. `to_pandas`
  rejects an unsupported frame the same way.
- **Named means required.** A column the caller names must be in the frame,
  or the boundary raises a `ValidationError` naming the parameter and the
  column. Nothing downstream probes the frame to compensate. Three optional
  columns have canonical defaults (`rs`, the fine-mapping `cs` and the PheWAS
  `category`): left at the default, a frame without the column draws without
  the feature; any other name is required. `validation.resolve_column` is
  the one owner of that rule. LD from a reference fileset needs SNP ids, so
  it requires the `rs_col` column. The same holds for option names: a config
  model rejects a field it does not declare instead of ignoring it.
- **The chromosome is a column role.** `ColumnConfig.chrom_col` defaults to
  `"chr"`, and a frame without it raises. `chrom_col=None` is the explicit
  opt-in to position-only selection, for a frame already scoped to one
  chromosome; the eQTL and fine-mapping inputs take the same field.
  `utils.normalize_chrom_series` is the one spelling of chromosome identity:
  a `chr` prefix is stripped, in any case, and integral floats read as
  integers. Genome-wide preparation projects the chromosome through it with
  the other roles, so a `chr`-prefixed frame lays out in species order.
- **One p-value policy table.** `_data.P_VALUE_POLICY` states, per family,
  whether an exact zero is a valid p-value and whether an invalid row (null,
  non-numeric, out of range) is dropped or rejected. The table keeps each
  family's 4.x outcome:

  | Family | Zero | Invalid rows |
  |--------|------|--------------|
  | Regional association, genome-wide (Manhattan, Miami, categorical) | valid | dropped |
  | QQ, regional eQTL | invalid | dropped |
  | PheWAS, colocalization | invalid | rejected |
  | Loaders | invalid | rejected |

  The plot families that draw many variants drop a few bad rows and warn;
  the families that draw one row per phenotype or merge two sources reject
  them, because a dropped row there is a missing result. `prepare_pvalue_data`
  applies the table for every plot family, and `validation.check` applies the
  loader row for every loader.
- **Column contracts are values.** `schemas.py` holds one `ColumnSpec`
  constant per fixed contract and one builder per contract over caller column
  names. Callers run them with `validation.check`. The registry and the
  wrappers are deleted.
- **One model per optional panel.** `PanelInputs` takes `eqtl=EqtlInput(...)`,
  `finemapping=FinemappingInput(...)` and `ld_heatmap=LDHeatmapInput(...)`,
  so an eQTL gene without an eQTL frame cannot be expressed. Thresholds and
  heights are bounded, `ld_heatmap.metric` is `"r2"` or `"dprime"`, and every
  frame field collects a Spark frame through `to_pandas`.
- **The LD lead rule holds per panel.** `StackedPlotConfig` resolves each
  panel's `LDConfig` from the broadcast and the per-panel lists, and every
  panel that computes LD from a fileset needs a lead position, as `plot()`
  already required.
- **No deprecation shim survives 5.0.** The pre-4.0 column aliases and the
  GWAS loader output-column parameters are deleted, following ADR-0004 and
  ADR-0008. `scripts/migrate_to_config_models.py` rewrites flat `PanelInputs`
  fields into the nested models.

## Consequences

- A typo in a named option fails at the call with the parameter's name,
  instead of drawing a plausible figure without the feature.
- A regional frame without a `chr` column needs `ColumnConfig(chrom_col=None)`
  to plot. The CHANGELOG lists this with every other input that now raises.
- A caller catching pydantic's `ValidationError` from a config constructor
  catches `pylocuszoom.ValidationError` instead.
- `validate_gwas_df`, `validate_genes_df`, `validate_phewas_df`,
  `validate_forest_df`, `validate_eqtl_df` and `validate_finemapping_df` leave
  the toolbox tier. The plot methods validate their own input.
- The loader zero rule and the Manhattan zero rule still differ; the table
  now says so in one place instead of two mechanisms implying it.
