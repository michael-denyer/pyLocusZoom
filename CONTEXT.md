# pyLocusZoom — Domain & Design Context

Shared vocabulary for architecture work on pyLocusZoom. Names here are the ones
to use in design discussion, reviews, and ADRs so the language stays consistent.

Architecture terms follow the deep-module vocabulary — **module, interface,
depth, seam, adapter, port, leverage, locality**. Domain terms name the GWAS
visualization concepts the library renders. Settled decisions live in
[docs/adr/](docs/adr/); do not re-litigate them.

## Domain language

- **Regional association plot** — the core figure: `-log10(p)` against genomic
  position for a locus, optionally stacked with other panels.
- **Panel** — one horizontal band of a regional figure (association scatter, gene
  track, eQTL, fine-mapping, LD heatmap). A figure is an ordered list of panels.
- **Lead SNP** — the anchor variant of a locus; LD and colouring are computed
  relative to it.
- **LD (linkage disequilibrium)** — `r^2` of each variant to the lead SNP, binned
  into colour bands.
- **Recombination map** — cM/Mb rate track drawn on a secondary axis; canine maps
  auto-download and are lifted over between genome builds.
- **Liftover** — mapping coordinates between genome builds (CanFam3.1 → CanFam4).
- **eQTL / colocalization / fine-mapping (credible set) / PheWAS / forest /
  Manhattan / QQ / Miami** — the plot families, each with its own data intake and
  semantic renderer.

## Current seams

- **`PlotBackend` primitive protocol** (`backends/base.py`) — the external
  rendering seam. Structural (`typing.Protocol`); matplotlib, Plotly, and Bokeh
  are adapters. Documented public extension point via `@register_backend`.
- **Semantic renderers** — own panel composition, axes, legends, and layout,
  translating backend-neutral figure intent into primitive calls
  (`ManhattanQQRenderer`, `MiamiRenderer`, `ColocRenderer`, `StatsRenderer`,
  `LDHeatmapRenderer`, `RegionalPlotComposer`).
- **`RegionalPlotComposer` + `RegionalFigurePlan`** (`_regional.py`) — the deep
  renderer for regional figures; accepts a typed, ordered panel plan and owns
  figure creation, axis assignment, and dispatch.
- **`DataFrameValidator`** (`validation.py`) — the fluent validation engine
  (`require_columns`/`require_numeric`/`require_range`/`require_not_null`/
  `require_ci_ordering`). Per-family `error_class`.
- **`prepare_pvalue_data`** (`_data.py`) — the p-value intake owner: filters
  invalid p-values and produces `neglog10p`, with `allow_zero` distinguishing the
  Manhattan convention (0 allowed, clipped) from the strict eQTL `(0, 1]` domain.
- **Download port** (`recombination._download_with_progress`) — the single
  network seam every recomb/liftover fetch routes through.
- **Ensembl `cache_dir` port** (`ensembl.py`) — injected cache directory; HTTP
  behind `_make_ensembl_request`.

### Load-bearing intake distinctions (not defects)

- **GWAS two-tier strictness** — load-time validation (`schemas.py`) is strict
  (numeric, `(0, 1]`, non-positive positions); plot-time validation
  (`utils.validate_gwas_df`) checks only columns + non-empty, deferring p-value
  policy to `prepare_pvalue_data` so the Manhattan `p=0` convention survives.
- **eQTL `gene` requirement** — required at load (multi-gene GTEx files keyed by
  gene) but optional at plot (frames may already be gene-filtered). Same
  intentional two-tier shape as GWAS.
- **Nulls are a load-time failure, not a plot-time one** — every `schemas.py`
  validator rejects nulls in the numeric columns it range-checks, so a malformed
  file fails at load with the offending column named. Plot-time intake stays
  lenient and filters nulls with a warning (`prepare_pvalue_data`,
  `prepare_eqtl_for_plotting`), because a frame reaching the plotter may have
  been assembled by the caller rather than loaded from a file.

## Deepenings from the design grill

Approved in a design grill and recorded in ADR-0004. The composition module, the
`CoordinateLifter` port, and the capability protocols have shipped; the
validation engine/spec split is still outstanding.

- **Composition module** (`backends/composition.py`, ADR-0004) — pure functions
  that drive backend primitives above the seam:
  - **`LegendEntry`** — backend-neutral legend spec (`label`, `color`, `marker`,
    `edgecolor`); drained through a single neutral `add_legend(entries)` primitive.
    Retires the five semantic legend methods.
  - **`SecondaryAxis` handle** — opaque per-backend return of `create_twin_axis`;
    all `*_secondary` primitives take it, letting one
    `render_recombination_overlay` replace the 61-line-×3 duplication.
  - **Capability protocols** — `SupportsSNPLabels`, `SupportsSecondaryAxis`
    (joining `SupportsRegionHighlight`) as the single capability-negotiation
    mechanism; `supports_hover` stays a boolean (rendering-quality flag).
- **`CoordinateLifter` port** (`_liftover.py`) — liftover seam with a production
  `pyliftover` adapter and an in-memory test adapter; pure `liftover_positions`.
  Plus pure `ensure_recomb_header` for header detection.
- **Validation engine/spec split** (outstanding) — `schemas.py` becomes per-family specs
  expressed on `DataFrameValidator` (deleting the parallel hand-rolled
  implementation), with a single `require_pvalue` domain helper so the strict
  `(0, 1]` policy has one owner shared with `prepare_pvalue_data`.
