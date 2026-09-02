# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Gene tracks for CanFam3.1, CanFam4 and FelCat9, sourced from UCSC.** Ensembl retired all three and its archive REST hosts redirect to a help page, so they have no Ensembl source at any URL. `LocusZoomPlotter` with `auto_genes=True` now fetches them from UCSC's `ncbiRefSeq` track, so the canine default build of CanFam3.1 and the feline default of FelCat9 both produce a gene track in the coordinates the rest of the plot is in. Any other build still comes from Ensembl. The whole policy is the `UCSC_BUILDS` table in `reference_genes.py`.
- **`get_genes_for_build()`**, the build-aware entry point that picks the source, alongside `fetch_genes_from_ucsc()`, `fetch_exons_from_ucsc()`, `get_genes_for_region_ucsc()` and `clear_ucsc_cache()`. `UCSCAPIError` mirrors `EnsemblAPIError`, and both subclass the new `ReferenceAPIError`, so a service failure stays distinguishable from a region that genuinely has no genes whichever source answered.
- **`plot()` accepts `eqtl_df`, `eqtl_gene`, `eqtl_threshold`, `finemapping_df` and `finemapping_cs_col`, and `plot_stacked()` accepts `auto_genes`.** The two methods were separate copies of one pipeline that had drifted: eQTL and fine-mapping panels existed only on `plot_stacked()`, and the constructor's `auto_genes` was honoured only by `plot()`. Both now run through one private pipeline, so every optional panel is available from either entry point. `plot_stacked(auto_genes=None)` inherits the constructor setting.
- **Gene and exon frames from Ensembl carry an `assembly` column** recording the assembly each row was served on. It is written to the disk cache, so the mismatch warning also fires on a cache hit in a later session rather than only on the original fetch. The column is optional on frames you supply yourself and is ignored when plotting.

### Fixed

- **Ensembl gene tracks no longer land silently in the wrong coordinate system.** Ensembl serves exactly one reference assembly per species and ignores a `coord_system_version` asking for any other, returning HTTP 200 either way. A canine plot built on CanFam3.1 therefore drew genes in ROS_Cfam_1.0 coordinates, which puts `ATP9B` on chr1 at 938,796 instead of 1,136,865, a shift of about 198 kb, with no error and nothing in the figure to show it. Ensembl 116 is the final release on the legacy platform, so CanFam3.1 will not return; feline moved from FelCat9 to F.catus_Fca126_mat1.0 on the same terms. `fetch_genes_from_ensembl`, `fetch_exons_from_ensembl` and `get_genes_for_region` now take `genome_build` and warn with a `UserWarning` naming both assemblies when they disagree, and `LocusZoomPlotter` passes its own `genome_build` through when `auto_genes` is enabled.
- **An interrupted download no longer leaves a truncated file the cache trusts.** The liftover chain and the recombination tarball were streamed straight to their destination, so a connection dropped mid-download left a partial `canFam3ToCanFam4.over.chain.gz` that passed the exists check on every later CanFam4 plot and made pyliftover fail with no route to recovery short of finding the cache directory. Downloads now stream to a `.part` sibling and are renamed into place only once the stream completes; a failure removes the partial file.
- **Bokeh text annotations honour every alignment pair.** `add_text` mapped four `(ha, va)` combinations and fell back to centre/bottom for the rest, so the gene track's "No genes" notice and the colocalisation r/p annotation sat at the wrong anchor on bokeh only. Horizontal and vertical alignment now map independently, as they already did on matplotlib and plotly.
- **`load_gwas` names the file in its unknown-format warning.** The message used a printf-style `%s` under loguru, which formats with braces, so the emitted line contained a literal `%s` and no path.
- **A gene-source outage warns instead of drawing an empty gene track.** With `auto_genes=True`, an Ensembl or UCSC failure came back as an empty frame and was logged at debug as "No genes found in region". `LocusZoomPlotter.plot()` now emits a `UserWarning` naming the region and the cause, and draws the plot without the gene track.

### Changed

- **Loaders that could not map their columns now raise instead of warning (breaking).** `load_gtex_eqtl`, `load_eqtl_catalogue`, `load_susie`, `load_finemap`, `load_polyfun`, `load_ensembl_genes` and `load_bed` used to log a warning, skip validation and return a frame missing the columns their docstring promises, so the failure surfaced later as a plot-time error naming the output columns rather than the source ones. All seven now raise `LoaderValidationError` at load time, where the candidate list and the file's actual columns are still known. `load_matrixeqtl` gains validation of the columns it does carry. `load_caviar` and `load_matrixeqtl` still load without a position column, which neither format provides.
- **`utils.filter_by_region` raises `ValidationError`, not `KeyError` (breaking).** A missing position column leaked pandas' vocabulary through a boundary that otherwise speaks `ValidationError`, and was reachable from `plot_stacked`. Code catching `KeyError` around it must catch `ValidationError`; both remain catchable as `PyLocusZoomError`.
- **`utils.validate_gwas_df` and `utils.validate_genes_df` moved to `validation.py`.** They were a second hand-rolled validation tier beside the `ColumnSpec` engine that already expressed them. `utils.validate_dataframe`, `schemas.validate_file_path`, `colors.get_eqtl_color_palette` and `finemapping.calculate_credible_set_coverage` are removed; none was exported or called anywhere in the library.
- **Manhattan and QQ share the one p-value intake.** `manhattan._filter_invalid_pvalues` and the inline QQ filter are gone; both now call `_data.prepare_pvalue_data`, which gained `out_col` and `on_empty` so each keeps the column name and the empty-input error it raised before. A p-value column of strings is now coerced and filtered rather than raising `TypeError`.
- **`BackendError` is removed (breaking).** It was exported and documented but nothing in the library ever raised it; code that imports it must drop the import.
- **`ReferenceAPIError` now subclasses `DataDownloadError`, not `ValidationError` (breaking).** An Ensembl or UCSC outage is a download failure, not an input validation failure, so it is no longer catchable as `ValueError`. `EnsemblAPIError` and `UCSCAPIError` follow their parent; `except PyLocusZoomError` still catches all three.
- **Recombination download failures raise `DataDownloadError`.** `ensure_recomb_header`, `download_canine_recombination_maps` and `download_liftover_chain` raised bare `RuntimeError` for a corrupt archive header, an incomplete map set and a missing map file, and the chain download re-raised `requests.HTTPError`, so `except PyLocusZoomError` missed every failure in the one subsystem that touches the network. All of them now raise `DataDownloadError`, which is still a `RuntimeError`.
- **`ensure_recomb_maps` has one failure policy.** A network failure used to degrade silently to no overlay while a corrupt archive crashed `plot()`. Every download failure now emits a `UserWarning` with the cause and the plot renders without the recombination overlay.

- **`add_heatmap` no longer takes `mask_upper`, a breaking change for custom backends.** Both in-tree callers passed `mask_upper=True`, so which half of a symmetric LD matrix to draw was never a backend decision. `composition.lower_triangle()` now masks the matrix above the seam and each backend draws whatever it is handed, which also removes three separate implementations of the same mask. A custom backend's `add_heatmap` must drop the parameter; a masked cell arrives masked, and an unmasked `NaN` still means missing data.
- **Internal: the interactive backends share their coercions out of matplotlib's vocabulary.** Figure sizing, marker area to diameter, and scalar broadcast were written twice, once per backend, and now live in `backends/_coerce.py`. Plotly's subplot geometry moves to `backends/plotly_layout.py` as a value type plus pure functions, which puts every file under `backends/` below 1000 lines. Bokeh's 63-line colour-palette reimplementation is replaced by matplotlib's colormap, which is byte-identical at the 256 colours the library asks for. No rendered output changes.
- **The `PlotBackend` contract shrank from 28 required methods to 24, a breaking change for custom backends.** `hide_spines` is gone: it was matplotlib styling pushed through the seam, stubbed as `pass` on plotly and bokeh and called 17 times by the renderers. Matplotlib now hides the top and right spines in `create_figure` and `create_figure_grid`, which is what plotly's template and bokeh's figure styling already did implicitly. `save`, `show` and `close` are gone too; they had no caller outside the backends, because the plotter layer returns the native figure and callers export with its own methods. A custom backend that still defines any of the four is unaffected, since the protocol is structural.
- **`SupportsSNPLabels` lost `adjust_snp_labels`, and `add_snp_labels` lost `genes_df` and `chrom`.** The two parameters travelled from `AssociationPanel` through the composer and the protocol into `labels.add_snp_labels`, which deleted them on its first line, so a custom backend had to accept two arguments nothing read. `adjust_snp_labels` advertised a deferred-adjust path no caller took. The module-level `labels.adjust_snp_labels` function is unchanged.
- **Two figures render differently.** The standalone matplotlib LD heatmap loses the top and right border of its box; `LDHeatmapRenderer` never hid a spine, so it was the one figure drawing a full frame, and plotly and bokeh already render it without one. `hide_yaxis` on matplotlib now also hides the left spine, matching its own docstring ("ticks, labels, and line") and the other two backends; the composed LD heatmap panel uses it in place of `set_yticks(ax, [], [])`, so that panel's y-axis line is hidden on every backend.
- **The two gene sources share one fetch-and-cache orchestration.** `ensembl.py` and `ucsc.py` held the same algorithm written twice around different parsers, and `reference_genes.py` sat on top as a two-way `if`. A `GeneSource` value now carries what differs between them and `get_genes_for_build()` holds the single copy of the cache lookup, the fetch, the error translation and the rule that a failed fetch is never cached. Every public name keeps its signature. The one behaviour change is that a UCSC gene track with `include_exons=True` costs one request rather than two: `ncbiRefSeq` rows already carry `exonStarts`/`exonEnds`, so genes and exons come out of the same payload.
- **The recombination downloads retry.** `download_liftover_chain()` and `download_canine_recombination_maps()` had their own HTTP transport with no retry, so the 50 MB tarball the first canine plot depends on got one attempt while a 5 KB JSON payload got three. Both now stream through the shared transport in `_http.py`, which retries a dropped connection and a 429 or 503 with doubling backoff. A 404 still fails on the first attempt.
- **Internal helpers moved to where their callers are.** `validate_plink_files()` moved from `utils.py` to `ld.py`, its only caller; `_normalize_build()` is deleted in favour of `utils.assembly_token()`, which does the same job for every build name the library recognises; and the empty `reference_data/` package is removed, along with the `ARCHITECTURE.md` claim that recombination maps were written into it. Maps have always gone to the platform cache directory.

- **The eQTL panel goes through `prepare_eqtl_for_plotting`, the exported intake, instead of a private copy with different rules (breaking).** Two behaviours change. An eQTL row with `p_value == 0` is dropped with the existing out-of-range warning, as the helper always did; the inline copy drew it. `eqtl_gene` naming a gene when the frame has no `gene` column raises `EQTLValidationError`; the inline copy logged a warning and drew every eQTL unfiltered, which is a silent failure dressed as a plot.
- **`start=0` is rejected by `plot()` and `plot_stacked()` (breaking).** Genomic coordinates are 1-based, which `lead_pos` already enforced; `RegionConfig.start` now carries the same `>= 1` bound instead of `>= 0`.
- **One Manhattan panel, three callers (internal only; no public API change).** The panel policy of scatter loop, significance line, limits, ticks, labels, title and corner label was written three times: once in `ManhattanQQRenderer`, once again for the categorical variant, and twice more inline in `MiamiRenderer`, which also carried its own copy of the chromosome-tick derivation. Changing tick size or spine policy meant a three-file edit, and the code nowhere said that a Miami plot is two Manhattan panels with the lower one mirrored. A frozen `ManhattanPanelSpec` now names what those copies actually varied on and `render_manhattan_panel` draws any of them, so Miami is two specs and the fontsize that used to be inferred from another fontsize is gone. Hover data is built at one gated site for the whole family rather than at a primitive plus a caller deciding separately. Every plotter signature and every rendered figure is unchanged.

- **The Ensembl gene cache key includes the genome build.** Two builds of the same region no longer share an entry. Every entry written before this release is orphaned rather than reused, because the assembly its coordinates are in cannot be recovered. Call `clear_ensembl_cache()` to reclaim the disk.

## [2.1.1] - 2026-07-25

A patch release of user-visible bug fixes. No API changes, no behaviour changes
for code that was already working.

### Fixed

- **A gene-sparse region no longer breaks on its second plot, and a failed Ensembl fetch no longer hides a region's genes.** Both defects surfaced as `pandas.errors.EmptyDataError` raised from inside pandas, with nothing naming pyLocusZoom or its cache. An empty `DataFrame` carries no columns, so it serialised to a one-byte CSV that `pd.read_csv` cannot parse back: caching empty results exists to spare the API on regions that genuinely have no genes, and it had never once round-tripped. Separately, a fetch that failed was cached as if it were an empty region, so a single dropped connection or Ensembl 503 permanently replaced that region's genes with "none" until the user wiped the cache by hand. Empty frames now carry their schema's columns, and `get_genes_for_region` declines to cache a fetch it could not complete. The corrupt-cache guard now also names `EmptyDataError`, which subclasses `ValueError` rather than `ParserError` and so never fired, recovering caches poisoned by earlier releases without user action.
- **PLINK 2 `--glm` output loads.** `_PLINK_SPEC` passed `comment="#"` to `read_csv`, and a PLINK 2 `--glm` header line starts with `#CHROM`. pandas discarded the header, promoted the first variant to column labels, and `load_plink_assoc` failed with `Missing columns: ['ps', 'p_wald']` listing that variant's own values as the available columns. A file also silently lost one variant. The `"#CHROM"` entry in the chromosome candidates was added to support these files and had never once matched. PLINK writes no comment lines of its own, so the prefix is simply dropped; REGENIE and GTF keep theirs, where `##` metadata lines are real. `load_gwas` now also detects `.glm.` filenames instead of warning and falling through to its PLINK default.
- **The Windows cache fallback returns to `~/AppData/Local`.** Unifying the two hand-rolled cache resolvers in 2.1.0 took the recombination fallback for an unset `%LOCALAPPDATA%` (`~/`) over the Ensembl one (`~/AppData/Local`), which moved both caches to the profile root and stopped matching the value the variable stands in for. `%LOCALAPPDATA%` is set on essentially every Windows install, so this branch is close to unreachable; where it is reached, the conventional location now wins. Windows installs with `%LOCALAPPDATA%` set are unaffected, as are macOS, Linux, and Databricks.
- **`StatsRenderer.render_forest` names both methods a backend is missing.** `SupportsBarCharts` requires `hbar` and `errorbar_h`, but the rejection said only "does not support error bars", so a backend implementing `errorbar_h` (the one method a forest plot actually calls) was refused with a message naming the thing it did support. It now names both, matching `LDHeatmapRenderer`.

### Added

- **`EnsemblAPIError`**, raised when the Ensembl REST API is unreachable or returns an error. It subclasses `ValidationError`, so existing `raise_on_error=True` callers are unaffected; it exists so a service failure stays distinguishable from a region that genuinely has no genes.

## [2.1.0] - 2026-07-25

**Upgrading from 2.0.** Two changes can alter working code, neither of which a
version constraint will catch for you.

1. `MiamiPlotter.backend_name`, `ColocPlotter.backend_name`,
   `LDHeatmapPlotter.backend_name`, and `LDHeatmapPlotter.species` no longer
   exist. Reading one now raises `AttributeError`. They were never documented as
   public attributes and nothing in the library read them; the backend name is
   private on every plotter now, as it already was on `LocusZoomPlotter`.
2. `ManhattanPlotter`, `MiamiPlotter`, and `StatsPlotter` now honour the
   `genomewide_threshold` you pass to their constructor. Before 2.1 they stored
   it and drew at 5e-8 regardless. **If you passed that argument, your
   significance line moves.** That is the fix, but it is a visible change: check
   any figure built by a plotter constructed with a non-default threshold. Code
   that never passed the argument is unaffected, and an explicit per-call
   threshold still wins over the constructor's.

Custom backends need no changes. `PlotBackend` sheds five methods to two new
capability protocols, but a 2.0 backend that implements them satisfies both
structurally, and no signature changed.

### Changed

- **Heatmap and bar-chart drawing are optional backend capabilities** (ADR-0005). `add_heatmap`, `add_colorbar`, and `highlight_heatmap_snp` move to a new `SupportsHeatmap` protocol; `hbar` and `errorbar_h` move to `SupportsBarCharts`. Both are `@runtime_checkable`, matching the `SupportsSNPLabels` / `SupportsSecondaryAxis` / `SupportsRegionHighlight` pattern from 2.0. `PlotBackend` drops from 33 required methods to 28, so a custom backend that draws neither LD heatmaps nor forest plots no longer has to stub five methods. **Method signatures are unchanged**, so a 2.0 backend implementing all five satisfies both new protocols with no edits. A backend that declines the capability is handled per caller: the regional heatmap panel is skipped with a debug log, while `LDHeatmapPlotter` and forest plots raise `TypeError` naming the backend and the missing protocol, because there the capability is the whole figure.
- **Ensembl gene and exon fetching share one overlap helper.** `fetch_genes_from_ensembl` and `fetch_exons_from_ensembl` were about 70% duplicate. They are now thin wrappers over a single `_fetch_overlap_features`, and `_make_ensembl_request` carries one raise-on-failure contract instead of a dual raise-or-return-`None` contract. Public signatures and return shapes are unchanged.
- **Internal: DataFrame validators are declarative `ColumnSpec`s.** The nine hand-written `DataFrameValidator` chains across `schemas.py` (GWAS, eQTL, fine-mapping, genes) and the feature modules (`coloc`, `eqtl`, `phewas`, `forest`, `finemapping`) are replaced by a frozen `ColumnSpec` (with `RangeRule`) run through a single `check(df, spec)` engine in `validation.py`. Each validator keeps its own dataset name, error class, rule set, and rule order, so accepted/rejected frames, exception types, and error messages are byte-for-byte unchanged (verified by an equivalence harness against the old chains). The intentional two-tier split (strict load-time `schemas.py` versus tolerant plot-time validators) is preserved, and the public validator function names are unchanged.
- **Internal: the file-format loaders are table-driven.** The loader functions in `loaders.py` shared one copy-pasted shape (read a delimited file, rename source columns to standard names, log a count, validate), with the p-value-preference branch hand-rolled four times and the first-present column loop repeated eleven times. That shape is now a frozen `LoaderSpec` plus one `_load_tabular` engine: each static format is a spec constant and a thin public wrapper, p-value preference is `p_candidates`, and first-present column selection is `col_candidates`. `load_gtf`, `load_bed`, `load_caviar`, and the GTEx `variant_id` split remain bespoke where the table did not fit. Every public loader keeps its name, signature, return value, and log message; there is no behavior change, verified byte-identical against a pre-refactor output snapshot.
- **Internal: heatmap-highlight geometry and hover tooltip formats have one owner each.** The bounds guard and row-then-column cell walk inside `highlight_heatmap_snp` were maintained separately in all three adapters; they are now `composition.heatmap_highlight_cells`. `hover.py`'s `build_plotly_template` and `build_bokeh_tooltips` had no production caller while each backend re-derived the same column-name-to-format heuristic inline in `scatter`; they are now module-level `plotly_hovertemplate` and `bokeh_tooltips`, called by both backends, verified byte-identical to the old inline logic across 6556 column-name combinations. `HoverDataBuilder` keeps `build_dataframe`, which is live. Also removed: the `_convert_label` identity wrappers in both interactive backends, and five hand-written copies of Bokeh's scalar-to-sequence broadcast.
- **Internal: `ColocPlotter.plot_coloc` is decomposed into named steps.** The method ran about 110 lines and thirteen decision points between building its config and calling the renderer. It is now `_merge_and_transform`, `_assign_colors`, and `_resolve_lead_idx`, module-level functions passing a frozen `_MergedColoc` that carries the merged frame beside its post-merge column names, so the suffix resolution happens once. `_resolve_lead_idx` no longer writes a `combined_score` column onto the frame handed to the renderer; nothing read it. Rendered output is unchanged.
- **Internal: every family renderer now has a `RecordingBackend` contract test.** `ColocRenderer`, `LDHeatmapRenderer`, and `StatsRenderer` had none, unlike `ManhattanQQRenderer` and `MiamiRenderer`. They also called `create_figure` positionally while the covered two used keywords; all five now use keywords, which is what let the recording double cover them. The double gained the missing required primitives, and a `FullCapabilityBackend` subclass supplies the optional heatmap and bar-chart methods, so plain `RecordingBackend` keeps doubling as the backend that declines every capability.
- **Internal: deduplicated shared drawing logic across three modules; no behavior change.** A Miami plot is a mirrored Manhattan, so the single/stacked Manhattan renderer (`_rendering.py`) and the Miami renderer (`_miami_renderer.py`) now share their per-chromosome scatter loop, y-limit padding, and cumulative-position x-limit padding through a new `_manhattan_panel.py`. In `colors.py` the credible-set and PheWAS colour functions route their cyclic palette indexing through one `_cyclic_color` helper (credible sets keep their 1-based `cs_id - 1` offset and `< 1` guard; PheWAS keeps its 0-based index). In `plotter.py`, `plot` and `plot_stacked` build the LD-heatmap and gene-track panels through shared `_build_heatmap_panel` and `_build_gene_panel` helpers. Public API and rendered output are unchanged.

### Removed

- **Public `backend_name` on `MiamiPlotter`, `ColocPlotter`, and `LDHeatmapPlotter`, and `LDHeatmapPlotter.species`.** Nothing read any of them; `LocusZoomPlotter` already kept the backend name private, so it is now private everywhere. `LDHeatmapPlotter` still accepts `species`, documented plainly as accepted for signature symmetry and ignored, since a heatmap renders a precomputed matrix. `MiamiPlotter.species` is unaffected: it feeds `prepare_manhattan_data` on both panels.

### Fixed

- **`ManhattanPlotter`, `MiamiPlotter`, and `StatsPlotter` honour their constructor `genomewide_threshold`.** All three accepted and stored the documented argument, but every plot method took its own threshold defaulting to the module constant, so the constructor value was discarded: `ManhattanPlotter(genomewide_threshold=1e-3).plot_manhattan(df)` drew the line at 5e-8. Only `LocusZoomPlotter` honoured it. Those method arguments now default to an internal `UNSET` sentinel that resolves to the plotter's threshold; passing `None` still means "draw no line", and an explicit p-value still wins. Default construction is unchanged, so plots change only for callers who passed the argument and were being ignored. Related: `StatsRenderer.render_phewas` raised `TypeError` on `-np.log10(None)` rather than skipping the line, so PheWAS never supported the `None` it documented.
- **Plotly LD heatmaps honour `show_colorbar` and the metric label.** `PlotlyBackend.add_colorbar` was a no-op on the reasoning that `add_heatmap` already set `showscale=True`, which made both of the caller's arguments dead. A Plotly heatmap rendered its colour scale even with `show_colorbar=False`, and always titled it `R²` even for a D' matrix. `add_heatmap` now leaves the scale off and `add_colorbar` turns it on with the caller's label and orientation, matching how matplotlib and Bokeh already treat the two calls. `add_heatmap` also returns the figure's own trace instead of the detached `go.Heatmap` it constructed, since `add_trace` stores a copy. Matplotlib and Bokeh are unaffected.
- **Recombination and Ensembl caches resolve through one shared root.** `recombination.get_default_data_dir` and `ensembl.get_ensembl_cache_dir` each hand-rolled platform detection for the same `~/.cache/pylocuszoom` root and disagreed. The Ensembl resolver ignored `$XDG_CACHE_HOME` on macOS and had no Databricks branch, so the two caches pointed to different roots on Databricks and diverged on macOS when `XDG_CACHE_HOME` was set. Both now derive from a single `utils._platform_cache_base()` and append their own leaf (`recombination_maps` / `ensembl`), so the Ensembl cache honours `$XDG_CACHE_HOME` everywhere and routes to `/dbfs/FileStore/reference_data/ensembl` on Databricks. This changes the Ensembl cache location in three cases: on macOS with `XDG_CACHE_HOME` set, on Databricks, and on Windows when `LOCALAPPDATA` is unset, where the fallback moves from `~/AppData/Local/pylocuszoom/ensembl` to `~/pylocuszoom/ensembl`. Windows installs with `LOCALAPPDATA` set (the normal case) are unaffected. A relocated cache is a cold cache, not a broken one.

## [2.0.0] - 2026-07-21

### Changed

- **BREAKING: the `PlotBackend` extension contract is primitives-only** (ADR-0004). Custom backends written against 1.x need migration; see [Custom backends in 2.0](docs/ARCHITECTURE.md#custom-backends-in-20). No compatibility shim is provided.
  - The five semantic legend methods (`add_ld_legend`, `add_effect_legend`, `add_eqtl_legend`, `add_finemapping_legend`, `add_simple_legend`) and the old `add_legend(handles, labels)` are replaced by one neutral `add_legend(ax, entries, loc, title)` taking backend-neutral `LegendEntry` values. Legend content is built once by pure functions in the new `backends/composition.py`.
  - `add_recombination_overlay` is removed. `composition.render_recombination_overlay()` composes it from primitives. `create_twin_axis(ax)` now returns an opaque `SecondaryAxis` handle that every `*_secondary` primitive takes as its first argument, replacing the divergent `ax + yaxis_name` convention.
  - The `supports_snp_labels` and `supports_secondary_axis` boolean properties are removed. Optional capabilities are negotiated with `@runtime_checkable` protocols (`SupportsRegionHighlight`, `SupportsSNPLabels`, `SupportsSecondaryAxis`) checked by `isinstance`. `supports_hover` stays a boolean.
- **Legend styling is uniform across plot families.** The five per-family legend styles are replaced by one. Visible differences from 1.x: the LD lead-SNP marker is 7pt rather than 6pt; eQTL, fine-mapping, and effect legend text is 9pt rather than 8pt, making those boxes slightly taller; and the effect legend now carries an "Effect" title. Legend content, colours, and ordering are unchanged.
- **`schemas.py` validators are expressed on `DataFrameValidator`** instead of a parallel hand-rolled implementation. Which frames are accepted or rejected is unchanged; error wording differs where the engine phrases a check differently (`Missing columns: [...]` for `Missing required column`, `values <= 0` for `non-positive`, `null values` for `missing values`). A frame with a missing column now also reports its dtype and range faults instead of stopping at the first phase.
- **BREAKING: load-time validation rejects nulls consistently.** `validate_gwas_dataframe` already required non-null positions and p-values; `validate_eqtl_dataframe`, `validate_finemapping_dataframe`, and `validate_genes_dataframe` did not, so a file with a blank `p_value`, `pip`, `start`, or `end` loaded silently and lost those rows later at plot time. All four now reject nulls in the numeric columns they range-check, naming the offending column. Plot-time intake is unchanged and still filters nulls with a warning, which is the intended two-tier split.
- **BREAKING: two never-read parameters are removed**: `validate_gwas_dataframe(rs_col=...)` and `validate_finemapping_dataframe(cs_col=...)`. Both were accepted and documented but never used.
- **Internal**: coordinate liftover moved behind a `CoordinateLifter` port with a pure `liftover_positions` (`_liftover.py`), and recombination header detection was extracted as the pure `ensure_recomb_header`. `liftover_recombination_map` keeps its signature.

### Fixed

- **Plotly and Bokeh legends honour `loc` and `edgecolor`.** Both accepted a matplotlib-style `loc` and a per-entry `edgecolor` and discarded them, hard-coding an upper-right legend with black swatch borders. `add_legend(..., loc="lower left")` now places the legend where asked, so it no longer covers plotted data. An unmapped `loc` (including matplotlib's `"best"`, which neither library has) falls back to upper right.
- **The LD legend title renders as mathtext again.** Hoisting the legends replaced `r"$r^2$"` with the literal string `"r²"`, losing matplotlib's italic superscript. Interactive backends convert the mathtext through `convert_latex_to_unicode`, whose table also silently dropped the exponent (`$r^2$` → `"r"`); it now yields `"r²"` and `"R²"`.
- **Non-numeric confidence-interval columns report a validation error.** `require_ci_ordering` compared the columns without checking whether `require_numeric` had already flagged one, so a forest DataFrame with a text `effect` or CI column surfaced pandas' `TypeError` instead of `ForestValidationError`.

## [1.6.0] - 2026-07-20

### Changed

- **Semantic rendering seams across all plotter families** (#22): every plotter (regional, stacked, Manhattan/QQ, Miami, PheWAS/forest, colocalisation, LD heatmap) now composes through explicit seams instead of one method doing data intake, layout, and rendering inline. New internal modules `_data.py` (p-value intake and the shared `P_VALUE_FLOOR`), `_rendering.py`, `_family_renderers.py`, and `_regional.py` own their respective concerns, and backends declare capability flags (e.g. `supports_secondary_axis`) through the base protocol. The compatibility adapter methods left on `LocusZoomPlotter` during the migration (`_plot_association`, `_create_figure_with_heatmap`, `_render_heatmap_panel`, `_highlight_heatmap_snp`, `_add_recombination_overlay`) were removed; callers use `RegionalPlotComposer` directly. The public plotting API and rendered output are unchanged.
- **CI tooling bumps**: GitHub Actions updated (checkout 6.0.2 → 7.0.1, setup-node 6.4.0 → 7.0.0, setup-uv 8.1.0 → 8.3.2, attest-build-provenance 4.1.0 → 4.1.1, lychee-action 2.8.0 → 2.9.0), and the redundant `include_fragments = false` was dropped from `lychee.toml` for lychee 0.24 compatibility.

### Security

- **Bumped vulnerable transitive dependencies** (#23): tornado 6.5.5 → 6.5.7, urllib3 2.6.3 → 2.7.0, and idna 3.11 → 3.18, resolving seven Dependabot advisories (four high, two medium, one low). All three are transitive (via bokeh and requests); this is a lockfile-only change within the parents' version constraints.

## [1.5.1] - 2026-07-20

### Changed

- **LD enrichment has one canonical path**: single and stacked regional plots now share the same LD calculation, merge, warning, and recovery logic in `_ld_plotting.py`. `LocusZoomPlotter` no longer duplicates that orchestration.
- **CODEMAP uses stable file links**: removed source anchor comments, volatile line-number links, the custom synchronization script, and its pre-commit hook. Architecture documentation no longer dictates production source layout.
- **Recombination maps publish as complete generations**: downloads stage and validate all 38 canine autosomal average maps, then atomically switch the active generation while preserving ancillary files such as the CanFam3-to-CanFam4 chain.

### Fixed

- **Empty LD output has a typed contract**: `parse_ld_output()` now raises `EmptyLDOutputError`, a `PlinkError` subclass. Plotting code catches that type directly instead of changing behavior based on exception message text, and stacked plots now recover consistently with single plots.
- **Canine recombination cache could never be complete**: the upstream archive contains average, female, and male maps for chromosomes 1-38 and no X map. The old downloader overwrote each chromosome three times based on archive order, then required 39 cached files. It now selects average maps explicitly and validates the exact 38-file manifest.
- **Partial recombination updates**: a failed download or incomplete archive can no longer leave a mixed old/new map set that later passes the cache check.

## [1.5.0] - 2026-07-20

### Removed

- **Deprecated matplotlib-only `plot_gene_track`** (pyLocusZoom-x61): removed `plot_gene_track` and `_draw_strand_arrows_matplotlib` from `gene_track.py`, along with the `plot_gene_track` export from `pylocuszoom.__init__`. Use `plot_gene_track_generic` (backend-agnostic) instead. The function has carried a `DeprecationWarning` since 1.3.0.

### Changed

- **CLAUDE.md and AGENTS.md are now gitignored**: these files hold agent-local instructions, not project docs. They should never have been committed; the repo now enforces this via `.gitignore` plus a `no-gitignored-files` pre-commit hook. Project-facing setup and release guidance lives in `CONTRIBUTING.md`, `README.md`, and `docs/`.
- **Empty PLINK LD output now raises `PlinkError`**: `parse_ld_output()` previously returned an empty DataFrame when PLINK exited 0 but produced no LD pairs (lead SNP monomorphic, MAF-filtered, or absent from the `.bim` file), leading to silently uncoloured plots with no diagnostic. It now raises `PlinkError` with an actionable message naming the likely cause and the lead SNP. `LocusZoomPlotter.plot()` narrows the catch: the "empty LD output" case is downgraded to a warning and the plot continues without LD colouring (common with singleton lead SNPs), while timeouts, non-zero exits, and missing output propagate to the caller.
- **Gene-track height calculation** (pyLocusZoom-9dm): extracted duplicated `1.0 + (n_rows - 1) * 0.5` logic from `LocusZoomPlotter.plot()` and `plot_stacked()` into a shared `calculate_gene_track_height()` helper in `_plotter_utils.py`.
- **Docs: CODEMAP anchor drift eliminated**: every `docs/CODEMAP.md` table row now points at a per-symbol `# [ID:Symbol]` anchor comment in the source; `scripts/sync-codemap.py` regenerates the line numbers and runs in pre-commit so the doc can never silently drift when code shifts. Module-level anchors previously placed above module docstrings are now below, so Sphinx/`help()`/IDE `__doc__` pickup is no longer suppressed. Directory tree in `docs/ARCHITECTURE.md` no longer lists the gitignored `CLAUDE.md`, and the broken `CONTRIBUTING.md#docstring-example` fragment in `docs/DEVELOPMENT.md` is replaced with a live source reference.

### Fixed

- **LD colouring disabled by duplicate lead SNP**: `parse_ld_output()` appended an explicit lead row (R2=1.0) on top of PLINK's own lead self-pair, so the lead SNP appeared twice in the returned DataFrame. The duplicate key then broke the `validate="many_to_one"` LD merge in `LocusZoomPlotter.plot()` with a `MergeError`, so every plot built from a real PLINK `--ld-snp` run (plink1.9 emits the self-pair) raised instead of colouring by r². The self-pair is now dropped before the canonical lead row is added. Existing `parse_ld_output` fixtures omitted the self-pair, so the regression was invisible to the suite; a new test asserts the returned `SNP` column is unique.
- **NaN p-values silently passing validation** (pyLocusZoom-3oy, pyLocusZoom-y0l): `_validate_coloc_df` and `validate_phewas_df` now reject NaN p-values rather than letting them slip through `require_range`'s `.dropna()` and blank downstream axes. Added `require_not_null([p_col])` between `require_numeric` and `require_range` in the three affected validators.
- **String p-values crashing range check** (pyLocusZoom-y23): `DataFrameValidator` now short-circuits `require_range` when `require_numeric` has already flagged the column, surfacing a structured `ValidationError` instead of a bare `TypeError` from the bound comparison.
- **Recombination `cM` column silent coercion** (pyLocusZoom-z7a): `load_recombination_map()` now logs a warning (with a sample of offending values) when `pd.to_numeric(errors="coerce")` drops non-numeric `cM` values, matching the existing behaviour for `pos`/`rate`.

- **Label backfill regression**: `add_snp_labels()` now filters near-lead non-lead SNPs *before* selecting the top N, so a strong peak no longer collapses to a single label when multiple top hits cluster around the lead.
- **Cross-chromosome lead in stacked plots**: `plot_stacked()` lead auto-detection now filters by chromosome via `filter_by_region`, recognizing both `chrom` and `chr` column conventions, preventing the diamond marker from being placed at a position from the wrong chromosome on multi-chromosome GWAS DataFrames.
- **eQTL p-value validation**: `prepare_eqtl_for_plotting()` now drops rows with NaN, non-positive, or `>1` p-values (with a warning) instead of producing negative or NaN `neglog10p` values that poisoned axis limits.
- **Gene track strand arrows on missing data**: NaN strand values now render with the neutral color and skip directional arrows, instead of silently picking the wrong direction.
- **Global matplotlib state leak**: Removed `plt.ioff()` / `plt.ion()` calls that flipped matplotlib's interactive mode globally, breaking auto-display in caller notebooks.
- **`lead_pos=0` silently disabling LD**: Truthiness checks (`if lead_pos`) replaced with `is not None`, so a valid lead position of 0 no longer skips LD calculation.
- **PLINK prefix with dots**: `validate_plink_files()` no longer truncates prefixes containing `.` (e.g. `ukbb.v3`); existence checks now use string concatenation rather than `Path.with_suffix()`.
- **Degenerate y-axis on flat data**: Manhattan and Miami plots clamp `ylim` to a minimum of 1.0 when every p-value rounds to 1, preventing matplotlib's singular-axis warning. The clamp now also applies to `plot_manhattan_stacked`, `plot_manhattan_qq`, and `plot_manhattan_qq_stacked`, which previously used a bare `y_max * 1.1`. Factored into a shared `_padded_ymax` helper.
- **Recombination map header detection**: `download_canine_recombination_maps()` detects headers by attempting numeric parsing of the first token rather than case-sensitive string matching, so mirrors shipping `Chr`/`CHR` headers no longer corrupt the saved file.
- **PheWAS NaN-effect rows dropped**: Variants with missing effect direction are now rendered as circles instead of being silently filtered out by `>= 0` / `< 0` comparisons.
- **Defensive copy in `_plot_association`**: Adds `df = df.copy()` before assigning `ld_bin`, removing a latent landmine where the function could mutate caller-owned DataFrames.

## [1.4.1] - 2026-04-14

### Security

- **Dependency upgrades** resolving 8 Dependabot alerts (pillow, tornado, orjson, pytest, requests, pygments) via `uv lock --upgrade`
- **SHA-pinned GitHub Actions** in CI and publish workflows to prevent tag-hijack supply-chain attacks
- **Sigstore build provenance attestations** added to PyPI publish workflow — consumers can now verify package provenance
- **Pinned `hatchling==1.29.0`** build backend to prevent compromised build-time dep execution
- Added Dependabot config for weekly github-actions updates

## [1.4.0] - 2026-04-02

### Added

- **Lead SNP proximity filtering for labels**: Non-lead SNPs within 5% of the region width of the lead SNP are excluded from labeling, eliminating the ugly triangle/fan of connector lines when multiple top SNPs cluster near the peak. The lead SNP is always labeled. Controlled by new `lead_pos`, `region_span`, and `min_label_distance` parameters on `add_snp_labels()`.

### Fixed

- **Plotly/Bokeh backend signature mismatch**: Updated `add_snp_labels` signatures on Plotly and Bokeh backends to match the protocol, preventing `TypeError` when using interactive backends
- **Silent no-op on partial parameters**: Now logs a warning when `lead_pos` is provided without a valid `region_span` instead of silently skipping filtering
- **Input validation for `min_label_distance`**: Raises `ValueError` if value is outside `[0, 1]` range

## [1.3.9] - 2026-04-02

### Fixed

- **adjustText arrow artifacts**: Hide connector lines by setting arrow color to `none`, eliminating gray line artifacts from repositioned labels
- **Lead SNP marker clipping**: Increase y-axis headroom from 5% to 15% so diamond markers and labels render without clipping

## [1.3.8] - 2026-04-01

### Fixed

- **Lead SNP marker clipping**: Set axis limits (`set_ylim`, `set_xlim`) before `add_snp_labels` so adjustText has finalized bounds for label positioning
- **adjustText arrow artifacts**: Switched to single-pass `adjust=True` instead of deferred two-pass adjustment, eliminating stale FancyArrowPatch remnants
- **Y-axis headroom**: Added 5% headroom above max data point to prevent top markers from being clipped

### Changed

- Removed `submit-bioconda-pr` job from publish workflow; BiocondaBot handles recipe updates automatically

## [1.3.7] - 2026-03-17

### Fixed

- **eQTL effect size bin gap**: Added missing bins for near-zero effects (0.0-0.1, -0.1-0.0) and fixed negative fallback returning most extreme instead of least extreme bin
- **Feline chromosome support**: `RegionConfig.chrom` now accepts string chromosomes (e.g., "A1", "X") with validation rejecting empty strings and non-positive integers
- **pyliftover ImportError**: Wrapped optional pyliftover import with actionable install message; plotter surfaces it as `UserWarning` visible in notebooks
- **eQTL significance line**: Added `eqtl_threshold` parameter to `plot_stacked()` so eQTL panel uses its own threshold instead of genome-wide significance
- **Stacked plot height**: `plot_stacked()` now uses `max(figsize[1], sum(height_ratios))` instead of ignoring user-provided height
- **LD file path safety**: Sanitize SNP IDs containing special characters (e.g., `chr1:12345:A:G`) for safe use in temp file paths
- **Global RNG pollution**: Categorical Manhattan jitter now uses local `numpy.random.default_rng(42)` instead of mutating global `numpy.random.seed(42)`
- **Plotly monkey-patching**: Replaced `fig._n_cols`/`fig._n_rows` with explicit tuple format; added length guard for unexpected tuple sizes
- **Empty DataFrame validation**: `validate_gwas_df()` now raises `ValidationError` early for empty DataFrames
- **Bokeh heatmap cell sizing**: LD heatmap uses per-cell midpoint-based sizing for correct rendering with non-uniform SNP spacing
- **NaN in range validation**: `DataFrameValidator.require_range()` now drops NaN before range checks to avoid false positives

## [1.3.6] - 2026-02-22

### Fixed

- PyPI project URLs now use correct repository casing (`pyLocusZoom`) for Trusted Publisher verification
- Added Issues and Changelog links to PyPI project metadata

## [1.3.5] - 2026-02-20

### Added

- 15 new Bokeh backend tests covering save(), panel labels, hover collisions, color parsing, heatmap batching, colorbar, and legend range

### Fixed

- **Bokeh `save()` polluted global state** via `output_file()` — now passes `filename` directly to `bokeh.io.save()` and raises `ValueError` for unsupported file formats
- **Bokeh `add_panel_label()` broken with `DataRange1d`** — label placed at data coordinates `(0, 0)` when range not yet resolved; now uses screen-unit pixel offsets
- **Bokeh `scatter()` hover column name collision** — hover data columns named `x`, `y`, `color`, or `size` overwrote internal scatter keys; now namespaced as `hover_<col>`
- **Bokeh `_create_color_palette()` only handled 6-digit hex** — 3-digit hex (`#F00`) and named CSS colors (`red`) caused silent errors; new `_parse_color_to_rgb()` handles all formats via matplotlib's `to_rgb()`
- **Bokeh `highlight_heatmap_snp()` created O(n) renderers** — each highlighted cell added a separate renderer; now batches all cells into at most 2 `rect()` calls with `ColumnDataSource`
- Bokeh `add_colorbar()` identity `orientation_map` removed — orientation passed directly
- Bokeh `_ensure_legend_range()` redundant local `ColumnDataSource` import removed — uses module-level import

## [1.3.4] - 2026-02-20

### Added

- Path traversal protection (`_safe_species_dir`) for Ensembl gene cache — prevents `../../` and absolute path escapes
- `_filter_invalid_pvalues()` DRY helper consolidating duplicate p-value filtering in Manhattan/categorical plots
- 22 new tests: path traversal (6), p-value filtering (8), heatmap coordinates (4), categorical integers (2), version (1), absolute path injection (1)

### Changed

- `_filter_invalid_pvalues()` always returns a copy (no inconsistent mutation) and raises `ValueError` when all rows are invalid (prevents blank plots)
- Bokeh/Plotly heatmaps use actual genomic coordinates instead of integer indices
- Version sourced from `importlib.metadata` instead of hardcoded string
- Default pytest addopts exclude integration tests (`-m 'not integration'`)

### Fixed

- Categorical Manhattan plots with integer category columns silently dropped all points (filtered on wrong column)
- Bokeh heatmap cell dimensions now scale with coordinate spacing instead of hardcoded `width=1`
- Plotly heatmap axes showed index values instead of genomic positions
- Misleading code comments (bokeh "lower triangle", manhattan "-1 for unknown")

## [1.3.3] - 2026-02-20

### Added

- Exception hierarchy: `PyLocusZoomError` base, with `PheWASValidationError`, `ForestValidationError`, `EQTLValidationError`, `LoaderValidationError`, `FinemappingValidationError`, `DataDownloadError`, `BackendError` — all backward-compatible (still catchable as `ValueError`/`RuntimeError`)
- `highlight_heatmap_snp()` added to `PlotBackend` protocol — all three backends implement SNP highlighting natively
- Schema validation edge case tests, loader format detection tests, LD test coverage improvements

### Changed

- `LD_BINS` and `EQTL_BINS` converted from bare tuples to `NamedTuple` (`LDBin`, `EQTLBin`) for self-documenting field access
- Extracted `_find_ld_bin()` helper in `colors.py` (DRY with existing `_find_eqtl_bin`)
- Extracted `_validate_or_warn()` helper in `loaders.py` replacing 7 copy-pasted validation blocks
- Extracted `_add_species_flags()` helper in `ld.py` deduplicating species flag dispatch
- `eqtl.py` migrated to `error_class=` validation pattern (matching finemapping/phewas/forest)
- Removed vestigial `_transform_pvalues` wrapper from `LocusZoomPlotter` — callers use `transform_pvalues()` directly
- Capped `pandas>=1.4.0,<3` to avoid pandas 3.0 breaking changes (StringDtype default, Copy-on-Write)

### Fixed

- CI flaky test failures on Python 3.11: loguru capture tests in `test_recombination.py` now use production logger API instead of raw handlers (xdist-safe)
- Consolidated 3 separate `except` blocks in `ensure_recomb_maps()` into one

### Removed

- Dead `AxesType`/`FigureType` type aliases in `base.py` (never referenced outside definition)
- Dead `TestStdlibWrapperDirect` test class (146 lines testing a locally-defined fake wrapper instead of the actual `_StdlibWrapper`)
- Dead `REQUIRED_EQTL_COLS`/`OPTIONAL_EQTL_COLS` constants from `eqtl.py`

## [1.3.2] - 2026-02-12

### Changed

- **BREAKING**: `get_backend()` now raises `ImportError` when plotly/bokeh is not installed instead of silently falling back to matplotlib
- **BREAKING**: Removed `add_recombination_overlay()` standalone function — use the backend's `add_recombination_overlay()` method instead (called automatically by the plotter)
- Recombination overlay rendering moved into backend protocol — all three backends (matplotlib, plotly, bokeh) now implement `add_recombination_overlay()` directly
- PLINK subprocess failures now log at ERROR level with full stderr (was WARNING with truncated output)

### Fixed

- `_plotter_utils.transform_pvalues()` now filters NaN and out-of-range p-values instead of warn-and-proceed (affected `ManhattanPlotter` and `StatsPlotter`)
- LD merge in `plot()` and `plot_stacked()` now uses `validate='many_to_one'` to prevent silent row duplication when PLINK returns duplicate SNP entries
- Ensembl API `raise_on_error` contract: now correctly raises `ValidationError` after retry exhaustion on repeated 429/503 responses (was silently returning `None`)
- Ensembl API `response.json()` now wrapped in try/except for malformed JSON responses
- `disable_logging()` no longer suppresses `logger.error()` calls — errors always reach stderr
- PLINK subprocess calls now have a 5-minute timeout (was unbounded)
- Narrowed `except Exception` in recombination map download to `(RequestException, OSError, TarError)`

### Removed

- Dead `_create_figure()` method from `LocusZoomPlotter` (48 lines, unreachable code)

### Security

- Upgraded Pillow 12.1.0 → 12.1.1 (fixes CVE for out-of-bounds write when loading PSD images)

## [1.3.1] - 2026-02-04

### Fixed

- SNP labels no longer extend outside plot bounds or have crossing connector lines
  - Root cause: `adjustText` was called before axis limits were set
  - Fix: Deferred `adjust_text()` call until after `finalize_layout()`
  - Added `adjust_snp_labels()` function for manual label adjustment

## [1.3.0] - 2026-02-02

### Added

- `MiamiPlotter` class for mirrored Manhattan plots comparing two GWAS datasets
  - Top panel shows -log10(p) ascending, bottom panel shows -log10(p) descending (inverted)
  - Consistent chromosome colors and shared x-axis across both panels
  - Per-panel significance thresholds (`top_threshold`, `bottom_threshold`)
  - Panel labels to identify datasets (`top_label`, `bottom_label`)
  - SNP annotations independent per panel (`top_snp_annotations`, `bottom_snp_annotations`)
  - Region highlighting across both panels (`highlight_regions`)
  - Interactive hover tooltips in plotly/bokeh backends
  - Full support for all three backends (matplotlib, plotly, bokeh)
- `LDHeatmapPlotter` class for triangular LD heatmap visualization
  - Displays pairwise LD (R² or D') as triangular heatmap
  - White-to-red color gradient for LD values
  - Lead SNP highlighting with visual emphasis
  - Colorbar legend with metric label
  - Full support for all three backends (matplotlib, plotly, bokeh)
- `calculate_pairwise_ld()` function for computing pairwise LD matrices via PLINK
- LD heatmap integration in `LocusZoomPlotter.plot()` and `plot_stacked()`
  - New parameters: `ld_heatmap_df`, `ld_heatmap_snp_ids`, `ld_heatmap_height`
  - Heatmap panel automatically aligns x-axis with regional association plot
- `ColocPlotter` class for GWAS-eQTL colocalization scatter plots
  - Scatter plot comparing GWAS vs eQTL -log10(p) values
  - Points colored by LD (R²) with lead SNP
  - Lead SNP labeled on plot
  - Pearson correlation coefficient and p-value displayed
  - Significance threshold reference lines
  - Optional effect direction coloring (green=congruent, red=incongruent)
  - Optional coloc H4 posterior probability display
  - Full support for all three backends (matplotlib, plotly, bokeh)
- `ColocConfig` Pydantic model for colocalization plot configuration

## [1.2.0] - 2026-02-02

### Added

- Property-based testing with Hypothesis library for improved test coverage
- `tests/strategies.py` module with reusable GWAS data generators
- Hypothesis test profiles (ci/dev/debug) for configurable test intensity
- Property tests for validation, colors, plotter, gene track, Manhattan, and QQ modules
- `ensure_recomb_maps()` function exported from package for pre-downloading recombination data
- `plot_finemapping()` function exported from package for standalone fine-mapping plots

### Changed

- **BREAKING**: Removed deprecated wrapper methods from `LocusZoomPlotter`
  - Removed: `plot_manhattan()`, `plot_qq()`, `plot_manhattan_stacked()`, `plot_manhattan_qq()`, `plot_manhattan_qq_stacked()` — use `ManhattanPlotter` instead
  - Removed: `plot_phewas()`, `plot_forest()` — use `StatsPlotter` instead

### Internal

- Test count increased from 667 to 690 with hypothesis property tests
- Removed unused `OPTIONAL_FINEMAPPING_COLS` constant from finemapping module

## [1.1.2] - 2026-01-30

### Fixed

- README documentation links now work on PyPI (use absolute GitHub URLs)

## [1.1.1] - 2026-01-30

### Fixed

- README images now display correctly on PyPI (use absolute GitHub URLs)

## [1.1.0] - 2026-01-30

### Added

- `plot_manhattan()` method for genome-wide Manhattan plots with chromosome coloring
- `plot_qq()` method for QQ plots with 95% confidence bands and genomic inflation factor (λ)
- `plot_manhattan_stacked()` method for comparing multiple GWAS studies in stacked Manhattan plots
- `plot_manhattan_qq()` method for side-by-side Manhattan and QQ plots in a single figure
- `plot_manhattan_qq_stacked()` method for multi-GWAS comparison with Manhattan+QQ pairs
- `create_figure_grid()` backend method for side-by-side subplot layouts
- `set_suptitle()` backend method for overall figure titles
- `manhattan` module for Manhattan plot data preparation (chromosome ordering, colors, cumulative positions)
- `qq` module for QQ plot data preparation (lambda calculation, confidence bands)
- `set_xticks()` method for all backends (matplotlib, plotly, bokeh)
- Categorical Manhattan plot support (PheWAS-style) via `category_col` parameter
- Species aliases support in Manhattan plots (dog→canine, cat→feline)
- `ManhattanPlotter` class for genome-wide Manhattan and QQ plots
- `StatsPlotter` class for PheWAS and forest plots
- `_plotter_utils.py` module with shared constants and helper functions

### Changed

- Manhattan and QQ plot styling: thinner edge linewidth (0.2) for cleaner appearance
- Manhattan plot colors: switched to colorcet glasbey_bw_minc_20_minl_30 palette
- `LocusZoomPlotter` now delegates Manhattan/QQ to `ManhattanPlotter` and PheWAS/forest to `StatsPlotter`
- Consolidated duplicate styling constants into `_plotter_utils.py` (DRY refactoring)

### Internal

- Test coverage improved from 78% to 83%
- Added comprehensive tests for `ManhattanPlotter`, `StatsPlotter`, and plotter utilities

## [1.0.2] - 2026-01-29

### Added

- `colorcet` as required dependency (for Manhattan plot chromosome colors)

### Changed

- README Quick Start example now shows recombination overlay and auto_genes
- README hero image updated to show recombination rate overlay

## [1.0.1] - 2026-01-29

### Changed

- Gene track font size increased from 7pt to 9pt for better readability
- Removed black connecting line between gene arrows in gene track
- Bokeh legend symbols increased from 10px to 14px (Lead SNP from 12px to 16px)
- Bokeh recombination secondary axis tick marks hidden for cleaner appearance

### Fixed

- Plotly recombination overlay now renders on correct panel (fixed secondary y-axis naming conflict with subplot axes)
- Example plots now show exons in all gene tracks

### Internal

- Backend style mappings moved to module-level constants (bokeh, plotly)
- Lazy imports for Bokeh I/O functions to reduce startup time

## [1.0.0] - 2026-01-28

### Added

- Unified exception hierarchy with `PyLocusZoomError` base class
- Custom exceptions: `ValidationError`, `DownloadError`, `LiftoverError`, `DataError`, `PLINKError`, `ConfigurationError`
- Internal Pydantic validation for plot parameters (validates kwargs at call time)
- Error path tests for download failures and validation edge cases
- CI ordering validation for forest plots (`ci_lower <= effect <= ci_upper`)
- P-value validation warnings for NaN and out-of-range values
- Vectorized eQTL/PheWAS scatter calls for better performance

### Changed

- All validation errors now raise `ValidationError` (also a `ValueError` for backward compatibility)
- Test randomization enabled via pytest-randomly (visible in CI output)
- Config classes (`PlotConfig`, `StackedPlotConfig`) are now internal implementation details, not part of public API
- Capped pytest-xdist workers at 8 to prevent terminal issues

### Fixed

- Recombination overlay now uses correct twin axis for matplotlib (no longer distorts GWAS y-limits)
- Mb formatting now applied to gene track axis for interactive backends (Plotly/Bokeh)
- Gene track row assignment algorithm now correctly prevents overlapping genes in same row
- Handle all-NaN p-values in stacked plot lead SNP detection
- Replaced broad `except Exception` blocks with specific exception types (only 1 justified fallback remains)
- Download error handling now catches specific HTTP/network errors

## [0.8.0] - 2026-01-28

### Added

- `set_yticks()` backend method for consistent y-axis labels across all backends
- Shared `convert_latex_to_unicode()` utility for interactive backends
- Automatic gene annotation fetching from Ensembl REST API (`auto_genes=True`)
- `get_genes_for_region()` function to fetch genes from Ensembl with disk caching
- `fetch_genes_from_ensembl()` and `fetch_exons_from_ensembl()` low-level API functions
- `clear_ensembl_cache()` utility to clear cached Ensembl data
- Support for human, mouse, rat, and any Ensembl species
- Retry logic with exponential backoff for Ensembl API resilience
- 5Mb region size validation (Ensembl API limit)
- `DataFrameValidator` builder class for consistent validation across modules
- `filter_by_region()` shared utility for chromosome/position filtering
- `HoverDataBuilder` for constructing hover tooltips across backends
- Backend capability system with `supports_*` properties for feature detection
- Backend registration system with `get_backend()` and automatic fallback
- Pre-commit hook for pytest with coverage enforcement (70% minimum)

### Changed

- Forest plot example now uses odds ratios with `null_value=1.0` (more representative)
- PheWAS and forest plot y-axis labels now work correctly in Plotly and Bokeh backends
- Gene track styling: arrows now 75% height and 10% wider for better proportions
- Gene track labels increased from 5.5pt to 7pt for improved readability
- Migrated eQTL, finemapping, phewas, and forest validation to `DataFrameValidator`
- Plotter now uses capability-based dispatch instead of backend name checks
- Removed empty `__init__` methods from backend classes
- Removed unused matplotlib imports from plotter (now backend-agnostic)

### Fixed

- `load_gwas()` now forwards `**kwargs` to format-specific loaders
- Forest plot validator now checks that effect and CI columns are numeric
- PheWAS validator now checks that p-values are numeric and within (0, 1] range

### Security

- Tar extraction now includes path traversal protection for recombination map downloads

## [0.7.0] - 2026-01-27

## [0.6.0] - 2026-01-27

### Added

- `plot_phewas()` method for phenome-wide association study plots
- `plot_forest()` method for forest plots (meta-analysis visualization)
- PheWAS category color palette with 12 distinct colors
- Forest plot and PheWAS validation utilities
- Backend methods: `axvline()`, `hbar()`, `errorbar_h()` for new plot types
- Example plots for PheWAS and forest plots
- Progress bars (tqdm) for recombination map and liftover chain downloads
- `requests` and `tqdm` as core dependencies for reliable downloads with progress
- `pytest-randomly` and `pytest-xdist` as dev dependencies for test randomization and parallel execution

### Changed

- Bumped minimum Plotly version to 5.15.0 (required for multiple legends feature)
- eQTL loaders now output `effect_size` column instead of `effect` for plotter compatibility
- Download functions now use `requests` with streaming and progress bars instead of `urllib`

### Fixed

- SAIGE loader now prefers SPA-adjusted p-values (`p.value.NA`) over raw p-values when both present
- BED loader now handles BED12 format and files with more than 6 columns
- eQTL panel in `plot_stacked()` now filters by chromosome in addition to position
- Validation errors for non-numeric p-values or positions now show clear "must be numeric" message instead of runtime errors

## [0.5.0] - 2026-01-27

### Added

- Hover tooltips for fine-mapping scatter plots (Plotly/Bokeh backends)
- Hover tooltips for eQTL scatter plots (Plotly/Bokeh backends)
- Interactive HTML example plots for eQTL and fine-mapping (Plotly/Bokeh)
- Comprehensive marker and hover data tests for interactive backends

### Changed

- Plotly/Bokeh backends now hide grid lines for cleaner LocusZoom appearance
- Plotly/Bokeh backends now show black axis lines (matching matplotlib style)
- Plotly/Bokeh gene track panels now hide y-axis (ticks, labels, line, grid)
- Plotly/Bokeh backends now hide minor ticks and zero lines

## [0.4.0] - 2026-01-26

### Added

- **File format loaders** for common GWAS, eQTL, and fine-mapping formats:
  - GWAS: `load_gwas`, `load_plink_assoc`, `load_regenie`, `load_bolt_lmm`, `load_gemma`, `load_saige`, `load_gwas_catalog`
  - eQTL: `load_gtex_eqtl`, `load_eqtl_catalogue`, `load_matrixeqtl`
  - Fine-mapping: `load_susie`, `load_finemap`, `load_caviar`, `load_polyfun`
  - Gene annotations: `load_gtf`, `load_bed`, `load_ensembl_genes`
- Pydantic validation for file loaders with detailed error messages
- `py.typed` marker for PEP 561 type checking support
- Pre-commit configuration for automated linting
- GitHub issue templates for bug reports and feature requests
- Codecov badge in README

### Changed

- eQTL and fine-mapping legends now route through backend protocol (works with all backends)
- Simplified backend code with reduced duplication
- Backend protocol class diagram added to ARCHITECTURE.md

### Fixed

- Additional robustness improvements for edge cases

## [0.3.0] - 2026-01-26

### Added

- Bioconda recipe for conda installation
- `adjustText` moved to default dependencies (was optional)
- **Interactive plotly backend** - use `backend="plotly"` for hover tooltips and pan/zoom
- **Interactive bokeh backend** - use `backend="bokeh"` for dashboard-ready plots

### Changed

- `plot()` and `plot_stacked()` now use backend protocol for all rendering (scatter, line, axes, layout)
- **Gene track now works with all backends** (plotly, bokeh, matplotlib)
- **Recombination overlay now works with all backends** - secondary y-axis with rate line and fill
- **LD legend now works with all backends** - r² color scale (lead SNP highlighted in plot, not legend)
- SNP labels remain matplotlib-only (interactive backends use hover tooltips instead)
- Default `genomewide_threshold` changed from 5e-7 to 5e-8 (standard GWAS significance)
- Gene track strand colors: forward strand now goldenrod (#DAA520), reverse strand light blue (#6BB3FF)
- Gene track directional arrows: black for forward, dark grey for reverse
- Added panel spacing (hspace=0.1) between stacked/fine-mapping panels for visual separation
- Tightened gene track internal spacing for more compact layout

### Fixed

- Bokeh backend `x_range=None` error when creating figures with shared x-axis
- Bokeh backend `legend_label=None` error in scatter plots
- Bokeh backend LD legend not rendering (empty scatter plots don't create legend glyphs)
- Bokeh backend deprecated `FuncTickFormatter` replaced with `CustomJSTickFormatter`
- Bokeh backend deprecated `circle()` method replaced with `scatter(marker=...)`
- Bokeh backend `FIXED_SIZING_MODE` validation warning in column layouts

## [0.2.0] - 2026-01-26

### Added

- Fine-mapping/SuSiE visualization with credible set coloring
- Example plots in `examples/` directory
- Plot generation script for documentation

### Fixed

- Ruff linting and formatting errors
- Bokeh security vulnerability (bumped to >= 3.8.2)
- `plot()` KeyError when `rs_col` column missing with `ld_reference_file` provided
- `plot_stacked()` now validates eQTL DataFrame columns before use
- `plot_stacked()` now validates list lengths for `lead_positions`, `panel_labels`, and `ld_reference_files`
- `calculate_ld()` docstring now documents `ValidationError` for missing PLINK files

### Changed

- Minimum Python version bumped to 3.10 (required by bokeh 3.8.2)
- Renamed species terminology: "dog" → "canine", "cat" → "feline"
- Clarified interactive backend status in README (coming soon)

## [0.1.0] - 2026-01-26

### Added

- Initial release of pyLocusZoom
- Regional association plots with LD coloring
- Gene and exon track visualization
- Recombination rate overlay (canine only)
- Automatic SNP labeling with adjustText
- Species support: Canine (CanFam3.1/CanFam4), Feline (FelCat9), custom
- CanFam4 coordinate liftover via pyliftover
- Stacked plots for multi-GWAS comparison
- eQTL overlay panel support
- PySpark DataFrame support
- Backend infrastructure for matplotlib, plotly, bokeh (matplotlib only active)
- Logging via loguru
- Comprehensive test suite

### Dependencies

- matplotlib >= 3.5.0
- pandas >= 1.4.0
- numpy >= 1.21.0
- loguru >= 0.7.0
- pyliftover >= 0.4
- plotly >= 5.0.0
- bokeh >= 3.8.2
- kaleido >= 0.2.0

[1.3.5]: https://github.com/michael-denyer/pyLocusZoom/compare/v1.3.4...v1.3.5
[1.3.4]: https://github.com/michael-denyer/pyLocusZoom/compare/v1.3.3...v1.3.4
[1.3.3]: https://github.com/michael-denyer/pyLocusZoom/compare/v1.3.2...v1.3.3
[1.3.2]: https://github.com/michael-denyer/pyLocusZoom/compare/v1.3.1...v1.3.2
[1.3.1]: https://github.com/michael-denyer/pyLocusZoom/compare/v1.3.0...v1.3.1
[1.3.0]: https://github.com/michael-denyer/pyLocusZoom/compare/v1.2.0...v1.3.0
[1.2.0]: https://github.com/michael-denyer/pyLocusZoom/compare/v1.1.2...v1.2.0
[1.1.2]: https://github.com/michael-denyer/pyLocusZoom/compare/v1.1.1...v1.1.2
[1.1.1]: https://github.com/michael-denyer/pyLocusZoom/compare/v1.1.0...v1.1.1
[1.1.0]: https://github.com/michael-denyer/pyLocusZoom/compare/v1.0.2...v1.1.0
[1.0.2]: https://github.com/michael-denyer/pyLocusZoom/compare/v1.0.1...v1.0.2
[1.0.1]: https://github.com/michael-denyer/pyLocusZoom/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/michael-denyer/pyLocusZoom/compare/v0.8.0...v1.0.0
[0.8.0]: https://github.com/michael-denyer/pyLocusZoom/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/michael-denyer/pyLocusZoom/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/michael-denyer/pyLocusZoom/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/michael-denyer/pyLocusZoom/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/michael-denyer/pyLocusZoom/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/michael-denyer/pyLocusZoom/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/michael-denyer/pyLocusZoom/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/michael-denyer/pyLocusZoom/releases/tag/v0.1.0
