# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- **`colors.py` owns the palette (breaking for two moved names).** Seventeen hex literals lived in seven other modules, five of them duplicating a constant `colors.py` already exported, so "what colours does this library use" was a grep and a colour-blind-safe theme would have been an edit in seven places. Every one moves in as a named export: `STRAND_COLORS` and `STRAND_ARROW_COLORS`, `GENE_LABEL_COLOR`, `RECOMB_COLOR`, `SNP_LABEL_COLOR`, `EQTL_MARKER_COLOR`, `QQ_POINT_COLOR` and `QQ_CI_COLOR`, `UNCATEGORISED_COLOR` and `FOREST_MARKER_COLOR`. `gene_track.STRAND_COLORS`, `recombination.RECOMB_COLOR`, `_plotter_utils.QQ_POINT_COLOR` and `_plotter_utils.QQ_CI_COLOR` no longer exist at those names; import them from `pylocuszoom.colors`. `backends/composition.py` imports `RECOMB_COLOR` at module level instead of lazily reaching up into `recombination` inside a function body. `tests/test_colors.py` fails if a six-digit hex literal appears in any `src/pylocuszoom/*.py` other than `colors.py`. Rendered output is unchanged.
- **One `#BEBEBE`, named for what it means.** The grey was written out three times in `_regional_panels.py` and stood for three unrelated "no data" cases. `colors.NO_DATA_COLOR` is the neutral name and `LD_NA_COLOR` stays as the LD-specific alias for the same value. `palette.get(bin_label, "#BEBEBE")` is now `palette[bin_label]`: `get_ld_bin` only ever returns a label the palette carries, so the fallback was unreachable and separately maintained. `LD_BINS` is a tuple, matching the eQTL bins, and `_find_ld_bin` and the positive eQTL branch share one `_first_at_or_above` lookup, which also removes the `StopIteration` a later bin-boundary edit could have caused.

- **One column vocabulary across the package, with a deprecation path (breaking).** `load_gwas` emitted `chr`/`ps`/`p_wald`/`rs`, `ColumnConfig` defaulted to the same, `ManhattanPlotter` to `chrom`/`pos`/`p`, `eqtl.py` to `pos`/`p_value`, and `plotter.py` sniffed `("chrom", "chr")` to cope, so `plot_manhattan(load_gwas(f))` could not work without renaming three columns. `schemas.Canonical` now names the contract: `chr`, `pos`, `p_value`, `rs`. `chr` was chosen over `chrom` because the gene-annotation, fine-mapping and eQTL families already emit it, and `p_value` over `p` because the eQTL and PheWAS families already use it, so the package has one spelling per concept rather than two. Every loader emits the canonical names, `ColumnConfig`, `GenomeWideConfig`, `ColocConfig` and `utils.filter_by_region` default to them, and the chromosome sniff is gone. `Canonical` is exported from `pylocuszoom`.

  Two aliases survive until **5.0.0**, each emitting a `DeprecationWarning` that names the canonical column:

  - Passing `pos_col=`, `p_col=` or `rs_col=` to a GWAS loader still renames the output column. Rename after loading instead.
  - Handing a plotter a frame that carries `ps` or `p_wald` and no canonical column, without naming your own columns, still plots. A frame carrying both is read as canonical, silently.

  ```python
  # before
  gwas_df = load_gwas("results.assoc.txt")        # chr, ps, p_wald, rs
  fig = plotter.plot_manhattan(gwas_df.rename(columns={"chr": "chrom", "ps": "pos", "p_wald": "p"}))

  # after
  gwas_df = load_gwas("results.assoc.txt")        # chr, pos, p_value, rs
  fig = plotter.plot_manhattan(gwas_df)
  ```

- **One helper strips the `chr` prefix from a column.** `.astype(str).str.replace("chr", "", regex=False)` was open-coded at five sites across `loaders.py`, `gene_track.py` and `utils.py`, beside a `normalize_chrom` that only took a scalar. `utils.normalize_chrom_series` is the column-level companion and every site now calls it. `gene_track.get_nearest_gene` searched a window with its own copy of the chromosome match and the overlap test that `filter_genes_by_region` already held; it now calls that function with the window as the region. Rendered output is unchanged.

- **`plot()` and `plot_stacked()` take the config models as values (breaking).** Each declared 26 and 29 keyword parameters, 24 of them identical, and every one was restated in the `PlotConfig.from_kwargs` call and the docstring while twelve were already fields of `PanelInputs`; `from_kwargs` then ran every nested validator twice. Both methods now take the region (`chrom`, `start`, `end`) plus `columns: ColumnConfig`, `display: DisplayConfig`, `ld: LDConfig` and `panels: PanelInputs`, each option declared once on the model that owns it, and `plot_stacked()` keeps its per-panel lists (`lead_positions`, `panel_labels`, `ld_reference_files`). The four models are exported from `pylocuszoom` and are frozen, so one built in a notebook serves many calls. `PlotConfig.from_kwargs` and `StackedPlotConfig.from_kwargs` no longer exist. No compatibility shim is provided; `scripts/migrate_to_config_models.py` rewrites `.py` and `.ipynb` callers in place. ADR-0008 records the decision.

  ```python
  # before
  fig = plotter.plot(gwas_df, chrom=1, start=1_000_000, end=2_000_000,
                     lead_pos=1_500_000, ld_reference_file="ref", genes_df=genes_df)
  fig = plotter.plot_stacked([gwas_a, gwas_b], chrom=1, start=1_000_000, end=2_000_000,
                             lead_positions=[1_500_000, 1_700_000], ld_col="R2",
                             panel_labels=["A", "B"], genes_df=genes_df, label_top_n=1)

  # after
  from pylocuszoom import DisplayConfig, LDConfig, PanelInputs
  fig = plotter.plot(gwas_df, chrom=1, start=1_000_000, end=2_000_000,
                     ld=LDConfig(lead_pos=1_500_000, ld_reference_file="ref"),
                     panels=PanelInputs(genes_df=genes_df))
  fig = plotter.plot_stacked([gwas_a, gwas_b], chrom=1, start=1_000_000, end=2_000_000,
                             lead_positions=[1_500_000, 1_700_000], panel_labels=["A", "B"],
                             ld=LDConfig(ld_col="R2"), panels=PanelInputs(genes_df=genes_df),
                             display=DisplayConfig(label_top_n=1))
  ```

- **`plot()` and `plot_stacked()` take a per-call `significance_threshold`.** Every other family already let a call override the plotter's `genomewide_threshold` or pass `None` for no line; the regional path had only the constructor value. Omit it to inherit the plotter, as before.
- **`auto_genes` is a `DisplayConfig` field, and `plot()` honours it too.** It was a `plot_stacked()` parameter special-cased in the plotter, the one regional option that was not a config field. `DisplayConfig(auto_genes=None)`, the default, inherits the constructor setting; the constructor keyword is unchanged.
- **The Manhattan family takes a `GenomeWideConfig` and validates at the boundary (breaking).** `ManhattanPlotter`'s five methods and `MiamiPlotter.plot_miami` each restated `chrom_col`, `pos_col`, `p_col` and `custom_chrom_order`, imported no config model and ran no boundary check, so a frame lacking a column reached the layout code before it was rejected and an empty frame reached matplotlib. All six now take `config: GenomeWideConfig` (defaults `chrom`, `pos`, `p`, species order) and hand every frame to `manhattan.prepare_genomewide_frames`, which runs `validate_gwas_df` against the configured names, and `rs_col` on `plot_miami`, before any layout. Two error surfaces change: a missing column raises `ValidationError` (`Missing columns: [...]`, still a `ValueError`) instead of the preparation step's `Column '...' not found`, and an empty frame raises `ValidationError` instead of matplotlib's axis-limit error. Every option after the frame is now keyword-only. `validate_gwas_df` gains `chrom_col`. `GenomeWideConfig` is exported from `pylocuszoom`.

  ```python
  # before
  fig = plotter.plot_manhattan(gwas_df, chrom_col="CHR", pos_col="BP", p_col="P",
                               significance_threshold=5e-8)

  # after
  from pylocuszoom import GenomeWideConfig
  fig = plotter.plot_manhattan(gwas_df, config=GenomeWideConfig(chrom_col="CHR", pos_col="BP", p_col="P"),
                               significance_threshold=5e-8)
  ```

- **`ColocPlotter` colours the merged frame inside the constructor that builds it.** `_MergedColoc` is a frozen value, and `plot_coloc` wrote a `color` column into its payload one line after construction. Internal; rendered output is unchanged.
- **`DisplayConfig.label_top_n` defaults to `None`, meaning the method default.** `plot()` labelled 5 SNPs and `plot_stacked()` 3 per panel, a difference the two signatures carried. With one model serving both, `None` keeps each method's default and an explicit value applies to whichever method receives it. Rendered output is unchanged.

- **`LDHeatmapRequest` is `LDHeatmapPanel` and draws itself (breaking).** `pylocuszoom._ld_heatmap_renderer` is now `pylocuszoom._ld_heatmap_panel`; the request is renamed `LDHeatmapPanel`, loses its `figsize` field (the figure's size belongs to the `FigurePlan` the plotter builds), and gains `draw(backend, ax)` holding what `render_ld_heatmap` drew. `render_ld_heatmap` no longer exists. Rendered output is unchanged.

- **`ColocRequest` is `ColocPanel` and draws itself (breaking).** `render_coloc` opened with `create_figure` and closed with `finalize_layout` around one panel's drawing. `pylocuszoom._coloc_renderer` is now `pylocuszoom._coloc_panel`, the request is renamed `ColocPanel` with the same fields, its `draw(backend, ax)` holds what `render_coloc` drew, and `ColocPlotter` puts it on a `FigurePlan`. Both grey threshold lines go through `add_significance_line`. `render_coloc` no longer exists. Rendered output is unchanged.

- **`StatsRenderer` is removed; PheWAS and forest are panels (breaking).** The class held only the backend, and its two methods shared no code: each grouped or sized its frame, created a figure, drew, and finalized. `pylocuszoom._stats_renderer` is now `pylocuszoom._stats_panels`, holding `PhewasPanel` and `ForestPanel`. Each is built through `from_frame`, which does the grouping and the marker sizing the renderer did at draw time, and draws itself through `draw(backend, ax)`; `StatsPlotter` puts one on a `FigurePlan`. `StatsPlotter._renderer` is gone, and the plotter no longer copies the frame before `prepare_pvalue_data` and `from_frame` copy it again. Rendered output is unchanged.
- **One helper draws every significance line.** `add_significance_line` drew the horizontal line for the regional and Manhattan panels; the PheWAS renderer wrote the vertical one inline with the same red dashed hairline. The helper takes `axis="x"` for a vertical line and a `color`, so the PheWAS line goes through it with its `alpha=0.7`. Every backend's `axvline` and `axhline` already defaulted `zorder` to the 1 the helper passes, so nothing draws differently.

- **The Miami figure is two panels on a `FigurePlan` (breaking).** `render_miami` created the figure, drew both halves through `render_manhattan_panel`, annotated them, highlighted regions across both axes, and finalized. `pylocuszoom._miami_renderer` is now `pylocuszoom._miami_panels`: `MiamiRequest` is unchanged, `miami_plan(request)` builds the plan, and a new `MiamiPanel` draws one mirrored half with its SNP annotations. The region highlights become the plan's `highlights`, resolved against the shared layout's offsets when the plan is built. `render_miami` no longer exists. One ordering detail moves: a half's SNP annotations are now drawn right after that half, rather than after both halves, so a Plotly Miami export that has both panel labels and SNP annotations lists its annotations in a different order; nothing drawn changes, and the example exports are unchanged.

- **`ManhattanQQRenderer` is removed; `ManhattanPlotter` builds figure plans (breaking).** The renderer's six methods each created a figure, built one or more panel specs, drew them, and finalized, and `render_manhattan` and `render_categorical` were the same eleven lines around a different spec. Each `plot_*` method now builds its specs and hands a `FigurePlan` to `render_figure`. `categorical_spec` joins `manhattan_spec` in `_manhattan_panel.py` as the one place the category-axis styling lives, and `stacked_manhattan_specs` moves there too. `ManhattanPanelSpec` and `QQPanelSpec` gain `draw(backend, ax)`, forwarding to `render_manhattan_panel` and `render_qq_panel`. The stacked Manhattan figure titles its top panel rather than the figure, so `FigurePlan` carries `first_panel_title` beside `suptitle`. `pylocuszoom._rendering` and `ManhattanQQRenderer` no longer exist, and `ManhattanPlotter._renderer` is gone. Rendered output is unchanged.

- **One figure plan and one renderer for every family (breaking).** A figure was composed three ways: regional plots through `RegionalFigurePlan` and a composer, Manhattan and QQ through six renderer methods that each hand-wrote the figure scaffolding, and the other families through a request value and a `render_*` function that did the same. `create_figure` was called at ten sites, `create_figure_grid` at two and `finalize_layout` at fifteen, so figure-level policy had no owner. A new `_figure.py` holds a frozen `FigurePlan` (the panels in cell order, `figsize`, `height_ratios`, `n_cols`, `width_ratios`, `sharex`, `xlabel`, `mb_xaxis`, `highlights`, `suptitle`, `top`, `hspace`) and `render_figure(backend, plan)`, the only caller of those four backend methods outside `backends/`. Every panel value has `draw(backend, ax)`, every plotter builds a plan, and `docs/adr/0007-one-figure-plan.md` records the decision. `RegionalFigurePlan` and `pylocuszoom._regional` are removed; `LocusZoomPlotter` builds a `FigurePlan` whose `xlabel` names the chromosome and whose `mb_xaxis` formats every axis in megabases. Rendered output is unchanged.

- **Prepared Manhattan and QQ data are typed values, not frames with a side channel (breaking).** `prepare_manhattan_frames` and `prepare_categorical_data` returned frames carrying their `GenomeLayout` or `CategoryLayout` in `DataFrame.attrs["layout"]`, and `prepare_qq_data` stored `lambda_gc` and `n_variants` the same way, so a typed value reached seven read sites through an untyped dict bolted to a frame, invisible to a type checker and failing as a `KeyError` at draw time. `prepare_manhattan_frames` now returns a list of `PreparedManhattan(frame, layout)`, `prepare_categorical_data` returns one, and `prepare_qq_data` returns `PreparedQQ(frame, lambda_gc, n_variants)`. `manhattan_spec` takes the prepared value, and `MiamiRequest.top` and `bottom` are prepared values. `prepare_manhattan_data` is removed: it was the one-element case of `prepare_manhattan_frames`, so call that with a one-element list and take the first result. None of the four names is exported from the package; code importing them from `pylocuszoom.manhattan` or `pylocuszoom.qq` reads `.frame`, `.layout`, `.lambda_gc` and `.n_variants` instead of `attrs`. Rendered output is unchanged.

- **Each regional panel draws itself; `RegionalPlotComposer` is removed (breaking).** The composer was a class holding the backend and the genome-wide threshold, with a `singledispatchmethod` whose five arms were each one call into a `draw_*` function in the module the panels already lived in, threading a `fig` no arm read. The five panel types now carry a `draw(backend, ax)` method holding what the matching `draw_*` function did, and the figure is rendered by `render_figure` (next entry). `AssociationPanel` gains `region` and `genomewide_threshold` fields, and `GenePanel` and `HeatmapPanel` gain `region`, so no draw needs the figure plan; the threshold is `Optional[float]` with `None` drawing no line, the same contract every other family's threshold has. The eQTL panel's significance line goes through `add_significance_line` like the association line: the helper's `zorder=1` is the value every backend already defaulted to, so the earlier note that routing it would change the stacking was wrong. `RegionalPlotComposer`, `draw_association`, `draw_finemapping`, `draw_eqtl`, `draw_genes` and `draw_heatmap` no longer exist, and `LocusZoomPlotter._regional_composer` is gone; a test that spied on `render_panel` patches `AssociationPanel.draw` instead. Rendered output is unchanged.

- **One function decides whether a region gets a recombination overlay, and one caller reports it (breaking).** Three layers made that call independently, so whether the user heard about it depended on which one fired: a download failure warned from inside `ensure_recomb_maps`, an unsupported species was silent apart from a debug log, and a missing map file warned again from the plotter. Both warnings also carried a `stacklevel` pointing at pyLocusZoom's own internals rather than the caller's line. A new `recomb_for_region` returns a `RecombResult` carrying a `RecombStatus` of `OK`, `NO_MAPS_FOR_SPECIES`, `NO_MAP_FOR_CHROMOSOME`, `DOWNLOAD_FAILED` or `LIFTOVER_UNAVAILABLE` and a one-sentence `detail`; it warns about none of them. `LocusZoomPlotter` emits exactly one `UserWarning` for whatever it gets back, at the same `stacklevel` as the gene-track warning beside it, so the file and line the user sees is their own `plot()` call. `recomb_for_region`, `RecombResult` and `RecombStatus` are exported. `ensure_recomb_maps` keeps its signature and its `Path | None` return, but a failed download now raises `DataDownloadError` or `OSError` where it previously warned and returned `None`; `recomb_for_region` is the entry point for callers that would rather degrade. `LocusZoomPlotter._ensure_recomb_maps` is removed and `_get_recomb_for_region` returns a `RecombResult` rather than a frame.

- **The recombination download is species-generic and the canine specifics are a table row (breaking).** `download_canine_recombination_maps` was 135 lines doing five jobs: directory resolution, a cache-hit check duplicating the one `ensure_recomb_maps` had just done, the HTTP download, tar extraction with an inline traversal guard, a three-arm glob guessing at three archive layouts for an archive with one, a four-stage string peel to name each chromosome, and the atomic publish. All of it sat inside the species-generic `RecombSource.download` port, so a second species could reuse none of it. `RecombSource` now carries `url`, `archive_glob` and `chrom_pattern` and no longer carries `download`; a new module-level `download_recombination_maps(source, output_path)` does the work for any source, with `_extract_archive` and `_stage_maps` split out and unit-tested without a network. The two speculative glob arms are gone. A file the glob matches but `chrom_pattern` cannot parse now raises `DataDownloadError`; the old peel wrote it out silently as `chr_recomb.tsv`. The cache-hit check lives in the two callers that own the policy. `download_canine_recombination_maps` keeps its signature, its `force` flag and its behaviour. Code building a `RecombSource` directly, or patching `CANINE_SOURCE.download` in a test, must patch `download_recombination_maps` instead.
- **One helper resolves the recombination map directory.** `download_canine_recombination_maps`, `load_recombination_map` and `ensure_recomb_maps` each spelled the same `data_dir is None` default differently, one of them rebinding a `str` parameter to a `Path` so the annotation was wrong for the rest of the function. All three call `_resolve_map_dir` now.
- **`filter_by_region` takes `chrom_col=None` for a frame with no chromosome column.** The recombination path passed `chrom_col=""`, which worked only because `"" in df.columns` is false. `None` says what is meant.

- **A species is one record, resolved once at the API boundary (breaking).** `species` was a bare string that five modules interpreted with five private tables: `ld` compared it literally against `"canine"`, `recombination` keyed on `"canine"`, `ensembl` also accepted `"dog"`, `manhattan` carried a second table of the same name mapping the opposite direction, and `plotter` an inline default-build dict. Adding a species was a five-place edit with nothing forcing the fifth. A new `pylocuszoom.species` module holds one frozen `Species` record per species, carrying the Ensembl name, the PLINK chromosome-set flags, the default genome build and the whole-genome chromosome order, and one `resolve_species` that folds case and aliases. `Species` and `resolve_species` are exported. `LocusZoomPlotter`, `ManhattanPlotter` and `MiamiPlotter` resolve at construction and store the record, so `plotter.species` is a `Species` and not a string; read `plotter.species.key` for the old value. They still accept a string. A name the table does not carry becomes an Ensembl-only record with no PLINK flags, no default build and no chromosome order, so `species="sus_scrofa"` keeps working for the gene track as the README promises, and the recombination path now warns that it has no maps for it. `manhattan.SPECIES_ALIASES`, `manhattan.CHROMOSOME_ORDERS`, `ensembl.SPECIES_ALIASES` and `LocusZoomPlotter._default_build` are removed. `get_chromosome_order` no longer takes a `Literal`; a species it knows but has no chromosome order for (mouse, rat) now says so instead of reporting the species unknown. `ld.build_ld_command`, `ld.build_pairwise_ld_command`, `calculate_ld`, `calculate_pairwise_ld`, `ensure_recomb_maps`, `load_recombination_map`, `get_recombination_rate_for_region` and `source_for` all accept a name or a record.

- **Four of the five optional backend protocols are folded into `PlotBackend` (breaking).** `SupportsHeatmap`, `SupportsErrorBars`, `SupportsSecondaryAxis` and `SupportsRegionHighlight` were implemented by all three shipped backends, so the six `isinstance` gates that negotiated them guarded branches no backend could reach, and the unreachable branches had grown three different policies for one question: the regional heatmap panel logged and skipped, the Miami highlight returned silently, and the forest plot and the standalone LD heatmap each raised a hand-written `TypeError`. Their methods are required on `PlotBackend` again. `pylocuszoom.backends.SupportsHeatmap`, `SupportsErrorBars`, `SupportsSecondaryAxis` and `SupportsRegionHighlight` no longer exist, `_ld_heatmap_renderer.require_heatmap_backend` is removed, and `LDHeatmapPlotter` and `StatsRenderer.render_forest` no longer raise `TypeError` for a backend missing those methods: an incomplete backend now fails with `AttributeError` where it draws, like any other missing primitive. `SupportsSNPLabels` is the one capability still optional, because adjustText has no plotly or bokeh equivalent and two of the three shipped backends really do decline it. A third-party backend that implemented every protocol needs no change; one that implemented only some must add the rest. Rendered output is unchanged.

- **`create_figure` derives the panel count, and `finalize_layout` drops the three margins nobody set (breaking).** `create_figure(n_panels, height_ratios, ...)` took two parameters encoding one fact, and the adapters disagreed about which was authoritative: Bokeh derived the count from `height_ratios` and returned one panel for `create_figure(3, [1.0], ...)` where matplotlib raised `ValueError` on the same call. Every one of the ten call sites passed both and they always agreed, so nothing was broken, but a contract with two readings of its own arguments eventually produces a wrong figure on exactly one backend. `n_panels` is gone; `len(height_ratios)` is the count. `finalize_layout` loses `left`, `right` and `bottom`, which no caller ever set, and keeps `top` and `hspace`, which eight of its fifteen call sites do set and which matplotlib and Plotly both honour: dropping them shifts eleven matplotlib PNGs and three Plotly exports, so they stay on the contract, documented as advisory in the same sense as `zorder`. Each backend now holds its own side and bottom margins. `add_region_highlight` drops `fig` and takes `axes` alone; the Plotly adapter reads the figure off the panel handle it is already given. A third-party backend must update all three signatures.

- **Every drawing primitive returns `None` (breaking).** `scatter`, `line`, `fill_between`, `axhline`, `axvline`, `add_text`, `add_rectangle`, `add_polygon`, `errorbar_h` and `add_legend` were annotated `-> Any`, which promises a value without saying what it is. No caller in the library bound one of the 40 returns, and the Plotly adapter had already drifted: five of those methods were annotated `-> Any` and returned `None`. The protocol and all three adapters now say `-> None`, so the annotation matches what the code does and a fourth backend is not asked to manufacture a handle. The matplotlib adapter drops the `(line,) = ax.plot(...)` unpack it did only to have something to hand back. One return still crosses the seam and is now named: `add_heatmap` returns a `Mappable` (an alias of `Any`, documented as opaque) that `add_colorbar` consumes; `add_colorbar` itself returns `None`. Code that bound the result of a drawing primitive gets `None` instead; read the artist back off the figure or the axes.

- **`hbar` is removed and `SupportsBarCharts` is now `SupportsErrorBars` (breaking).** The horizontal-bar primitive had no caller anywhere in the library: the forest plot, the sole reason the protocol exists, draws with `errorbar_h`, `scatter` and `axvline`. It cost 106 lines across the protocol and the three adapters, and a required member of an opt-in protocol has to be implemented correctly by anyone writing a fourth backend for a figure that is never drawn. `hbar` leaves `base.py`, `MatplotlibBackend`, `PlotlyBackend` and `BokehBackend`, and the protocol is renamed to the one method it declares. A third-party backend implementing `hbar` still works, but the name is no longer part of any contract, and code importing `pylocuszoom.backends.SupportsBarCharts` must import `SupportsErrorBars`. `StatsRenderer.render_forest` still raises `TypeError` for a backend without `errorbar_h`; the message now says "does not support error bars" and names `errorbar_h` alone.

- **Miami, colocalization and the LD heatmap take one request value instead of a keyword list (breaking).** Each plotter received a bag of keywords, forwarded them verbatim, and the renderer re-declared them under the same names, so every new option cost three edits in lockstep and no type checker could catch the middle one drifting. `plot_miami` forwarded 20 keywords into a 16-parameter `render`, and `plot_coloc` packed 14 arguments into a `ColocConfig` and then unpacked six of them again. The plotter now builds one frozen `MiamiRequest`, `ColocRequest` or `LDHeatmapRequest` and calls a module-level `render_miami`, `render_coloc` or `render_ld_heatmap`. `MiamiRenderer`, `ColocRenderer` and `LDHeatmapRenderer` are removed; each had one public method and one attribute. `LDHeatmapRenderer.__init__` was where a backend without `SupportsHeatmap` was refused, so that check is now `require_heatmap_backend`, which `LDHeatmapPlotter.__init__` calls; the `TypeError` and its message are unchanged, and so is the point in the call where it fires. `ManhattanQQRenderer` and `StatsRenderer` stay classes because each has more than one entry point. The public signatures of `plot_miami`, `plot_coloc` and `plot_ld_heatmap` are unchanged.

- **The gene track is laid out once.** `GenePanel.from_genes` filtered the genes to the region and ran the greedy row assignment to size the panel, then threw the layout away and kept only its maximum; the drawing filtered and assigned again. One `plot()` with a gene track ran the region filter three times and the row assignment twice, and the height and the drawing could disagree if either side's sort or filter changed. `GenePanel` now carries `genes` (region-filtered, sorted by start), `rows`, `exons` (region-filtered) and `height`, and the drawing iterates them. `GenePanel.data` and `GenePanel.exons_df` are renamed to `genes` and `exons`, and the field order changed, both breaks for code building a panel directly.

- **`plot_finemapping` is removed (breaking).** It was a 104-line public function in the validation layer that took `backend` and `ax` as its first two positional parameters and drove `PlotBackend` primitives directly, which the layering rule reserves for the composition layer. Its one caller in the library also re-derived the credible-set list the function already computed internally. The drawing is now `_regional_panels.draw_finemapping`, the credible sets are resolved once in `FinemappingPanel.from_frame` through `get_credible_sets`, and the three scatter branches are one loop over styled point groups. `pylocuszoom.plot_finemapping` no longer exists; use `LocusZoomPlotter.plot(finemapping_df=...)`. Every other name in `finemapping.py` is unchanged, including `prepare_finemapping_for_plotting`, `get_credible_sets` and `filter_by_credible_set`. Rendered output is unchanged.

- **Every regional panel resolves its mode when it is built and draws through one function.** `RegionalPlotComposer._render_eqtl` and `_render_association` decided the panel's mode and assembled its hover contract at draw time from sixteen `in df.columns` probes, and two of the four optional panels drew inside the composer while the other two drew in a domain module. `EqtlPanel` now carries `effect_col` and `hover`, `FinemappingPanel` carries its resolved `cs_col` and `credible_sets`, and `AssociationPanel` carries `hover` and an `ld_col` that is a column of `data` whenever it is not None. The five panel types, their constructors and one `draw_*` function each live in `_regional_panels.py`, so each `render_panel` arm is a single call and "how is this panel drawn" has one answer per panel. `gene_track.plot_gene_track_generic` and `_draw_strand_arrows_generic` are removed in favour of `_regional_panels.draw_genes`, and `gene_track._compute_arrow_geometry` is now `compute_arrow_geometry`; `gene_track.py` keeps the region filter, the row layout and the arrow geometry and calls no backend. `AssociationPanel` gains a required `hover` field and the panel types now import from `pylocuszoom._regional_panels`, both breaks for code building a panel directly. Rendered output is unchanged.

- **One function draws the significance line for the association panel (breaking).** The regional composer drew its own dashed red line, a near-copy of `_plotter_utils.add_significance_line` differing only in opacity. The helper now takes an `alpha` keyword and the composer calls it. `RegionalPlotComposer` takes `genomewide_threshold` (a p-value) instead of `genomewide_line` (an already-transformed value), which is a break for code constructing the composer directly. The eQTL panel keeps drawing its own line: the shared helper forces `zorder=1` and the eQTL line does not set one, so routing it would change how the line stacks against the eQTL points. `plot()` and `plot_stacked()` still take no threshold argument; adding one is a feature, not this refactor.
- **Four unused colour names are removed (breaking).** `colors.get_credible_set_color_palette`, `colors.get_eqtl_bin`, `colors.LD_HEATMAP_CMAP_NAME` and `colors.LD_HEATMAP_MISSING_COLOR` had no caller in the library and no mention in the documentation. `plotter.DEFAULT_GENOMEWIDE_LINE` goes with them for the same reason. None was exported from the `pylocuszoom` package.
- **The eQTL bin lookup drops its boundary special case.** `_find_eqtl_bin` matched a half-open interval and then re-tested the outermost edge in a second clause, in both the positive and the negative loop. It now takes the first bin the effect reaches, so an effect past the outermost edge lands in that edge's bin without a special case. Verified identical to the old lookup on 24023 values covering every bin edge, both infinities, and the neighbouring floats of each boundary.
- **An LD heatmap that cannot be drawn says so (breaking).** `HeatmapPanel.from_matrix` returned `None` after a log warning when the GWAS frame had no SNP id column or when no heatmap SNP fell inside the region, and the figure was then drawn without the panel the caller asked for. Both are faults in what the caller supplied, on the same footing as an `ld_heatmap_df` with no `ld_heatmap_snp_ids`, which already raised. All three now raise `ValueError` with the messages the first two used to log.
- **The recombination overlay is skipped through one path.** A missing map file was logged, a missing pyliftover was passed to `warnings.warn`, and a download failure was warned about upstream in `ensure_recomb_maps`. Three routes to the same outcome with two mechanisms. The plotter now warns once, in the same "Recombination overlay skipped" shape `ensure_recomb_maps` uses, and returns None. A missing map file now raises a `UserWarning` where it previously only wrote a log line. An `ImportError` that is not about pyliftover still propagates, because a broken environment is not a missing overlay.
- **`add_snp_labels` no longer logs about a missing `region_span`.** The branch was unreachable from the library: the only caller computes the span from a region whose end is validated to exceed its start. Passing `lead_pos` without `region_span` still skips proximity filtering, silently now.
- **A stacked figure knows how many panels it has.** `plot_stacked` checked the length of `lead_positions`, `panel_labels` and `ld_reference_files` in three near-identical blocks because the config had no panel count to check them against. `StackedPlotConfig` now takes `n_panels` and one validator loops over the three lists, raising the same messages. `StackedPlotConfig` is constructed with `n_panels` from now on, which is a break for code building one directly.
- **The stacked LD rule no longer overrides its parent by name.** `StackedPlotConfig.validate_ld_requires_lead_pos` relied on pydantic replacing `PlotConfig`'s validator because the two methods happened to share a name, so renaming either would have run both and rejected valid stacked input. There is now one validator on `PlotConfig` over an overridable `_lead_is_set`, which the stacked config widens to accept a `lead_positions` list. Both error messages are unchanged.
- **A stacked plot picks its lead from the prepared frame.** `_auto_lead_positions` re-implemented the `[0, 1]` p-value domain rule that `_data.prepare_pvalue_data` owns, on the raw frame, before the frame was validated. The lead is now the strongest `neglog10p` in the region of the frame the pipeline already prepared, so the domain rule has one owner. The chromosome guard is unchanged. One logging difference: a region with no rows at all now also logs "No valid p-values in region, cannot determine lead SNP", where before only a region with rows but no usable p-value did.
- **Every column contract lives in one table (breaking).** Eight `ColumnSpec` values were spread over four modules, two of them named `GENES_SPEC` with different rules, so which contract applied to a genes frame depended on which module you imported it from. `schemas.py` is now the single registry: `spec(family, tier)` returns the contract for a `Family` at a `Tier`, strict at `Tier.LOAD` for a frame a loader just parsed and permissive at `Tier.PLOT` for one the caller assembled. `validation.py` is the engine and knows no family. Removed with no replacement, all of them with no caller inside the library: `validation.DataFrameValidator`, whose fluent `require_*` chain is now plain accumulation inside `check`; `schemas.validate_gwas_dataframe`, `validate_eqtl_dataframe`, `validate_finemapping_dataframe` and `validate_genes_dataframe`, each of which was one `check` call that returned its argument; and the constants `validation.gwas_spec`, `validation.GENES_SPEC`, `schemas.EQTL_SPEC`, `schemas.FINEMAPPING_SPEC` and `schemas.GENES_SPEC`, which `spec()` replaces. `validation.validate_gwas_df` and `validation.validate_genes_df` moved to `schemas.py`. No rule changed, and no validation error message changed.
- **The three single-contract modules are gone (breaking).** `phewas.py`, `forest.py` and `coloc.py` held one `ColumnSpec` each and nothing else, and `coloc.py` spent two public wrappers on a difference of one string. `validate_phewas_df` and `validate_forest_df` moved to `schemas.py` with their signatures unchanged and are still exported as `pylocuszoom.validate_phewas_df` and `pylocuszoom.validate_forest_df`. `validate_coloc_gwas_df` and `validate_coloc_eqtl_df` are replaced by one `schemas.validate_coloc_df(df, name, pos_col, p_col, rs_col=None)`. Code that imported `pylocuszoom.phewas`, `pylocuszoom.forest` or `pylocuszoom.coloc` as a module must import from `pylocuszoom` or `pylocuszoom.schemas` instead.
- **The optional-panel arguments are validated by the config layer.** `plot()` and `plot_stacked()` declared the same thirteen panel keywords, passed them through `_render_regional` untouched, and enforced the one cross-field rule among them as a hand-written raise deep in the pipeline, while the equivalent LD rule was a `PlotConfig` validator. `PlotConfig` now carries a `panels` field, a `PanelInputs` model that declares each of them once, and `PlotConfig.from_kwargs` routes them like every other keyword. The public signatures of `plot()` and `plot_stacked()` are unchanged. `ld_heatmap_df` without `ld_heatmap_snp_ids` still raises the same message, now as a `ValidationError` from `PanelInputs` rather than a bare `ValueError`, and it is raised when the config is built rather than partway through rendering. A panel argument that is neither a DataFrame nor `None` is now rejected there too, instead of failing later inside the panel it fed.
- **The liftover chain has its own cache directory.** It was downloaded into the recombination-map directory, which made that directory something other than a replaceable set of maps and cost 54 lines of symlinked-generation publishing to work around. Chains now go to a `liftover` leaf beside `recombination_maps`, and a downloaded map set is swapped in with one rename. An existing chain is downloaded once more into the new location; the old copy can be deleted. A map directory published by an earlier release still reads, and is replaced correctly along with the generation it links to when it is next downloaded.
- **Recombination sources are a table.** `ensure_recomb_maps`, `get_recombination_rate_for_region` and `load_recombination_map` each compared the species string to `"canine"` in their own way. A `RecombSource` row now carries the map set, the build it is published in, the liftover chains it offers and its downloader, so a second species is one row rather than three edits.

- **`get_genes_for_build` takes a source and returns one type (breaking).** It took a `species` that a caller passing `source` had to supply anyway, and returned either a frame or a tuple of two depending on `include_exons`, so every caller and every type checker branched on a flag. It is now `get_genes_for_build(source, chrom, start, end, *, cache_dir=None, use_cache=True)` and always returns a `GeneAnnotations` named tuple of `genes` and `exons`. Use `source_for(species, genome_build)` to get the source; that function is now the only place a species or a build is interpreted.
- **The gene cache stores the exons too.** It saved only the gene frame, so a cache hit that also wanted exons went back to the network, which for UCSC meant re-requesting the whole `ncbiRefSeq` track the genes had just come from. Both frames are now written and read as one entry, and half an entry counts as a miss, which retires the gene-only entries earlier releases wrote.
- **The per-source gene wrappers are removed (breaking).** `get_genes_for_region`, `fetch_genes_from_ensembl`, `fetch_exons_from_ensembl`, `clear_ensembl_cache`, `get_genes_for_region_ucsc`, `fetch_genes_from_ucsc`, `fetch_exons_from_ucsc` and `clear_ucsc_cache` were eight names for two operations that `get_genes_for_build` and `clear_gene_cache` already performed, and each source also carried a private copy of the same failure-to-empty helper. `get_genes_for_build` and `clear_gene_cache` are exported in their place, along with `source_for`, which returns the `GeneSource` for a species and build so a caller can still address one source directly. `ensembl.get_ensembl_cache_dir` and `ucsc.get_ucsc_cache_dir` are gone too; `_gene_cache.cache_root` was always what they returned.
- **A gene-source failure always raises (breaking).** `raise_on_error` is removed from every signature. It defaulted to False, which turned an outage into an empty frame indistinguishable from a region with no genes, and the one caller in the library opted out of it. `LocusZoomPlotter` still catches `ReferenceAPIError` under `auto_genes=True`, warns, and draws the plot without a gene track.
- **`GeneSource` and the frame schema move to `_gene_source.py`.** `ensembl.py` and `ucsc.py` imported the shared vocabulary from `reference_genes.py` while `reference_genes.py` imported the sources back inside its function bodies to break the cycle. The vocabulary now lives in a leaf that imports neither source, each source exports one `GeneSource` constructor, and the routing module imports both at the top of the file. `GeneSource.error_cls` is gone: both source errors subclass `ReferenceAPIError`, which the orchestration catches once.

### Fixed

- **A species with no recombination maps no longer skips the overlay in silence.** `ensure_recomb_maps` returned `None` for an unsupported species with only a `logger.debug` line, so `LocusZoomPlotter(species="human")` drew a plot with no overlay and said nothing, while a download failure for the same plot warned. Every reason the overlay is skipped now reaches the user as one `UserWarning` naming the cause.
- **`species="dog"` no longer degrades silently in three subsystems.** The gene layer advertised `dog`, `cat`, `canis_familiaris` and case folding as aliases and resolved them correctly. The same string then reached PLINK, which compared it against `"canine"` and so omitted `--dog`, reading a 38-autosome canine fileset under PLINK's default human chromosome set; reached the recombination source table, which keyed on `"canine"` and returned no maps with only a debug log; and reached the default-build lookup, which returned `None` and thereby suppressed the assembly-mismatch warning, so CanFam3.1 data was served ROS_Cfam_1.0 gene coordinates without comment. One alias, four wrong results, no error. All of them now read the same `Species` record. `"Canine"`, `"cat"` and `"canis_familiaris"` were broken the same way and are fixed the same way.
- **`plot()` auto-detects the lead SNP, as `plot_stacked` already did.** `plot()` always handed the composer a one-element list holding `lead_pos`, so with no lead given it drew no lead diamond, while `plot_stacked([df])` on the same frame and region auto-detected the strongest in-region p-value and drew one. Both entry points now auto-detect when no lead is supplied. A regional plot that previously showed no diamond gains one at its strongest hit.
- **`calculate_pairwise_ld(metric=...)` rejects anything but `r2` and `dprime`.** The command builder chose `--r2` for every spelling other than the exact string `dprime`, so `metric="Dprime"` returned an r² matrix under a D' label. The metric is now a `Literal["r2", "dprime"]` looked up in one flag table, and any other value raises `ValidationError` before PLINK is resolved.
- **Bokeh keeps the legend's y-range when a recombination axis is added.** `create_twin_axis` replaced the figure's `extra_y_ranges` wholesale, discarding the `legend_range` that `add_legend` had registered; the regional plot only survived because the overlay happened to be drawn before the LD legend. The twin axis is now added alongside the existing ranges, in either order.
- **`validate_phewas_df` no longer takes a `category_col` it ignored.** The parameter was declared and documented but never reached the column spec, so a mistyped category name passed validation. The category column is optional at render time, so the validator does not accept the argument at all; passing it is a `TypeError`.
- **A missing `pyliftover` is reported by type, not by message.** The plotter decided between skipping the recombination overlay and crashing by searching the `ImportError` text for the word `pyliftover`. `liftover_recombination_map` now raises `OptionalDependencyMissing`, a new exception that is both a `PyLocusZoomError` and an `ImportError`; the plotter catches that type and lets any other `ImportError` propagate.
- **Plotly megabase ticks follow the axis limits, not the panel's traces.** `PlotlyBackend.format_xaxis_mb` derived its tick range from the panel's own traces whenever the panel had no `range` of its own, and a shared-x figure only sets a range on the panels that call `set_xlim`. The eQTL panel does not, so `eqtl_plotly.html` carried one tick ladder on the association and gene panels (1.00 to 2.00 Mb every 0.25) and another on the eQTL panel (1.50 and 1.60, the span of its points). The range is now read from the axis, or from any axis in its `matches` group, and a panel with no limits anywhere is left to plotly's automatic ticks. The eQTL and finemapping plotly exports change: every stacked panel now shares the region's ladder.

- **A missing recombination map named a function that does not exist.** `load_recombination_map` built its remedy as `download_{species}_recombination_maps()`, so a feline caller was told to run `download_feline_recombination_maps()`. The message now names `ensure_recomb_maps(species=...)` for a species with built-in maps, and says there are none and to supply `data_dir` for a species without.

- **An automatic gene track now draws exon structure.** `LocusZoomPlotter` with `auto_genes=True` asked its gene source for genes only, so every automatically fetched gene was drawn as a plain rectangle. On the Ensembl path the exons were unusable anyway: exon features name their transcript and nothing else, so every fetched exon carried an empty `gene_name` and the renderer, which matches exons to genes by that column, never found one. Ensembl is now asked for genes, transcripts and exons in one request and the transcript joins each exon to its gene symbol, and the plotter asks for exons. An exon whose gene is absent from the response is dropped rather than returned unattached.
- **A PheWAS point with a null category no longer repaints the whole plot.** `render_phewas` looped over `df[category_col].unique()`, and on an object column that keeps `None` rather than `NaN`, so the `cat is None` arm meant for "this frame has no category column at all" matched a real null category. It drew every row of the frame in the no-category blue, on top of the categorised points. A null category is now its own group, `Uncategorised`, and takes the next palette colour.
- **Miami and stacked Manhattan panels line up on x.** Each panel's chromosome offsets came from that panel's own maximum position per chromosome, so a study whose last chr1 variant sat 7 Mb short of the other study's shifted every later chromosome by 7 Mb. The panels still shared one x range and took their ticks from the first frame, so the misaligned panel was labelled with the other panel's chromosome boundaries and nothing in the figure showed it. `prepare_manhattan_frames` now computes one `GenomeLayout` from every frame in the figure and prepares all of them against it. Stacked Manhattan and stacked Manhattan+QQ figures over datasets with different chromosome extents move, and so do their x tick positions.
- **`MiamiPlotter(species=...)` orders the chromosomes again.** `plot_miami` passed the union of both frames' chromosomes as `custom_chrom_order`, which takes precedence over `species`, so the species table was never read and the axis was ordered numerically then alphabetically. On a feline plot that put MT before X and Y. An unknown species also produced a figure instead of the `ValueError` `ManhattanPlotter` raises for it. The union is gone; the species order, restricted to the chromosomes the frames carry, drives the axis.

- **The lead-SNP outline on a regional LD heatmap is drawn in genomic coordinates.** `highlight_heatmap_snp` received only a SNP index and a SNP count, so every adapter outlined unit cells at matrix indices while the regional heatmap panel is drawn at base-pair positions. The outline landed between x = -0.5 and x = 2.5 on an axis spanning the whole region, and on Plotly the shapes carried no axis reference, which bound them to the first subplot. Both renderers now compute the rectangles with `composition.heatmap_highlight_rects(snp_idx, x_coords, y_coords)` and draw them through `add_rectangle`, which targets the panel it is given. Standalone LD heatmaps are unaffected: their coordinates are matrix indices, so the geometry is what it was.
- **Plotly megabase ticks read the panel they are given.** With no explicit range on a panel, `format_xaxis_mb` took the minimum and maximum x of every trace on the figure, whichever subplot each belonged to. A panel showing 50 to 51 Mb beside one showing 1 to 2 Mb was labelled from 5 Mb to 50 Mb. The fallback now considers only traces bound to the panel's own x-axis. The eQTL and fine-mapping panels reach this path, so their tick labels change wherever panels span different regions.
- **Bokeh Manhattan+QQ figures keep their title.** `set_suptitle` looked only at the layout's first child, which on a `create_figure_grid` layout is a `Row` with no title, so `plot_manhattan_qq(title=...)` and `plot_manhattan_qq_stacked(title=...)` silently dropped the title on bokeh. The backend now walks to the first plot in the layout tree.
- **`load_gwas` auto-detects GEMMA `.assoc.txt` files.** The `.assoc` PLINK hint matched first, so every GEMMA default filename went through the PLINK loader. That passed by coincidence when the output column names were left at their GEMMA defaults and raised `LoaderValidationError` for `load_gwas("output.assoc.txt", pos_col="pos")`. Format hints are now tried longest first.
- **Importing pyLocusZoom no longer removes the host application's loguru handlers.** The logging wrapper called `loguru.logger.remove()` with no argument at import, which deletes every handler, not just loguru's default. It now removes only the default stderr handler, so sinks an application configured before `import pylocuszoom` keep receiving messages.

### Added

- **`ColocPlotter(genomewide_threshold=..., eqtl_threshold=...)`.** The colocalization plotter was the only threshold-bearing plotter with no constructor default to inherit.

### Changed

- **`ColocPlotter` thresholds join the `UNSET` model, so `None` means no line.** `plot_coloc(gwas_threshold=..., eqtl_threshold=...)` were plain floats defaulting to `5e-8` and `1e-5`, matched by `gt=0` fields on `ColocConfig`, so a caller could not ask for a plot with no significance lines and there was no plotter-level default to inherit. Both arguments now default to `UNSET` and resolve against the constructor, `None` draws no line, and the `5e-8` and `1e-5` literals are `DEFAULT_GENOMEWIDE_THRESHOLD` and the new `DEFAULT_EQTL_THRESHOLD`. Passing an explicit float still behaves as before.
- **Internal: `render_phewas` groups once instead of branching three ways.** The category loop had a `cat is None` arm, a `pd.isna(cat)` arm and a real-category arm, and inside each of them a duplicated scatter call that was the one-subset case of its own marker loop. It now assigns a `_group` column once, groups by it, and always loops over `_effect_subsets`, which returns one `(mask, marker)` pair when there is no effect column.
- **Internal: `manhattan_spec` names the policy fields it forwards.** It took `**policy: Any` and handed it to the dataclass, so `panel_lable=` was a runtime `TypeError` invisible to a type checker. The nine fields the three call sites vary are now keyword-only parameters, and a contract test fails if their defaults drift from `ManhattanPanelSpec`.
- **Internal: one way to draw a QQ panel.** `_render_qq_panel` drew the points and `_set_qq_labels_and_title` took nine parameters, including a `stacked` flag that only chose between two title strings and a `title_fontsize or fontsize + 2` fallback. The three call sites each passed a different combination. A frozen `QQPanelSpec` in the new `_qq_panel.py` now names what they vary on, `render_qq_panel` draws it, and the pure `qq_title` builds the title up front. Figures are unchanged.
- **Internal: one name for the transformed p-value.** The Manhattan family wrote `_neg_log_p` while `_data.prepare_pvalue_data`, the regional pipeline, the label placer, the stats renderer and the matplotlib backend all used `neglog10p`, so a prepared Manhattan frame could not feed any other family's panel. The Manhattan column is now `neglog10p` and `manhattan._prepare_pvalues`, which existed only to rename it, is gone. Colocalization still needs two transformed columns in one frame and keeps `neglog10_gwas` and `neglog10_eqtl`, but it builds them with `prepare_pvalue_data` before the merge instead of its own clip-and-log, so it emits the same NaN and out-of-range warnings as every other family. Figures are unchanged.
- **`SupportsHeatmap` drops `highlight_heatmap_snp` (breaking for custom backends).** The protocol is now `add_heatmap` and `add_colorbar`. Split cell geometry caused the defect above: `add_heatmap` knew the coordinates and `highlight_heatmap_snp` did not, so three adapters each invented a third rule. A custom backend can delete the method; `runtime_checkable` protocols check names, so one that keeps it still satisfies `SupportsHeatmap`.
- **The Plotly backend hands out `_Panel` values, not tuples.** `create_figure` returned `(fig, row)` pairs and `create_figure_grid` returned `(fig, row, col, n_cols)` quadruples, and every drawing method reparsed the handle through `_Panel.of` before using it. The constructors now return `_Panel` instances directly, `_Panel.of` and its length check are gone, and each signature says `ax: _Panel` where 25 of them claimed `Tuple[go.Figure, int]`, which was false for every grid handle. `_Panel` is still a `NamedTuple`, so code indexing `ax[0]` for the figure keeps working. `create_twin_axis` returns a `_SecondaryAxis` on Plotly and on Bokeh, which lets Bokeh label and scale the axis it created rather than searching the figure for it.
- **`SupportsSecondaryAxis` drops `line_secondary` and `fill_between_secondary` (breaking for custom backends).** Both were near-copies of `line` and `fill_between` in every adapter, and matplotlib's were already one-line delegations. The handle `create_twin_axis` returns is now a drawing target, so `render_recombination_overlay` calls `backend.fill_between(secondary, ...)` and `backend.line(secondary, ...)`, and the protocol is `create_twin_axis`, `set_secondary_ylim`, `set_secondary_ylabel`. A custom backend moves whatever its two methods did into `line` and `fill_between`, which must accept the twin handle. Plotly figures are unchanged in content; primary line and filled traces now name their axes in the trace rather than through `add_trace(row=, col=)`, which reorders keys in exported JSON.
- **`scatter` and `line` drop their `label` argument (breaking for custom backends).** No renderer passed it. It was also a back door around ADR-0004, which moved every legend decision above the seam: bokeh turned a `label` into `legend_label`, so a caller using it would have got a second, auto-built legend beside the one `add_legend` renders. `PlotBackend` now documents `zorder` as advisory, honoured by matplotlib and accepted in call order by the interactive backends.
- **`add_heatmap` requires `cmap_colors` (breaking for custom backends).** The white-to-red default was written out in all three adapters and in the protocol docstring, and both callers pass `LD_HEATMAP_COLORS` explicitly, so the default and its three `is None` branches served nobody.
- **Internal: one axis-style block per backend.** Plotly wrote its panel style twice, and the two copies had drifted: the single-column constructor set `minor_ticks=""` and the grid constructor did not. Both now use one `_AXIS_STYLE` constant without it, because plotly's default for minor ticks is already off, so exported figures lose an inert `minor` key and render the same. Bokeh's six-line panel styling was verbatim in both constructors and is now `_style_panel`. The `left` coercion in `hbar` was byte-identical in both interactive backends and is now `_coerce.per_point`.
- **Internal: `_coerce` names the types it branches on.** `marker_colors` and `broadcast` detected a Series with `hasattr(value, "tolist")` and `hasattr(value, "values")` while their signatures already named the union, and `broadcast` returned `Any`.
- **`add_rectangle` accepts `facecolor=None` for an outline with no fill.** Matplotlib draws it with `fill=False`, Plotly with a transparent `fillcolor`, Bokeh with `fill_color=None`.

## [3.0.0] - 2026-09-02

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

- **`BUILTIN_BACKENDS` is exported from `pylocuszoom.backends`.** The tuple of shipped backend names, derived from the `BackendType` literal so the two cannot drift. Custom-backend authors running the contract tests can parametrize over it instead of writing `["matplotlib", "plotly", "bokeh"]` out again.
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
