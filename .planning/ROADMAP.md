# Roadmap: pyLocusZoom

## Milestones

- **v1.0 Refactoring** - Phases 1-5 (shipped 2026-01-28)
- **v1.1 Visualization Expansion** - Phases 6-9 (complete)

## Phases

<details>
<summary>v1.0 Refactoring (Phases 1-5) - SHIPPED 2026-01-28</summary>

See git history for v1.0 milestone details. Delivered:
- DataFrameValidator builder pattern
- Backend Adapter with capability flags
- Pydantic configuration objects
- Unified exception hierarchy
- Test randomization and coverage improvements

</details>

### v1.1 Visualization Expansion (Complete)

**Milestone Goal:** Add Miami plots, LD heatmaps, and colocalization plots to enable comparative GWAS visualization and LD analysis.

- [x] **Phase 6: Miami Plot** - Mirrored Manhattan plots comparing two GWAS datasets
- [x] **Phase 7: LD Heatmap Core** - Triangular LD heatmap with pairwise LD calculation
- [x] **Phase 8: LD Heatmap Integration** - LD heatmap as panel below regional plots
- [x] **Phase 9: Colocalization Plot** - GWAS vs eQTL scatter with LD coloring

## Phase Details

### Phase 6: Miami Plot
**Goal**: Users can compare two GWAS datasets in a mirrored Manhattan layout
**Depends on**: Phase 5 (v1.0 complete)
**Requirements**: MIAMI-01, MIAMI-02, MIAMI-03, MIAMI-04, MIAMI-05, MIAMI-06, MIAMI-07, MIAMI-08, MIAMI-09, MIAMI-10
**Success Criteria** (what must be TRUE):
  1. User can create Miami plot from two GWAS DataFrames with top panel showing -log10(p) and bottom panel showing inverted -log10(p)
  2. User can see consistent chromosome colors across both panels with shared x-axis and chromosome labels between panels
  3. User can add genome-wide significance lines (configurable per-panel thresholds) on both panels
  4. User can add panel labels identifying each dataset and annotate SNPs independently on each panel
  5. User can highlight genomic regions across both panels simultaneously and view interactive tooltips (plotly/bokeh)
**Plans**: 2 plans (complete)

Plans:

- [x] 06-01-PLAN.md — Core MiamiPlotter with TDD (mirrored panels, shared x-axis, consistent colors)
- [x] 06-02-PLAN.md — Features (significance lines, labels, annotations, highlighting, tooltips)

### Phase 7: LD Heatmap Core
**Goal**: Users can create standalone triangular LD heatmaps from pairwise LD data
**Depends on**: Phase 6
**Requirements**: LDHM-01, LDHM-02, LDHM-03, LDHM-04, LDHM-05, LDHM-06, INFRA-01, INFRA-03
**Success Criteria** (what must be TRUE):
  1. User can generate pairwise LD matrix via PLINK from variant list and reference panel
  2. User can create triangular LD heatmap displaying R-squared (or D') values with white-to-red color gradient
  3. User can see colorbar legend showing LD metric scale
  4. User can identify lead SNP row/column via visual highlighting in the heatmap
  5. Heatmap renders correctly in all three backends (matplotlib, plotly, bokeh)
**Plans**: 3 plans

Plans:
- [x] 07-01-PLAN.md — Pairwise LD matrix computation (TDD, extends ld.py)
- [x] 07-02-PLAN.md — Backend protocol extension (add_heatmap, add_colorbar, color constants)
- [x] 07-03-PLAN.md — LDHeatmapPlotter class with tests and export

### Phase 8: LD Heatmap Integration
**Goal**: LD heatmap integrates as optional panel below regional association plots
**Depends on**: Phase 7
**Requirements**: LDHM-03, LDHM-07
**Success Criteria** (what must be TRUE):
  1. User can add LD heatmap panel below regional plot with x-axis alignment (SNP positions match)
  2. User can view combined regional + LD heatmap visualization in all three backends
**Plans**: 1 plan

Plans:

- [x] 08-01-PLAN.md — LD heatmap integration with plot() and plot_stacked() (TDD, coordinate transformation)

### Phase 9: Colocalization Plot
**Goal**: Users can visualize GWAS-eQTL colocalization with LD-colored scatter plots
**Depends on**: Phase 8
**Requirements**: COLOC-01, COLOC-02, COLOC-03, COLOC-04, COLOC-05, COLOC-06, COLOC-07, COLOC-08, INFRA-02
**Success Criteria** (what must be TRUE):
  1. User can create scatter plot comparing GWAS -log10(p) vs eQTL -log10(p) for matched SNPs
  2. User can see points colored by LD (R-squared) with lead SNP and lead SNP labeled
  3. User can see significance threshold reference lines and Pearson correlation coefficient with p-value
  4. User can optionally color points by effect direction agreement and display coloc H4 posterior probability
  5. Colocalization plot renders correctly in all three backends (matplotlib, plotly, bokeh)
**Plans**: 2 plans

Plans:
- [x] 09-01-PLAN.md — Core ColocPlotter with TDD (scatter, LD coloring, correlation, multi-backend)
- [x] 09-02-PLAN.md — Optional features (effect direction coloring, H4 display, ColocConfig)

## Progress

**Execution Order:** Phases execute in numeric order: 6 -> 7 -> 8 -> 9

| Phase | Milestone | Plans Complete | Status | Completed |
|-------|-----------|----------------|--------|-----------|
| 6. Miami Plot | v1.1 | 2/2 | Complete | 2026-02-02 |
| 7. LD Heatmap Core | v1.1 | 3/3 | Complete | 2026-02-02 |
| 8. LD Heatmap Integration | v1.1 | 1/1 | Complete | 2026-02-02 |
| 9. Colocalization Plot | v1.1 | 2/2 | Complete | 2026-02-02 |

### Phase 10: Fix all open bugs from codebase review

**Goal:** Fix 23 bugs from codebase review: P1 performance/correctness issues, error handling hardening, DRY violations, dead code removal, and test coverage gaps
**Depends on:** Phase 9
**Plans:** 2/6 plans executed

Plans:
- [ ] 10-01-PLAN.md — P1 bugs: Manhattan vectorization + silent LD failure + PlinkError
- [ ] 10-02-PLAN.md — Recombination module: HTTPError, tar traversal, vague catches, silent coerce, normalize_build
- [ ] 10-03-PLAN.md — Plotter DRY: deduplicate _transform_pvalues, gene track height, species flags, LDConfig validator
- [ ] 10-04-PLAN.md — Dead code removal + silent defaults: gene_track, schemas, eqtl, format detection, DataFrameLike, platform detection
- [ ] 10-05-PLAN.md — Backend protocol: heatmap highlighting, eQTL bin consolidation, LD_BINS NamedTuple, PlotBackend types
- [ ] 10-06-PLAN.md — Test coverage: schemas.py, loaders.py, logging.py + exception hierarchy completion

---
*Roadmap created: 2026-02-02*
*Last updated: 2026-02-19 — Phase 10 planned (6 plans in 2 waves)*
