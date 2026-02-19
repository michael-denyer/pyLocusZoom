# Project State

**Updated:** 2026-02-19 (10-01 complete)

## Project Reference

See: .planning/PROJECT.md (updated 2026-01-28)

**Core value:** Clean, maintainable codebase where adding new backends or validation rules requires minimal code changes.
**Current focus:** Phase 10 - Fix all open bugs from codebase review

## Current Position

Phase: 10 (Fix All Open Bugs from Codebase Review)
Plan: 5 of 6 complete in current phase
Status: In progress
Last activity: 2026-02-19 - Completed 10-01-PLAN.md

Progress: [████████░░] 83%

## Milestone History

| Version | Name | Status | Shipped |
|---------|------|--------|---------|
| v1.0 | Library Refactoring | Shipped | 2026-01-28 |
| v1.1 | Visualization Expansion | Complete | 2026-02-02 |

## Performance Metrics

**Velocity:**
- Total plans completed: 10
- Average duration: 10m 31s
- Total execution time: 1.75 hours

**By Phase:**

| Phase | Plans | Total | Avg/Plan |
|-------|-------|-------|----------|
| 06 | 2 | 10m 37s | 5m 19s |
| 07 | 3 | 30m 47s | 10m 16s |
| 08 | 1 | 14m 3s | 14m 3s |
| 09 | 2 | 17m 4s | 8m 32s |
| 10 | 2 | 57m 18s | 28m 39s |

**Recent Trend:**
- Last 5 plans: 4m 51s, 14m 3s, 8m 29s, 8m 35s
- Trend: Steady execution times

*Updated after each plan completion*

## Accumulated Context

### Decisions

See: .planning/PROJECT.md (Key Decisions table)

Recent decisions affecting current work:

- [v1.1 Planning]: Miami plot extends ManhattanPlotter (research validated)
- [v1.1 Planning]: LD heatmap requires backend protocol extension (heatmap, colorbar methods)
- [v1.1 Planning]: Colocalization uses scatter with LD coloring (not side-by-side regional)
- [06-01]: Y-axis inversion via set_ylim(max, 0) - works across all backends without protocol extension
- [06-01]: Chromosome union computed from both DataFrames for consistent alignment
- [06-02]: SNP annotations use basic add_text() without adjustText collision avoidance
- [06-02]: Region highlighting uses backend-specific code rather than extending PlotBackend protocol
- [07-01]: Return (DataFrame, list[str]) tuple for pairwise LD matrix and SNP IDs
- [07-01]: Raise ValidationError when requested SNPs missing from reference panel
- [07-02]: Plotly add_colorbar is no-op (configured in add_heatmap via showscale)
- [07-03]: SNP highlighting uses rectangle patches/shapes, not protocol extension
- [08-01]: Heatmap at very bottom of panel stack (GWAS -> finemapping -> eQTL -> genes -> LD)
- [08-01]: Use first GWAS DataFrame for SNP-to-position mapping in plot_stacked()
- [08-01]: Extent-based alignment in matplotlib for genomic coordinate heatmaps
- [09-01]: Inner join on position column for GWAS-eQTL merge
- [09-01]: Auto-select lead SNP as highest combined -log10(p) when ld_col provided
- [09-01]: Drop rows with NaN p-values after merge (cannot compute -log10 of NaN)
- [09-02]: Effect direction coloring uses green (congruent) and red (incongruent)
- [09-02]: H4 posterior displayed in bottom-right corner (data coordinates)
- [10-04]: Kept plot_gene_track body intact with deprecation wrapper rather than removing
- [10-04]: Removed all 4 RowModel classes from schemas.py (none exported or instantiated)
- [10-01]: Used map().fillna(0) instead of apply(axis=1) for vectorized Manhattan cumulative position
- [10-01]: PlinkError inherits both PyLocusZoomError and RuntimeError for dual catch capability
- [10-01]: Kept parse functions returning empty DataFrames for missing files (legitimate no-data case)

### Pending Todos

None.

### Blockers/Concerns

None.

### Roadmap Evolution

- Phase 10 added: Fix all open bugs from codebase review

### Quick Tasks Completed

| # | Description | Date | Commit | Directory |
|---|-------------|------|--------|-----------|
| 1 | Refactor recomb overlay dispatch into backend protocol | 2026-02-12 | 30451f9 | [1-refactor-recomb-overlay-dispatch-into-ba](./quick/1-refactor-recomb-overlay-dispatch-into-ba/) |
| 2 | Fix p-value validation: filter invalid instead of warn-and-proceed | 2026-02-12 | 1bf6552 | [2-fix-out-of-range-p-value-policy-clip-or-](./quick/2-fix-out-of-range-p-value-policy-clip-or-/) |

## Session Continuity

Last session: 2026-02-19
Stopped at: Completed 10-01-PLAN.md
Resume file: None

---
*State file created: 2026-01-27*
*Last updated: 2026-02-19 - Plan 10-01 complete (Manhattan vectorization, PlinkError on LD failures)*
