# ADR 0002: Share Regional Composition and Plot-Data Intake Policy

- Status: accepted
- Date: 2026-07-20

## Context

Single and stacked regional plots repeated the same association-panel policy,
while p-value filtering and `-log10` preparation were implemented separately
for GWAS and eQTL paths. `PlotConfig`, `StackedPlotConfig`, and `ColocConfig`
also validated arguments but their results were discarded by plotters.

## Decision

- Represent each prepared regional panel with a typed panel specification and
  route complete ordered figure plans through `RegionalPlotComposer`. The
  composer owns figure creation, axis assignment, panel dispatch, labels, and
  final layout. Delete the former plotter rendering adapters.
- Centralize p-value intake in `_data.prepare_pvalue_data`, with an explicit
  zero-value mode to preserve Manhattan compatibility while keeping eQTL's
  strict `(0, 1]` domain.
- Assign validated config objects in all three plotter entry points and use
  their normalized nested values for subsequent orchestration.

## Consequences

Association-panel behavior and optional-panel ordering now have one owner.
Future regional modes can submit prepared panels without duplicating axis
indices or layout branches. Data-family-specific validation remains at each
domain boundary, but the numerical p-value transformation no longer drifts.
Config models are executable boundary objects rather than validation-only
calls.
