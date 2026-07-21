# ADR 0004: Complete the Rendering Seam and Unify Capability Negotiation

- Status: accepted
- Date: 2026-07-21
- Target: 2.0 (breaking change to the `PlotBackend` extension contract)

## Context

ADR-0001 moved panel composition, axes, legends, and layout above the rendering
seam into semantic renderers, but preserved the existing primitive contract for
custom backends "during migration". Two forms of composition were left below the
seam, duplicated in every adapter:

- `add_recombination_overlay` is 61 identical lines in each of the matplotlib,
  Plotly, and Bokeh backends. Its five supporting secondary-axis primitives
  (`create_twin_axis`, `line_secondary`, `fill_between_secondary`,
  `set_secondary_ylim`, `set_secondary_ylabel`) have no caller anywhere except
  those three identical bodies.
- The five semantic legends (`add_ld_legend`, `add_effect_legend`,
  `add_eqtl_legend`, `add_finemapping_legend`, `add_simple_legend`) reimplement
  the same bin iteration and swatch construction in each adapter. The generic
  `add_legend(handles, labels)` primitive leaks native handle objects and has
  zero callers.

Capability negotiation drifted to three mechanisms: boolean `supports_*`
properties, `isinstance` against `SupportsRegionHighlight`, and `hasattr` on
`add_recombination_overlay` (which the protocol declares mandatory).
Recombination was gated by two of them at once.

## Decision

Complete ADR-0001: legends and the recombination overlay move above the
primitive contract into a new `backends/composition.py` of pure functions that
drive backend primitives.

- Replace the five semantic legends with a backend-neutral `LegendEntry`
  (`label`, `color`, `marker`, `edgecolor`) and a single neutral
  `add_legend(entries, loc, title)` primitive. Legend content is built once by
  pure builders in `composition.py`; each backend renders the entries natively.
  `add_simple_legend`'s reliance on native labeled-artist collection is dropped
  in favour of an explicit entry.
- Normalize the secondary-axis primitives: `create_twin_axis` returns an opaque
  per-backend `SecondaryAxis` handle, and every `*_secondary(handle, ...)` takes
  it. `render_recombination_overlay` composes them once; the matplotlib-only
  spine-hide tail uses the neutral `hide_spines` primitive.
- Negotiate optional capabilities with `@runtime_checkable` capability protocols
  only: add `SupportsSNPLabels` and `SupportsSecondaryAxis`, keep
  `SupportsRegionHighlight`, and delete the `hasattr` guard. `supports_hover`
  remains a boolean property because it is a rendering-quality flag with no
  method to key on.

The required `PlotBackend` protocol therefore sheds the six composite methods and
the moved optional methods, retaining only true primitives. This breaks the
documented custom-backend extension contract and ships in 2.0. No compatibility
shim is added: a shim would recreate the layer #22 just deleted, and custom
backends must change for the neutral `add_legend` signature regardless.

## Considered options

- Keep the composites as deprecated no-op stubs on each backend for one release:
  rejected because it recreates the shim layer #22 deleted and buys nothing, as
  the `add_legend` signature change forces custom backends to migrate anyway.
- Unify capabilities on boolean properties instead of protocols: rejected
  because a boolean flag keeps optional methods on the required protocol (a
  method that returns `False`), whereas a protocol lets a custom backend omit
  the capability structurally, which is the documented extension model.

## Consequences

The recombination overlay and legend composition have one owner each, tested at
the renderer seam with a recording backend rather than duplicated across three
adapters. The primitive protocol shrinks toward genuine primitives. Custom
backends must migrate to the neutral `add_legend`, the `SecondaryAxis` handle
convention, and the capability protocols; this is a 2.0 breaking change recorded
in the changelog.
