"""Significance-threshold defaults and resolution shared by the plotters.

Internal module - not part of public API. Drawing constants and the
significance line live with the panels, in ``panels._shared``.
"""

from typing import Optional, Union

# Significance thresholds
DEFAULT_GENOMEWIDE_THRESHOLD = 5e-8
DEFAULT_EQTL_THRESHOLD = 1e-5


class _Unset:
    """Type of the ``UNSET`` sentinel."""

    def __repr__(self) -> str:
        return "UNSET"


# Marks a threshold argument the caller did not supply, so it falls back to the
# plotter's genomewide_threshold. None cannot carry that meaning: it already
# means "draw no significance line", and a caller must keep being able to say so.
UNSET = _Unset()

# A significance threshold as a caller may express it: a p-value, None for no
# line, or UNSET to inherit the plotter's own threshold.
ThresholdArg = Union[float, None, _Unset]


def resolve_threshold(
    supplied: ThresholdArg, plotter_default: Optional[float]
) -> Optional[float]:
    """Pick the threshold to draw at, honouring an explicit None.

    Args:
        supplied: The caller's threshold argument.
        plotter_default: The plotter's own ``genomewide_threshold``.

    Returns:
        The p-value to draw the line at, or None to draw no line.
    """
    return plotter_default if isinstance(supplied, _Unset) else supplied


# Space left between chromosomes on a genome-wide x axis, in base pairs.
CHROMOSOME_GAP = 1_000_000
