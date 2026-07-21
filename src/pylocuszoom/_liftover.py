"""Coordinate-liftover port and the pure liftover transform.

The liftover math is a pure per-position loop; its only external dependency is a
``CoordinateLifter`` that maps a ``(chrom, pos)`` to lifted coordinates.
Production uses pyliftover (``PyLiftOverLifter``); tests use an in-memory mapping
(``InMemoryLifter``), so the coordinate math is exercised without network access
or the external library.
"""

from typing import List, Optional, Protocol, Tuple

import pandas as pd

from .logging import logger


class CoordinateLifter(Protocol):
    """Maps a genomic coordinate to lifted coordinates in another build."""

    def convert(self, chrom: str, pos: int) -> List[Tuple]:
        """Return pyliftover-style ``[(chrom, pos, strand, score), ...]`` or ``[]``."""
        ...


class PyLiftOverLifter:
    """Production adapter around ``pyliftover.LiftOver``."""

    def __init__(self, chain_path: object) -> None:
        from pyliftover import LiftOver

        self._lo = LiftOver(str(chain_path))

    def convert(self, chrom: str, pos: int) -> List[Tuple]:
        return self._lo.convert_coordinate(chrom, pos)


class InMemoryLifter:
    """Test adapter mapping ``(chrom, pos)`` to a lifted position via a dict."""

    def __init__(self, mapping: dict) -> None:
        self._mapping = mapping

    def convert(self, chrom: str, pos: int) -> List[Tuple]:
        new_pos = self._mapping.get((chrom, pos))
        if new_pos is None:
            return []
        return [(chrom, new_pos, "+", 0)]


def liftover_positions(
    recomb_df: pd.DataFrame,
    lifter: CoordinateLifter,
    chrom: Optional[int] = None,
) -> pd.DataFrame:
    """Liftover every position in ``recomb_df`` through ``lifter``.

    Positions that fail to map are dropped; the result is sorted by position.
    Requires either a ``chr`` column or the ``chrom`` argument.

    Raises:
        ValueError: If neither a ``chr`` column nor ``chrom`` is provided.
    """
    if "chr" in recomb_df.columns:
        chroms = recomb_df["chr"].astype(str)
    elif chrom is not None:
        chroms = pd.Series([str(chrom)] * len(recomb_df))
    else:
        raise ValueError("Either 'chr' column or chrom parameter required")

    new_positions = []
    keep_mask = []
    for chr_val, pos in zip(chroms, recomb_df["pos"]):
        chr_str = f"chr{chr_val}" if not str(chr_val).startswith("chr") else chr_val
        result = lifter.convert(chr_str, int(pos))
        if result and len(result) > 0:
            _, new_pos, _, _ = result[0]
            new_positions.append(int(new_pos))
            keep_mask.append(True)
        else:
            new_positions.append(None)
            keep_mask.append(False)

    result_df = recomb_df.copy()
    result_df["pos"] = new_positions
    result_df = result_df[keep_mask].copy()

    unmapped = len(recomb_df) - len(result_df)
    if unmapped > 0:
        logger.debug(f"Dropped {unmapped} positions that failed to liftover")

    return result_df.sort_values("pos").reset_index(drop=True)
