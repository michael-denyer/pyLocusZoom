"""LD enrichment for regional association plots."""

from typing import Optional

import pandas as pd

from .exceptions import EmptyLDOutputError
from .ld import calculate_ld
from .logging import logger


def enrich_with_ld(
    df: pd.DataFrame,
    *,
    reference_file: Optional[str],
    lead_index: Optional[int],
    ld_col: Optional[str],
    rs_col: Optional[str],
    start: int,
    end: int,
    plink_path: Optional[str],
    species: str,
    context: str = "plot",
) -> tuple[pd.DataFrame, Optional[str]]:
    """Assign LD values by SNP ID without changing the selected rows or their index.

    ``lead_index`` identifies a row already selected at the regional boundary.
    The helper never infers a variant ID from a potentially ambiguous position.
    The regional boundary has already required ``rs_col`` whenever
    ``reference_file`` is set.
    """
    if not reference_file or lead_index is None or ld_col is not None:
        return df, ld_col

    lead_snp_id = df.at[lead_index, rs_col]
    logger.debug(f"Calculating LD for lead SNP {lead_snp_id}")
    try:
        ld_df = calculate_ld(
            bfile_path=reference_file,
            lead_snp=lead_snp_id,
            window_kb=max((end - start) // 1000, 500),
            plink_path=plink_path,
            species=species,
        )
    except EmptyLDOutputError as exc:
        logger.warning(
            f"LD calculation skipped for {context}: {exc}. "
            "Proceeding without LD coloring."
        )
        return df, ld_col

    lookup = ld_df.set_index("SNP", verify_integrity=True)["R2"]
    return df.assign(R2=df[rs_col].map(lookup)), "R2"
