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
    lead_pos: Optional[int],
    ld_col: Optional[str],
    pos_col: str,
    rs_col: str,
    start: int,
    end: int,
    plink_path: Optional[str],
    species: str,
    context: str = "plot",
) -> tuple[pd.DataFrame, Optional[str]]:
    """Calculate and merge LD data using one recovery policy."""
    if not reference_file or lead_pos is None or ld_col is not None:
        return df, ld_col

    if rs_col not in df.columns:
        logger.warning(
            f"Cannot calculate LD for {context}: column '{rs_col}' not found "
            "in GWAS data. Provide rs_col or add SNP IDs to the DataFrame."
        )
        return df, ld_col

    lead_snp_row = df[df[pos_col] == lead_pos]
    if lead_snp_row.empty:
        logger.warning(
            f"Lead SNP at position {lead_pos} not found for {context}. "
            "LD coloring will be skipped."
        )
        return df, ld_col

    lead_snp_id = lead_snp_row[rs_col].iloc[0]
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

    enriched = df.merge(
        ld_df,
        left_on=rs_col,
        right_on="SNP",
        how="left",
        validate="many_to_one",
    )
    return enriched, "R2"
