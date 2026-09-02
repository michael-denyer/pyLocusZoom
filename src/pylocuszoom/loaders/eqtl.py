"""eQTL result loaders: GTEx, the eQTL Catalogue, and MatrixEQTL.

Every eQTL loader takes ``filepath`` plus an optional ``gene`` to filter to,
and returns a DataFrame with columns pos, p_value, gene, effect_size.

Raises:
    LoaderValidationError: If the format's columns cannot be mapped or the
        mapped values fail the format's schema.
"""

from pathlib import Path
from typing import Optional, Union

import pandas as pd

from ..schemas import Family, Tier, spec
from ._engine import LoaderSpec, _load_tabular


def _gtex_pos(df: pd.DataFrame, out_cols: dict[str, str]) -> pd.DataFrame:
    """Derive GTEx position from the variant_id, else fall back to tss/POS."""
    if "variant_id" in df.columns:
        df["pos"] = df["variant_id"].str.split("_").str[1].astype(int)
    elif "pos" not in df.columns:
        for col in ("tss_distance", "POS"):
            if col in df.columns:
                df = df.rename(columns={col: "pos"})
                break
    return df


_GTEX_SPEC = LoaderSpec(
    log_fmt="Loaded GTEx eQTL file with {n} associations",
    read={"sep": "\t"},
    col_candidates={
        "p_value": ("pval_nominal", "p_value", "pvalue", "P"),
        "gene": ("gene_id", "gene_name", "phenotype_id"),
        "effect_size": ("slope", "beta", "effect_size"),
    },
    transform=_gtex_pos,
    gene_filter="contains",
    schema=lambda out_cols: spec(Family.EQTL, Tier.LOAD),
)


def load_gtex_eqtl(
    filepath: Union[str, Path],
    gene: Optional[str] = None,
) -> pd.DataFrame:
    """Load GTEx eQTL significant pairs format.

    Derives position from the variant_id. ``gene`` matches as a substring, so
    either an ENSG ID or a gene symbol works. See the module docstring for the
    shared eQTL loader arguments and return value.

    Example:
        >>> eqtl_df = load_gtex_eqtl("GTEx_Analysis.signif_pairs.txt.gz", gene="BRCA1")
    """
    return _load_tabular(filepath, _GTEX_SPEC, gene=gene)


_EQTL_CATALOGUE_SPEC = LoaderSpec(
    log_fmt="Loaded eQTL Catalogue file with {n} associations",
    read={"sep": "\t"},
    col_map={
        "position": "pos",
        "pvalue": "p_value",
        "gene_id": "gene",
        "beta": "effect_size",
        "chromosome": "chr",
    },
    gene_filter="contains",
    schema=lambda out_cols: spec(Family.EQTL, Tier.LOAD),
)


def load_eqtl_catalogue(
    filepath: Union[str, Path],
    gene: Optional[str] = None,
) -> pd.DataFrame:
    """Load eQTL Catalogue format.

    See the module docstring for the shared eQTL loader arguments and return value.
    """
    return _load_tabular(filepath, _EQTL_CATALOGUE_SPEC, gene=gene)


_MATRIXEQTL_SPEC = LoaderSpec(
    log_fmt="Loaded MatrixEQTL file with {n} associations",
    read={"sep": "\t"},
    col_map={
        "SNP": "rs",
        "gene": "gene",
        "p-value": "p_value",
        "pvalue": "p_value",
        "beta": "effect_size",
        "t-stat": "t_stat",
    },
    gene_filter="exact",
    schema=lambda out_cols: spec(Family.EQTL, Tier.LOAD),
    # MatrixEQTL carries no position; callers merge a SNP annotation to add it.
    schema_requires=("pos",),
)


def load_matrixeqtl(
    filepath: Union[str, Path],
    gene: Optional[str] = None,
) -> pd.DataFrame:
    """Load MatrixEQTL output format.

    ``gene`` matches exactly here, unlike the substring match the other eQTL
    loaders use. MatrixEQTL carries no position, so merge a SNP annotation
    file to add ``pos`` before plotting. See the module docstring for the
    shared eQTL loader arguments and return value.
    """
    return _load_tabular(filepath, _MATRIXEQTL_SPEC, gene=gene)
