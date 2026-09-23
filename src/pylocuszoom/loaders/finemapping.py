"""Fine-mapping result loaders: SuSiE, FINEMAP, CAVIAR, and PolyFun.

Every fine-mapping loader takes ``filepath`` plus ``cs_col`` (output column
name for a source-supplied credible set, default "cs"). Loaders preserve PIPs
and any reported set membership; they do not infer credible sets. Position
and membership columns are included only when supplied by the source.

Raises:
    LoaderValidationError: If the format's columns cannot be mapped or the
        mapped values fail the format's schema.
"""

from pathlib import Path
from typing import Union

import pandas as pd

from ..schemas import FINEMAPPING_LOAD
from ._engine import LoaderSpec, _load_tabular


def _susie_cs(df: pd.DataFrame, out_cols: dict[str, str]) -> pd.DataFrame:
    """Standardize SuSiE credible set: -1 or NA (not in a set) become 0."""
    cs_col = out_cols["cs_col"]
    if cs_col in df.columns:
        df[cs_col] = df[cs_col].fillna(0).astype(int)
        df.loc[df[cs_col] < 0, cs_col] = 0
    return df


_SUSIE_SPEC = LoaderSpec(
    log_fmt="Loaded SuSiE file with {n} variants",
    read={"sep": "\t"},
    col_candidates={
        "pos": ("pos", "position", "BP", "bp", "POS"),
        "pip": ("pip", "PIP", "posterior_prob", "prob"),
        "cs_col": ("cs", "CS", "credible_set", "cs_index", "L"),
        "rs": ("snp", "SNP", "variant_id", "rsid"),
    },
    transform=_susie_cs,
    schema=lambda out_cols: FINEMAPPING_LOAD,
)


def load_susie(
    filepath: Union[str, Path],
    cs_col: str = "cs",
) -> pd.DataFrame:
    """Load SuSiE fine-mapping results.

    Supports both R susieR output (saved as TSV) and SuSiE-inf output. A
    credible set of -1 or NA (not in a set) becomes 0. See the module
    docstring for the shared fine-mapping arguments and return value.

    Example:
        >>> fm_df = load_susie("susie_results.tsv")
        >>> fig = plotter.plot_stacked([gwas_df], ..., finemapping_df=fm_df)
    """
    return _load_tabular(filepath, _SUSIE_SPEC, cs_col=cs_col)


_FINEMAP_SPEC = LoaderSpec(
    log_fmt="Loaded FINEMAP file with {n} variants",
    read={"sep": r"\s+"},
    col_map={
        "position": "pos",
        "prob": "pip",
        "rsid": "rs",
        "chromosome": "chr",
        "cs": "cs_col",
    },
    schema=lambda out_cols: FINEMAPPING_LOAD,
)


def load_finemap(
    filepath: Union[str, Path],
    cs_col: str = "cs",
) -> pd.DataFrame:
    """Load FINEMAP results (.snp output file).

    PIPs retain their input order. Credible-set membership is preserved if
    supplied, and otherwise omitted. See the module docstring for the shared
    fine-mapping arguments and return value.

    Example:
        >>> fm_df = load_finemap("results.snp")
    """
    return _load_tabular(filepath, _FINEMAP_SPEC, cs_col=cs_col)


_CAVIAR_SPEC = LoaderSpec(
    log_fmt="Loaded CAVIAR file with {n} variants",
    # CAVIAR .set files are headerless: SNP_ID Causal_Post_Prob.
    read={"sep": r"\s+", "header": None, "names": ["rs", "pip"]},
    schema=lambda out_cols: FINEMAPPING_LOAD,
    # CAVIAR carries no position; callers merge a SNP annotation to add it.
    schema_requires=("pos",),
)


def load_caviar(
    filepath: Union[str, Path],
    cs_col: str = "cs",
) -> pd.DataFrame:
    """Load CAVIAR results (.set output file).

    No credible set is inferred from the PIPs. CAVIAR carries no position, so
    merge a SNP annotation file to add ``pos`` before plotting.
    See the module docstring for the shared fine-mapping arguments and return
    value.
    """
    return _load_tabular(filepath, _CAVIAR_SPEC, cs_col=cs_col)


def _polyfun_cs(df: pd.DataFrame, out_cols: dict[str, str]) -> pd.DataFrame:
    """Normalize the credible-set column to non-null integers."""
    cs_col = out_cols["cs_col"]
    if cs_col in df.columns:
        df[cs_col] = df[cs_col].fillna(0).astype(int)
    return df


_POLYFUN_SPEC = LoaderSpec(
    log_fmt="Loaded PolyFun file with {n} variants",
    read={"sep": r"\s+"},
    col_map={
        "BP": "pos",
        "PIP": "pip",
        "SNP": "rs",
        "CHR": "chr",
        "CREDIBLE_SET": "cs_col",
    },
    transform=_polyfun_cs,
    schema=lambda out_cols: FINEMAPPING_LOAD,
)


def load_polyfun(
    filepath: Union[str, Path],
    cs_col: str = "cs",
) -> pd.DataFrame:
    """Load PolyFun/SuSiE fine-mapping results.

    See the module docstring for the shared fine-mapping arguments and return value.
    """
    return _load_tabular(filepath, _POLYFUN_SPEC, cs_col=cs_col)
