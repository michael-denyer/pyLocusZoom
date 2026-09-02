"""The column contract for every DataFrame family, at both validation tiers.

``Tier.LOAD`` is the strict tier a loader applies to a file it just parsed.
``Tier.PLOT`` is the permissive tier a plotter applies to a frame the caller
assembled. The split is deliberate: tightening the plot tier would reject
input that plots correctly today.
"""

from enum import Enum
from typing import Any, Callable, Dict, Optional

import pandas as pd

from .exceptions import (
    EQTLValidationError,
    FinemappingValidationError,
    ForestValidationError,
    LoaderValidationError,
    PheWASValidationError,
)
from .validation import ColumnSpec, RangeRule, check


class Canonical:
    """The column names every loader emits and every plotter defaults to.

    One vocabulary across the package: the GWAS, gene-annotation,
    fine-mapping and eQTL families all name their chromosome ``chr`` and
    their position ``pos``, and every family carrying a single p-value names
    it ``p_value``. A colocalization frame is the one exception, because it
    merges two p-value columns and cannot spell both of them the same way.
    """

    CHROM = "chr"
    POS = "pos"
    P = "p_value"
    RS = "rs"


# The pre-4.0 GEMMA spellings, keyed by the canonical name that replaced them.
# Read by the loaders (which no longer emit them) and by the plotters (which
# still accept a frame carrying them). Removed in DEPRECATED_ALIAS_REMOVED_IN.
DEPRECATED_COLUMN_ALIASES: Dict[str, str] = {
    Canonical.POS: "ps",
    Canonical.P: "p_wald",
}

DEPRECATED_ALIAS_REMOVED_IN = "5.0.0"


class Family(Enum):
    """A DataFrame family with a validation contract."""

    GWAS = "gwas"
    GENES = "genes"
    EQTL = "eqtl"
    FINEMAPPING = "finemapping"
    PHEWAS = "phewas"
    FOREST = "forest"
    COLOC = "coloc"


class Tier(Enum):
    """Strict load-time validation, or permissive plot-time validation."""

    LOAD = "load"
    PLOT = "plot"


def _gwas_load(
    pos_col: str = Canonical.POS,
    p_col: str = Canonical.P,
    rs_col: Optional[str] = None,
) -> ColumnSpec:
    return ColumnSpec(
        name="GWAS",
        required=(pos_col, p_col) if rs_col is None else (pos_col, p_col, rs_col),
        numeric=(pos_col, p_col),
        not_null=(pos_col, p_col),
        ranges=(RangeRule(pos_col, min_val=0, exclusive_min=True),),
        pvalue=p_col,
        error_class=LoaderValidationError,
    )


def _gwas_plot(
    pos_col: str = Canonical.POS,
    p_col: str = Canonical.P,
    rs_col: Optional[str] = None,
    chrom_col: Optional[str] = None,
) -> ColumnSpec:
    optional = tuple(col for col in (rs_col, chrom_col) if col is not None)
    return ColumnSpec(
        name="gwas_df",
        required=(pos_col, p_col, *optional),
        non_empty=True,
    )


def _genes_load() -> ColumnSpec:
    return ColumnSpec(
        name="Gene annotation",
        required=(Canonical.CHROM, "start", "end", "gene_name"),
        numeric=("start", "end"),
        not_null=("start", "end"),
        ranges=(RangeRule("start", min_val=0),),
        ordering=(("start", "end"),),
        error_class=LoaderValidationError,
    )


def _genes_plot() -> ColumnSpec:
    return ColumnSpec(
        name="genes_df", required=(Canonical.CHROM, "start", "end", "gene_name")
    )


def _eqtl_load() -> ColumnSpec:
    return ColumnSpec(
        name="eQTL",
        required=(Canonical.POS, Canonical.P, "gene"),
        numeric=(Canonical.POS, Canonical.P),
        not_null=(Canonical.POS, Canonical.P),
        ranges=(RangeRule(Canonical.POS, min_val=0, exclusive_min=True),),
        pvalue=Canonical.P,
        error_class=LoaderValidationError,
    )


def _eqtl_plot(pos_col: str = Canonical.POS, p_col: str = Canonical.P) -> ColumnSpec:
    return ColumnSpec(
        name="eQTL DataFrame",
        required=(pos_col, p_col),
        numeric=(p_col,),
        error_class=EQTLValidationError,
    )


def _finemapping_load() -> ColumnSpec:
    return ColumnSpec(
        name="Fine-mapping",
        required=(Canonical.POS, "pip"),
        numeric=(Canonical.POS, "pip"),
        not_null=(Canonical.POS, "pip"),
        ranges=(
            RangeRule(Canonical.POS, min_val=0, exclusive_min=True),
            RangeRule("pip", min_val=0, max_val=1),
        ),
        error_class=LoaderValidationError,
    )


def _finemapping_plot(pos_col: str = Canonical.POS, pip_col: str = "pip") -> ColumnSpec:
    return ColumnSpec(
        name="Fine-mapping DataFrame",
        required=(pos_col, pip_col),
        numeric=(pip_col,),
        ranges=(RangeRule(pip_col, min_val=0, max_val=1),),
        error_class=FinemappingValidationError,
    )


def _phewas_plot(
    phenotype_col: str = "phenotype", p_col: str = Canonical.P
) -> ColumnSpec:
    return ColumnSpec(
        name="PheWAS DataFrame",
        required=(phenotype_col, p_col),
        numeric=(p_col,),
        not_null=(p_col,),
        pvalue=p_col,
        error_class=PheWASValidationError,
    )


def _forest_plot(
    study_col: str = "study",
    effect_col: str = "effect",
    ci_lower_col: str = "ci_lower",
    ci_upper_col: str = "ci_upper",
) -> ColumnSpec:
    return ColumnSpec(
        name="Forest plot DataFrame",
        required=(study_col, effect_col, ci_lower_col, ci_upper_col),
        numeric=(effect_col, ci_lower_col, ci_upper_col),
        ordering=(
            (ci_lower_col, effect_col),
            (effect_col, ci_upper_col),
            (ci_lower_col, ci_upper_col),
        ),
        error_class=ForestValidationError,
    )


def _coloc_plot(
    name: str, pos_col: str, p_col: str, rs_col: Optional[str] = None
) -> ColumnSpec:
    required = (pos_col, p_col) if rs_col is None else (pos_col, p_col, rs_col)
    return ColumnSpec(
        name=name,
        required=required,
        numeric=(pos_col, p_col),
        not_null=(p_col,),
        pvalue=p_col,
    )


_SPECS: Dict[Family, Dict[Tier, Callable[..., ColumnSpec]]] = {
    Family.GWAS: {Tier.LOAD: _gwas_load, Tier.PLOT: _gwas_plot},
    Family.GENES: {Tier.LOAD: _genes_load, Tier.PLOT: _genes_plot},
    Family.EQTL: {Tier.LOAD: _eqtl_load, Tier.PLOT: _eqtl_plot},
    Family.FINEMAPPING: {
        Tier.LOAD: _finemapping_load,
        Tier.PLOT: _finemapping_plot,
    },
    Family.PHEWAS: {Tier.PLOT: _phewas_plot},
    Family.FOREST: {Tier.PLOT: _forest_plot},
    Family.COLOC: {Tier.PLOT: _coloc_plot},
}


def spec(family: Family, tier: Tier, **columns: Any) -> ColumnSpec:
    """Return the column contract for one family at one tier.

    Args:
        family: The DataFrame family.
        tier: ``Tier.LOAD`` for the strict loader contract, ``Tier.PLOT`` for
            the permissive plot-time one.
        **columns: Caller column names the contract is built over. Which are
            accepted depends on the family; a load-time contract fixes its own
            column names and takes none.

    Raises:
        KeyError: If the family has no contract at that tier.
    """
    tiers = _SPECS[family]
    if tier not in tiers:
        raise KeyError(f"{family.value} has no {tier.value}-time contract")
    return tiers[tier](**columns)


def validate_gwas_df(
    df: pd.DataFrame,
    pos_col: str = Canonical.POS,
    p_col: str = Canonical.P,
    rs_col: Optional[str] = None,
    chrom_col: Optional[str] = None,
) -> None:
    """Validate a GWAS DataFrame at the permissive plot-time tier.

    Args:
        df: GWAS results DataFrame.
        pos_col: Column name for position.
        p_col: Column name for p-values.
        rs_col: Column name for SNP IDs (optional).
        chrom_col: Column name for chromosome, required by the genome-wide
            families (optional).

    Raises:
        ValidationError: If required columns are missing or the frame is empty.
    """
    check(
        df,
        spec(
            Family.GWAS,
            Tier.PLOT,
            pos_col=pos_col,
            p_col=p_col,
            rs_col=rs_col,
            chrom_col=chrom_col,
        ),
    )


def validate_genes_df(df: pd.DataFrame) -> None:
    """Validate a gene annotations DataFrame at the plot-time tier.

    Args:
        df: Gene annotations DataFrame.

    Raises:
        ValidationError: If required columns are missing.
    """
    check(df, spec(Family.GENES, Tier.PLOT))


def validate_phewas_df(
    df: pd.DataFrame,
    phenotype_col: str = "phenotype",
    p_col: str = Canonical.P,
) -> None:
    """Validate PheWAS DataFrame has required columns and types.

    The category column is optional at render time, so it is not part of
    this contract.

    Args:
        df: PheWAS results DataFrame.
        phenotype_col: Column name for phenotype names.
        p_col: Column name for p-values.

    Raises:
        PheWASValidationError: If required columns are missing or have invalid types.
    """
    check(df, spec(Family.PHEWAS, Tier.PLOT, phenotype_col=phenotype_col, p_col=p_col))


def validate_forest_df(
    df: pd.DataFrame,
    study_col: str = "study",
    effect_col: str = "effect",
    ci_lower_col: str = "ci_lower",
    ci_upper_col: str = "ci_upper",
) -> None:
    """Validate forest plot DataFrame has required columns and types.

    Args:
        df: Forest plot data DataFrame.
        study_col: Column name for study/phenotype names.
        effect_col: Column name for effect sizes (beta, OR, HR).
        ci_lower_col: Column name for lower confidence interval.
        ci_upper_col: Column name for upper confidence interval.

    Raises:
        ForestValidationError: If required columns are missing or have invalid types.
    """
    check(
        df,
        spec(
            Family.FOREST,
            Tier.PLOT,
            study_col=study_col,
            effect_col=effect_col,
            ci_lower_col=ci_lower_col,
            ci_upper_col=ci_upper_col,
        ),
    )


def validate_coloc_df(
    df: pd.DataFrame,
    name: str,
    pos_col: str,
    p_col: str,
    rs_col: Optional[str] = None,
) -> None:
    """Validate a DataFrame for a colocalization plot.

    Args:
        df: Results DataFrame.
        name: Name for error messages (e.g. "GWAS DataFrame").
        pos_col: Column name for genomic positions.
        p_col: Column name for p-values.
        rs_col: Optional column name for SNP IDs.

    Raises:
        ValidationError: If required columns are missing or have invalid types.
    """
    check(
        df,
        spec(
            Family.COLOC,
            Tier.PLOT,
            name=name,
            pos_col=pos_col,
            p_col=p_col,
            rs_col=rs_col,
        ),
    )
