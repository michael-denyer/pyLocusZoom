"""The column contract for every DataFrame family, as ``ColumnSpec`` values.

Each contract is a value a caller hands to :func:`validation.check`: a
constant where the column names are fixed, a builder where they are the
caller's. The ``*_LOAD`` contracts are the strict tier a loader applies to a
file it just parsed; the ``*_plot_spec`` builders are the permissive tier a
plotter applies to a frame the caller assembled. The split is deliberate:
tightening the plot tier would reject input that plots correctly today.
"""

from typing import Dict, Optional

from .exceptions import (
    EQTLValidationError,
    FinemappingValidationError,
    ForestValidationError,
    LoaderValidationError,
    PheWASValidationError,
)
from .validation import ColumnSpec, RangeRule


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


def gwas_load_spec(
    pos_col: str = Canonical.POS, p_col: str = Canonical.P
) -> ColumnSpec:
    """Return the strict contract for a GWAS file a loader just parsed."""
    return ColumnSpec(
        name="GWAS",
        required=(pos_col, p_col),
        numeric=(pos_col,),
        not_null=(pos_col,),
        ranges=(RangeRule(pos_col, min_val=0, exclusive_min=True),),
        pvalue=p_col,
        error_class=LoaderValidationError,
    )


def gwas_plot_spec(
    pos_col: str = Canonical.POS,
    p_col: str = Canonical.P,
    rs_col: Optional[str] = None,
    chrom_col: Optional[str] = None,
) -> ColumnSpec:
    """Return the plot-time contract for a GWAS frame the caller assembled.

    Args:
        pos_col: Column name for position.
        p_col: Column name for p-values.
        rs_col: Column name for SNP ids, required when given.
        chrom_col: Column name for chromosome, required when given.
    """
    optional = tuple(col for col in (rs_col, chrom_col) if col is not None)
    return ColumnSpec(
        name="gwas_df",
        required=(pos_col, p_col, *optional),
        non_empty=True,
    )


GENES_LOAD = ColumnSpec(
    name="Gene annotation",
    required=(Canonical.CHROM, "start", "end", "gene_name"),
    numeric=("start", "end"),
    not_null=("start", "end"),
    ranges=(RangeRule("start", min_val=0),),
    ordering=(("start", "end"),),
    error_class=LoaderValidationError,
)

GENES_PLOT = ColumnSpec(
    name="genes_df", required=(Canonical.CHROM, "start", "end", "gene_name")
)

EXONS_PLOT = ColumnSpec(
    name="exons_df", required=(Canonical.CHROM, "start", "end", "gene_name")
)

EQTL_LOAD = ColumnSpec(
    name="eQTL",
    required=(Canonical.POS, Canonical.P, "gene"),
    numeric=(Canonical.POS,),
    not_null=(Canonical.POS,),
    ranges=(RangeRule(Canonical.POS, min_val=0, exclusive_min=True),),
    pvalue=Canonical.P,
    error_class=LoaderValidationError,
)


def eqtl_plot_spec(
    pos_col: str = Canonical.POS, p_col: str = Canonical.P
) -> ColumnSpec:
    """Return the plot-time contract for an eQTL frame."""
    return ColumnSpec(
        name="eQTL DataFrame",
        required=(pos_col, p_col),
        error_class=EQTLValidationError,
    )


FINEMAPPING_LOAD = ColumnSpec(
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


def finemapping_plot_spec(
    pos_col: str = Canonical.POS, pip_col: str = "pip"
) -> ColumnSpec:
    """Return the plot-time contract for a fine-mapping frame."""
    return ColumnSpec(
        name="Fine-mapping DataFrame",
        required=(pos_col, pip_col),
        numeric=(pip_col,),
        ranges=(RangeRule(pip_col, min_val=0, max_val=1),),
        error_class=FinemappingValidationError,
    )


def phewas_plot_spec(
    phenotype_col: str = "phenotype", p_col: str = Canonical.P
) -> ColumnSpec:
    """Return the contract for a PheWAS frame.

    The p-values themselves are checked by the ``"phewas"`` row of
    ``_data.P_VALUE_POLICY``, and the category column is optional, so this
    contract rules only on the two columns' presence.
    """
    return ColumnSpec(
        name="PheWAS DataFrame",
        required=(phenotype_col, p_col),
        error_class=PheWASValidationError,
    )


def forest_plot_spec(
    study_col: str = "study",
    effect_col: str = "effect",
    ci_lower_col: str = "ci_lower",
    ci_upper_col: str = "ci_upper",
) -> ColumnSpec:
    """Return the contract for a forest-plot frame."""
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


def coloc_plot_spec(name: str, pos_col: str, p_col: str) -> ColumnSpec:
    """Return the contract for one colocalization source frame.

    Args:
        name: The frame's name in error messages, e.g. "GWAS DataFrame".
        pos_col: Column name for genomic positions.
        p_col: Column name for p-values, whose values the ``"coloc"`` row of
            ``_data.P_VALUE_POLICY`` checks.
    """
    return ColumnSpec(name=name, required=(pos_col, p_col), numeric=(pos_col,))
