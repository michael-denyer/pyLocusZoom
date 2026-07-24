"""File format loaders for common GWAS, eQTL, and fine-mapping outputs.

Convenience functions to load data from standard file formats into
DataFrames ready for use with LocusZoomPlotter.

GWAS formats:
- PLINK (.assoc, .assoc.linear, .assoc.logistic, .qassoc)
- REGENIE (.regenie)
- BOLT-LMM (.stats)
- GEMMA (.assoc.txt)
- SAIGE (.txt)
- Generic TSV/CSV

eQTL formats:
- GTEx significant pairs format
- eQTL Catalogue format
- MatrixEQTL output

Fine-mapping formats:
- SuSiE (susieR output)
- FINEMAP (.snp output)
- CAVIAR (.set output)

Gene annotation formats:
- GTF/GFF3
- BED (4-column: chr, start, end, name)
"""

from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from .logging import logger
from .schemas import (
    LoaderValidationError,
    validate_eqtl_dataframe,
    validate_finemapping_dataframe,
    validate_genes_dataframe,
    validate_gwas_dataframe,
)


def _validate_or_warn(
    df: pd.DataFrame,
    required: list[str],
    loader_name: str,
    validate_fn,
    **validate_kwargs,
) -> None:
    """Run validation if required columns are present, otherwise warn."""
    if all(col in df.columns for col in required):
        validate_fn(df, **validate_kwargs)
    else:
        missing = [c for c in required if c not in df.columns]
        logger.warning(
            f"{loader_name} loader could not map columns: {missing}. "
            f"Validation skipped. Available columns: {list(df.columns)}"
        )


# =============================================================================
# Table-driven loader engine
# =============================================================================
#
# Most loaders share one shape: read a delimited file, rename source columns to
# standard names, log a count, then validate. `LoaderSpec` captures the parts
# that differ between formats, and `_load_tabular` is the single engine that
# runs the shared steps. A new static format is one `LoaderSpec` constant plus a
# thin wrapper, not another copy-pasted function body.
#
# Output-column tokens: a spec target equal to one of the `**out_cols` keys
# passed by the wrapper (`pos_col`, `p_col`, `rs_col`, `cs_col`) is substituted
# with the caller's chosen name; any other target is used literally. This is how
# the caller-configurable column names flow into an otherwise constant spec.


@dataclass(frozen=True)
class LoaderSpec:
    """Declarative description of how to load one tabular file format.

    Attributes:
        sep: Field separator passed to ``pandas.read_csv``.
        log_fmt: Debug log template; ``{n}`` is filled with the row count.
        validate: Callable ``(df, out_cols)`` that runs the format's validation.
        comment: Comment-line prefix for ``read_csv`` (e.g. ``"#"``).
        col_map: Static ``source -> target`` renames. Targets may be tokens.
        p_candidates: P-value source columns; the first present maps to ``p_col``.
        col_candidates: ``target -> candidates``; first present candidate maps to
            the (possibly token) target. Replaces hand-rolled first-match loops.
        transform: Optional ``(df, out_cols) -> df`` applied after renaming, for
            format-specific reshaping (e.g. credible-set assignment).
        gene_filter: ``"contains"`` or ``"exact"`` to enable eQTL gene filtering.
    """

    sep: str
    log_fmt: str
    validate: Callable[[pd.DataFrame, dict[str, str]], None]
    comment: Optional[str] = None
    col_map: dict[str, str] = field(default_factory=dict)
    p_candidates: tuple[str, ...] = ()
    col_candidates: dict[str, tuple[str, ...]] = field(default_factory=dict)
    transform: Optional[Callable[[pd.DataFrame, dict[str, str]], pd.DataFrame]] = None
    gene_filter: Optional[str] = None


def _resolve(target: str, out_cols: dict[str, str]) -> str:
    """Substitute an output-column token, or return the literal target."""
    return out_cols.get(target, target)


def _first_present(df: pd.DataFrame, candidates: tuple[str, ...]) -> Optional[str]:
    """Return the first candidate column present in ``df``, else None."""
    for col in candidates:
        if col in df.columns:
            return col
    return None


def _filter_gene(df: pd.DataFrame, gene: str, mode: str) -> pd.DataFrame:
    """Filter eQTL rows to a gene by substring ("contains") or exact match."""
    if mode == "contains":
        return df[df["gene"].str.contains(gene, case=False, na=False)]
    return df[df["gene"] == gene]


def _load_tabular(
    filepath: Union[str, Path],
    spec: LoaderSpec,
    *,
    gene: Optional[str] = None,
    **out_cols: str,
) -> pd.DataFrame:
    """Load a tabular file per ``spec``: read, map, rename, transform, validate."""
    df = pd.read_csv(filepath, sep=spec.sep, comment=spec.comment)

    col_map = {src: _resolve(dst, out_cols) for src, dst in spec.col_map.items()}

    match = _first_present(df, spec.p_candidates)
    if match is not None:
        col_map[match] = out_cols["p_col"]

    for target, candidates in spec.col_candidates.items():
        match = _first_present(df, candidates)
        if match is not None:
            col_map[match] = _resolve(target, out_cols)

    df = df.rename(columns=col_map)

    if spec.transform is not None:
        df = spec.transform(df, out_cols)

    if spec.gene_filter is not None and gene is not None and "gene" in df.columns:
        df = _filter_gene(df, gene, spec.gene_filter)

    logger.debug(spec.log_fmt.format(n=len(df)))
    spec.validate(df, out_cols)
    return df


def _validate_gwas(df: pd.DataFrame, out_cols: dict[str, str]) -> None:
    """Strict GWAS validation using the caller's position/p-value columns."""
    validate_gwas_dataframe(df, pos_col=out_cols["pos_col"], p_col=out_cols["p_col"])


def _no_validate(df: pd.DataFrame, out_cols: dict[str, str]) -> None:
    """No-op validator for formats that cannot be fully validated on load."""


def _warn_validator(
    required: list[str],
    loader_name: str,
    validate_fn,
) -> Callable[[pd.DataFrame, dict[str, str]], None]:
    """Build a validator that validates when required columns are present, else warns."""

    def _validate(df: pd.DataFrame, out_cols: dict[str, str]) -> None:
        _validate_or_warn(df, required, loader_name, validate_fn)

    return _validate


# =============================================================================
# GWAS Loaders
# =============================================================================


_PLINK_SPEC = LoaderSpec(
    sep=r"\s+",
    comment="#",
    log_fmt="Loaded PLINK file with {n} variants",
    col_candidates={
        "pos_col": ("BP", "POS", "bp", "pos"),
        "p_col": ("P", "P_BOLT_LMM", "p", "PVAL", "pval", "P_LINREG"),
        "rs_col": ("SNP", "ID", "rsid", "RSID", "MarkerName", "variant_id"),
        "chr": ("CHR", "chr", "CHROM", "chrom", "#CHROM"),
    },
    validate=_validate_gwas,
)


def load_plink_assoc(
    filepath: Union[str, Path],
    pos_col: str = "ps",
    p_col: str = "p_wald",
    rs_col: str = "rs",
) -> pd.DataFrame:
    """Load PLINK association results (.assoc, .assoc.linear, .assoc.logistic, .qassoc).

    Automatically detects PLINK format variant and maps columns to standard names.

    Args:
        filepath: Path to PLINK association file.
        pos_col: Output column name for position. Default "ps".
        p_col: Output column name for p-value. Default "p_wald".
        rs_col: Output column name for SNP ID. Default "rs".

    Returns:
        DataFrame with standardized column names.

    Example:
        >>> gwas_df = load_plink_assoc("results.assoc.linear")
        >>> fig = plotter.plot(gwas_df, chrom=1, start=1e6, end=2e6)
    """
    return _load_tabular(
        filepath, _PLINK_SPEC, pos_col=pos_col, p_col=p_col, rs_col=rs_col
    )


def _regenie_pvalue(df: pd.DataFrame, out_cols: dict[str, str]) -> pd.DataFrame:
    """Derive REGENIE p-values: prefer computed LOG10P, else rename P."""
    p_col = out_cols["p_col"]
    if "LOG10P" in df.columns:
        df[p_col] = 10 ** (-df["LOG10P"])
    elif "P" in df.columns:
        df = df.rename(columns={"P": p_col})
    return df


_REGENIE_SPEC = LoaderSpec(
    sep=r"\s+",
    comment="#",
    log_fmt="Loaded REGENIE file with {n} variants",
    col_map={"GENPOS": "pos_col", "ID": "rs_col", "CHROM": "chr"},
    transform=_regenie_pvalue,
    validate=_validate_gwas,
)


def load_regenie(
    filepath: Union[str, Path],
    pos_col: str = "ps",
    p_col: str = "p_wald",
    rs_col: str = "rs",
) -> pd.DataFrame:
    """Load REGENIE association results (.regenie).

    Args:
        filepath: Path to REGENIE results file.
        pos_col: Output column name for position. Default "ps".
        p_col: Output column name for p-value. Default "p_wald".
        rs_col: Output column name for SNP ID. Default "rs".

    Returns:
        DataFrame with standardized column names.

    Example:
        >>> gwas_df = load_regenie("results.regenie")
    """
    return _load_tabular(
        filepath, _REGENIE_SPEC, pos_col=pos_col, p_col=p_col, rs_col=rs_col
    )


_BOLT_LMM_SPEC = LoaderSpec(
    sep="\t",
    log_fmt="Loaded BOLT-LMM file with {n} variants",
    col_map={"BP": "pos_col", "SNP": "rs_col", "CHR": "chr"},
    p_candidates=("P_BOLT_LMM", "P_BOLT_LMM_INF"),
    validate=_validate_gwas,
)


def load_bolt_lmm(
    filepath: Union[str, Path],
    pos_col: str = "ps",
    p_col: str = "p_wald",
    rs_col: str = "rs",
) -> pd.DataFrame:
    """Load BOLT-LMM association results (.stats).

    Args:
        filepath: Path to BOLT-LMM stats file.
        pos_col: Output column name for position. Default "ps".
        p_col: Output column name for p-value. Default "p_wald".
        rs_col: Output column name for SNP ID. Default "rs".

    Returns:
        DataFrame with standardized column names.

    Example:
        >>> gwas_df = load_bolt_lmm("results.stats")
    """
    return _load_tabular(
        filepath, _BOLT_LMM_SPEC, pos_col=pos_col, p_col=p_col, rs_col=rs_col
    )


_GEMMA_SPEC = LoaderSpec(
    sep="\t",
    log_fmt="Loaded GEMMA file with {n} variants",
    col_map={"ps": "pos_col", "rs": "rs_col", "chr": "chr"},
    p_candidates=("p_wald", "p_lrt", "p_score"),
    validate=_validate_gwas,
)


def load_gemma(
    filepath: Union[str, Path],
    pos_col: str = "ps",
    p_col: str = "p_wald",
    rs_col: str = "rs",
) -> pd.DataFrame:
    """Load GEMMA association results (.assoc.txt).

    Args:
        filepath: Path to GEMMA association file.
        pos_col: Output column name for position. Default "ps".
        p_col: Output column name for p-value. Default "p_wald".
        rs_col: Output column name for SNP ID. Default "rs".

    Returns:
        DataFrame with standardized column names.

    Example:
        >>> gwas_df = load_gemma("output.assoc.txt")
    """
    return _load_tabular(
        filepath, _GEMMA_SPEC, pos_col=pos_col, p_col=p_col, rs_col=rs_col
    )


_SAIGE_SPEC = LoaderSpec(
    sep="\t",
    log_fmt="Loaded SAIGE file with {n} variants",
    col_map={"POS": "pos_col", "MarkerID": "rs_col", "CHR": "chr"},
    p_candidates=("p.value.NA", "p.value"),
    validate=_validate_gwas,
)


def load_saige(
    filepath: Union[str, Path],
    pos_col: str = "ps",
    p_col: str = "p_wald",
    rs_col: str = "rs",
) -> pd.DataFrame:
    """Load SAIGE association results.

    Args:
        filepath: Path to SAIGE results file.
        pos_col: Output column name for position. Default "ps".
        p_col: Output column name for p-value. Default "p_wald".
        rs_col: Output column name for SNP ID. Default "rs".

    Returns:
        DataFrame with standardized column names.

    Example:
        >>> gwas_df = load_saige("results.txt")
    """
    return _load_tabular(
        filepath, _SAIGE_SPEC, pos_col=pos_col, p_col=p_col, rs_col=rs_col
    )


_GWAS_CATALOG_SPEC = LoaderSpec(
    sep="\t",
    log_fmt="Loaded GWAS Catalog file with {n} variants",
    col_map={
        "base_pair_location": "pos_col",
        "variant_id": "rs_col",
        "chromosome": "chr",
        "p_value": "p_col",
    },
    validate=_validate_gwas,
)


def load_gwas_catalog(
    filepath: Union[str, Path],
    pos_col: str = "ps",
    p_col: str = "p_wald",
    rs_col: str = "rs",
) -> pd.DataFrame:
    """Load GWAS Catalog summary statistics format.

    Args:
        filepath: Path to GWAS Catalog file.
        pos_col: Output column name for position. Default "ps".
        p_col: Output column name for p-value. Default "p_wald".
        rs_col: Output column name for SNP ID. Default "rs".

    Returns:
        DataFrame with standardized column names.
    """
    return _load_tabular(
        filepath, _GWAS_CATALOG_SPEC, pos_col=pos_col, p_col=p_col, rs_col=rs_col
    )


# =============================================================================
# eQTL Loaders
# =============================================================================


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
    sep="\t",
    log_fmt="Loaded GTEx eQTL file with {n} associations",
    col_candidates={
        "p_value": ("pval_nominal", "p_value", "pvalue", "P"),
        "gene": ("gene_id", "gene_name", "phenotype_id"),
        "effect_size": ("slope", "beta", "effect_size"),
    },
    transform=_gtex_pos,
    gene_filter="contains",
    validate=_warn_validator(
        ["pos", "p_value", "gene"], "GTEx eQTL", validate_eqtl_dataframe
    ),
)


def load_gtex_eqtl(
    filepath: Union[str, Path],
    gene: Optional[str] = None,
) -> pd.DataFrame:
    """Load GTEx eQTL significant pairs format.

    Args:
        filepath: Path to GTEx eQTL file (e.g., signif_variant_gene_pairs.txt.gz).
        gene: Optional gene to filter to (ENSG ID or gene symbol).

    Returns:
        DataFrame with columns: pos, p_value, gene, effect_size.

    Example:
        >>> eqtl_df = load_gtex_eqtl("GTEx_Analysis.signif_pairs.txt.gz", gene="BRCA1")
    """
    return _load_tabular(filepath, _GTEX_SPEC, gene=gene)


_EQTL_CATALOGUE_SPEC = LoaderSpec(
    sep="\t",
    log_fmt="Loaded eQTL Catalogue file with {n} associations",
    col_map={
        "position": "pos",
        "pvalue": "p_value",
        "gene_id": "gene",
        "beta": "effect_size",
        "chromosome": "chr",
    },
    gene_filter="contains",
    validate=_warn_validator(
        ["pos", "p_value", "gene"], "eQTL Catalogue", validate_eqtl_dataframe
    ),
)


def load_eqtl_catalogue(
    filepath: Union[str, Path],
    gene: Optional[str] = None,
) -> pd.DataFrame:
    """Load eQTL Catalogue format.

    Args:
        filepath: Path to eQTL Catalogue file.
        gene: Optional gene to filter to.

    Returns:
        DataFrame with columns: pos, p_value, gene, effect_size.
    """
    return _load_tabular(filepath, _EQTL_CATALOGUE_SPEC, gene=gene)


_MATRIXEQTL_SPEC = LoaderSpec(
    sep="\t",
    log_fmt="Loaded MatrixEQTL file with {n} associations",
    col_map={
        "SNP": "rs",
        "gene": "gene",
        "p-value": "p_value",
        "pvalue": "p_value",
        "beta": "effect_size",
        "t-stat": "t_stat",
    },
    gene_filter="exact",
    validate=_no_validate,
)


def load_matrixeqtl(
    filepath: Union[str, Path],
    gene: Optional[str] = None,
) -> pd.DataFrame:
    """Load MatrixEQTL output format.

    Args:
        filepath: Path to MatrixEQTL output file.
        gene: Optional gene to filter to.

    Returns:
        DataFrame with columns: pos, p_value, gene, effect_size.

    Note:
        MatrixEQTL output doesn't include position by default.
        You may need to merge with a SNP annotation file.
    """
    return _load_tabular(filepath, _MATRIXEQTL_SPEC, gene=gene)


# =============================================================================
# Fine-mapping Loaders
# =============================================================================


def _susie_cs(df: pd.DataFrame, out_cols: dict[str, str]) -> pd.DataFrame:
    """Standardize SuSiE credible set: -1 or NA (not in a set) become 0."""
    cs_col = out_cols["cs_col"]
    if cs_col in df.columns:
        df[cs_col] = df[cs_col].fillna(0).astype(int)
        df.loc[df[cs_col] < 0, cs_col] = 0
    return df


_SUSIE_SPEC = LoaderSpec(
    sep="\t",
    log_fmt="Loaded SuSiE file with {n} variants",
    col_candidates={
        "pos": ("pos", "position", "BP", "bp", "POS"),
        "pip": ("pip", "PIP", "posterior_prob", "prob"),
        "cs_col": ("cs", "CS", "credible_set", "cs_index", "L"),
        "rs": ("snp", "SNP", "variant_id", "rsid"),
    },
    transform=_susie_cs,
    validate=_warn_validator(["pos", "pip"], "SuSiE", validate_finemapping_dataframe),
)


def load_susie(
    filepath: Union[str, Path],
    cs_col: str = "cs",
) -> pd.DataFrame:
    """Load SuSiE fine-mapping results.

    Supports both R susieR output (saved as TSV) and SuSiE-inf output.

    Args:
        filepath: Path to SuSiE results file.
        cs_col: Output column name for credible set. Default "cs".

    Returns:
        DataFrame with columns: pos, pip, cs.

    Example:
        >>> fm_df = load_susie("susie_results.tsv")
        >>> fig = plotter.plot_stacked([gwas_df], ..., finemapping_df=fm_df)
    """
    return _load_tabular(filepath, _SUSIE_SPEC, cs_col=cs_col)


def _finemap_cs(df: pd.DataFrame, out_cols: dict[str, str]) -> pd.DataFrame:
    """Assign a 95% credible set from cumulative PIP (FINEMAP has none)."""
    cs_col = out_cols["cs_col"]
    if cs_col not in df.columns and "pip" in df.columns:
        df = df.sort_values("pip", ascending=False)
        df["cumsum_pip"] = df["pip"].cumsum()
        df[cs_col] = (df["cumsum_pip"] <= 0.95).astype(int)
        df = df.drop(columns=["cumsum_pip"])
    return df


_FINEMAP_SPEC = LoaderSpec(
    sep=r"\s+",
    log_fmt="Loaded FINEMAP file with {n} variants",
    col_map={"position": "pos", "prob": "pip", "rsid": "rs", "chromosome": "chr"},
    transform=_finemap_cs,
    validate=_warn_validator(["pos", "pip"], "FINEMAP", validate_finemapping_dataframe),
)


def load_finemap(
    filepath: Union[str, Path],
    cs_col: str = "cs",
) -> pd.DataFrame:
    """Load FINEMAP results (.snp output file).

    Args:
        filepath: Path to FINEMAP .snp output file.
        cs_col: Output column name for credible set. Default "cs".

    Returns:
        DataFrame with columns: pos, pip, cs.

    Example:
        >>> fm_df = load_finemap("results.snp")
    """
    return _load_tabular(filepath, _FINEMAP_SPEC, cs_col=cs_col)


def load_caviar(
    filepath: Union[str, Path],
    cs_col: str = "cs",
) -> pd.DataFrame:
    """Load CAVIAR results (.set output file).

    Args:
        filepath: Path to CAVIAR output file.
        cs_col: Output column name for credible set. Default "cs".

    Returns:
        DataFrame with columns: pos, pip, cs.
    """
    # CAVIAR .set file format: SNP_ID Causal_Post_Prob
    df = pd.read_csv(filepath, sep=r"\s+", header=None, names=["rs", "pip"])

    # CAVIAR doesn't include position - user needs to merge
    logger.warning(
        "CAVIAR output doesn't include positions. "
        "Merge with SNP annotation file to add 'pos' column."
    )

    # Assign credible set based on PIP threshold
    df = df.sort_values("pip", ascending=False)
    df["cumsum_pip"] = df["pip"].cumsum()
    df[cs_col] = (df["cumsum_pip"] <= 0.95).astype(int)
    df = df.drop(columns=["cumsum_pip"])

    logger.debug(f"Loaded CAVIAR file with {len(df)} variants")

    # CAVIAR doesn't have pos - can't fully validate
    if "pip" in df.columns:
        if ((df["pip"] < 0) | (df["pip"] > 1)).any():
            raise LoaderValidationError("PIP values must be in range [0, 1]")

    return df


def _polyfun_cs(df: pd.DataFrame, out_cols: dict[str, str]) -> pd.DataFrame:
    """Normalize the credible-set column to non-null integers."""
    cs_col = out_cols["cs_col"]
    if cs_col in df.columns:
        df[cs_col] = df[cs_col].fillna(0).astype(int)
    return df


_POLYFUN_SPEC = LoaderSpec(
    sep=r"\s+",
    log_fmt="Loaded PolyFun file with {n} variants",
    col_map={
        "BP": "pos",
        "PIP": "pip",
        "SNP": "rs",
        "CHR": "chr",
        "CREDIBLE_SET": "cs_col",
    },
    transform=_polyfun_cs,
    validate=_warn_validator(["pos", "pip"], "PolyFun", validate_finemapping_dataframe),
)


def load_polyfun(
    filepath: Union[str, Path],
    cs_col: str = "cs",
) -> pd.DataFrame:
    """Load PolyFun/SuSiE fine-mapping results.

    Args:
        filepath: Path to PolyFun output file.
        cs_col: Output column name for credible set. Default "cs".

    Returns:
        DataFrame with columns: pos, pip, cs.
    """
    return _load_tabular(filepath, _POLYFUN_SPEC, cs_col=cs_col)


# =============================================================================
# Gene Annotation Loaders
# =============================================================================


def load_gtf(
    filepath: Union[str, Path],
    feature_type: str = "gene",
) -> pd.DataFrame:
    """Load gene annotations from GTF/GFF3 file.

    Args:
        filepath: Path to GTF or GFF3 file (can be gzipped).
        feature_type: Feature type to extract ("gene", "exon", "transcript").
            Default "gene".

    Returns:
        DataFrame with columns: chr, start, end, gene_name, strand.

    Example:
        >>> genes_df = load_gtf("genes.gtf", feature_type="gene")
        >>> exons_df = load_gtf("genes.gtf", feature_type="exon")
    """
    # GTF columns: seqname, source, feature, start, end, score, strand, frame, attributes
    df = pd.read_csv(
        filepath,
        sep="\t",
        comment="#",
        header=None,
        names=[
            "chr",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attributes",
        ],
    )

    # Filter to requested feature type
    df = df[df["feature"] == feature_type].copy()

    # Parse gene_name from attributes
    def extract_gene_name(attrs: str) -> str:
        """Extract gene_name or gene_id from GTF attributes."""
        for attr in attrs.split(";"):
            attr = attr.strip()
            if attr.startswith("gene_name"):
                # gene_name "BRCA1" or gene_name=BRCA1
                return attr.split('"')[1] if '"' in attr else attr.split("=")[1]
            if attr.startswith("gene_id"):
                return attr.split('"')[1] if '"' in attr else attr.split("=")[1]
        return ""

    df["gene_name"] = df["attributes"].apply(extract_gene_name)

    # Clean chromosome names
    df["chr"] = df["chr"].astype(str).str.replace("chr", "", regex=False)

    # Select and return relevant columns
    result = df[["chr", "start", "end", "gene_name", "strand"]].copy()
    logger.debug(f"Loaded {len(result)} {feature_type} features from GTF")
    validate_genes_dataframe(result)
    return result


def load_bed(
    filepath: Union[str, Path],
    has_header: bool = False,
) -> pd.DataFrame:
    """Load gene annotations from BED file.

    Supports BED4+ format (chr, start, end, name, ...).

    Args:
        filepath: Path to BED file.
        has_header: Whether file has header row. Default False.

    Returns:
        DataFrame with columns: chr, start, end, gene_name.

    Example:
        >>> genes_df = load_bed("genes.bed")
    """
    header = 0 if has_header else None
    df = pd.read_csv(filepath, sep="\t", header=header)

    # Assign column names if no header
    if not has_header:
        n_cols = len(df.columns)
        # Standard BED column names (up to BED12)
        bed_col_names = [
            "chr",
            "start",
            "end",
            "gene_name",
            "score",
            "strand",
            "thickStart",
            "thickEnd",
            "itemRgb",
            "blockCount",
            "blockSizes",
            "blockStarts",
        ]
        # Use standard names for known columns, generic for extras
        if n_cols <= len(bed_col_names):
            df.columns = bed_col_names[:n_cols]
        else:
            # More columns than BED12 - use known names + generic
            extra_cols = [f"col{i}" for i in range(len(bed_col_names), n_cols)]
            df.columns = bed_col_names + extra_cols

    # Standardize column names if header was present
    col_map = {
        "chrom": "chr",
        "chromStart": "start",
        "chromEnd": "end",
        "name": "gene_name",
    }
    df = df.rename(columns=col_map)

    # Clean chromosome names
    if "chr" in df.columns:
        df["chr"] = df["chr"].astype(str).str.replace("chr", "", regex=False)

    logger.debug(f"Loaded {len(df)} features from BED")

    _validate_or_warn(
        df, ["chr", "start", "end", "gene_name"], "BED", validate_genes_dataframe
    )

    return df


def _ensembl_strand(df: pd.DataFrame, out_cols: dict[str, str]) -> pd.DataFrame:
    """Convert Ensembl 1/-1 strand encoding to +/- symbols."""
    if "strand" in df.columns:
        df["strand"] = df["strand"].map({1: "+", -1: "-", "+": "+", "-": "-"})
    return df


_ENSEMBL_SPEC = LoaderSpec(
    sep="\t",
    log_fmt="Loaded {n} genes from Ensembl export",
    col_map={
        "Chromosome/scaffold name": "chr",
        "Gene start (bp)": "start",
        "Gene end (bp)": "end",
        "Gene name": "gene_name",
        "Strand": "strand",
        # Alternative column names
        "chromosome_name": "chr",
        "start_position": "start",
        "end_position": "end",
        "external_gene_name": "gene_name",
    },
    transform=_ensembl_strand,
    validate=_warn_validator(
        ["chr", "start", "end", "gene_name"], "Ensembl genes", validate_genes_dataframe
    ),
)


def load_ensembl_genes(
    filepath: Union[str, Path],
) -> pd.DataFrame:
    """Load Ensembl BioMart gene export.

    Args:
        filepath: Path to BioMart export file (TSV).

    Returns:
        DataFrame with columns: chr, start, end, gene_name, strand.
    """
    return _load_tabular(filepath, _ENSEMBL_SPEC)


# =============================================================================
# Generic Loader
# =============================================================================


def load_gwas(
    filepath: Union[str, Path],
    format: Optional[str] = None,
    pos_col: str = "ps",
    p_col: str = "p_wald",
    rs_col: str = "rs",
    **kwargs,
) -> pd.DataFrame:
    """Load GWAS results with automatic format detection.

    Args:
        filepath: Path to GWAS results file.
        format: File format. If None, auto-detects from extension.
            Options: "plink", "regenie", "bolt", "gemma", "saige", "catalog".
        pos_col: Output column name for position. Default "ps".
        p_col: Output column name for p-value. Default "p_wald".
        rs_col: Output column name for SNP ID. Default "rs".
        **kwargs: Additional arguments passed to format-specific loader.

    Returns:
        DataFrame with standardized column names.

    Example:
        >>> # Auto-detect format
        >>> gwas_df = load_gwas("results.assoc.linear")
        >>>
        >>> # Explicit format
        >>> gwas_df = load_gwas("results.txt", format="regenie")
    """
    filepath = Path(filepath)
    name = filepath.name.lower()

    # Auto-detect format from filename
    if format is None:
        if ".assoc" in name or ".qassoc" in name:
            format = "plink"
        elif ".regenie" in name:
            format = "regenie"
        elif ".stats" in name:
            format = "bolt"
        elif "gemma" in name or name.endswith(".assoc.txt"):
            format = "gemma"
        elif "saige" in name:
            format = "saige"
        else:
            logger.warning(
                "Could not detect GWAS file format for '%s'. "
                "Defaulting to 'plink'. Specify format= parameter explicitly.",
                filepath,
            )
            format = "plink"

    loaders = {
        "plink": load_plink_assoc,
        "regenie": load_regenie,
        "bolt": load_bolt_lmm,
        "gemma": load_gemma,
        "saige": load_saige,
        "catalog": load_gwas_catalog,
    }

    if format not in loaders:
        raise ValueError(f"Unknown format '{format}'. Options: {list(loaders.keys())}")

    return loaders[format](
        filepath, pos_col=pos_col, p_col=p_col, rs_col=rs_col, **kwargs
    )
