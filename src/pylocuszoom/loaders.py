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

Every ``load_*`` GWAS loader takes the same arguments and returns the same
shape, so the per-loader docstrings below cover only what is specific to their
format:

Args:
    filepath: Path to the results file.
    pos_col: Output column name for position. Default "ps".
    p_col: Output column name for p-value. Default "p_wald".
    rs_col: Output column name for SNP ID. Default "rs".

Returns:
    DataFrame with standardized column names.

Raises:
    LoaderValidationError: If the format's columns cannot be mapped or the
        mapped values fail the format's schema.

The eQTL loaders take ``filepath`` plus an optional ``gene`` to filter to, and
return a DataFrame with columns pos, p_value, gene, effect_size. The
fine-mapping loaders take ``filepath`` plus ``cs_col`` (output column name for
the credible set, default "cs") and return columns pos, pip, cs.
"""

from collections.abc import Callable
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any, Optional, Union

import pandas as pd

from .logging import logger
from .schemas import Family, Tier, spec
from .validation import ColumnSpec, check

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
        log_fmt: Debug log template; ``{n}`` is filled with the row count.
        read: Keyword arguments forwarded to ``pandas.read_csv``.
        schema: Validation to run after mapping. Takes the resolved output
            columns and returns the format's ``ColumnSpec``. None skips
            validation for a format that carries nothing checkable.
        schema_requires: Columns the format may structurally lack. An absent
            one is dropped from the schema's rules rather than failing, which
            is how CAVIAR and MatrixEQTL load without a position column.
        col_map: Static ``source -> target`` renames. Targets may be tokens.
        p_candidates: P-value source columns; the first present maps to ``p_col``.
        col_candidates: ``target -> candidates``; first present candidate maps to
            the (possibly token) target. Replaces hand-rolled first-match loops.
        transform: Optional ``(df, out_cols) -> df`` applied after renaming, for
            format-specific reshaping (e.g. credible-set assignment).
        gene_filter: ``"contains"`` or ``"exact"`` to enable eQTL gene filtering.
        clean_chrom: Strip a ``chr`` prefix from the ``chr`` column.
    """

    log_fmt: str
    read: dict[str, Any] = field(default_factory=dict)
    schema: Optional[Callable[[dict[str, str]], ColumnSpec]] = None
    schema_requires: tuple[str, ...] = ()
    col_map: dict[str, str] = field(default_factory=dict)
    p_candidates: tuple[str, ...] = ()
    col_candidates: dict[str, tuple[str, ...]] = field(default_factory=dict)
    transform: Optional[Callable[[pd.DataFrame, dict[str, str]], pd.DataFrame]] = None
    gene_filter: Optional[str] = None
    clean_chrom: bool = False


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


def _relax(spec: ColumnSpec, absent: tuple[str, ...]) -> ColumnSpec:
    """Drop every rule naming a column the file structurally lacks."""
    if not absent:
        return spec
    return replace(
        spec,
        required=tuple(c for c in spec.required if c not in absent),
        numeric=tuple(c for c in spec.numeric if c not in absent),
        not_null=tuple(c for c in spec.not_null if c not in absent),
        ranges=tuple(r for r in spec.ranges if r.column not in absent),
    )


def _validate(df: pd.DataFrame, spec: LoaderSpec, out_cols: dict[str, str]) -> None:
    """Validate ``df`` against the format's schema.

    Every format validates strictly, so a failed column mapping is reported
    at load time where the source columns are still known. A column named in
    ``schema_requires`` is exempt only when the file genuinely lacks it.
    """
    if spec.schema is None:
        return
    absent = tuple(c for c in spec.schema_requires if c not in df.columns)
    check(df, _relax(spec.schema(out_cols), absent))


def _load_tabular(
    filepath: Union[str, Path],
    spec: LoaderSpec,
    *,
    gene: Optional[str] = None,
    **out_cols: str,
) -> pd.DataFrame:
    """Load a tabular file per ``spec``: read, map, rename, transform, validate."""
    df = pd.read_csv(filepath, **spec.read)

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

    if spec.clean_chrom and "chr" in df.columns:
        df["chr"] = df["chr"].astype(str).str.replace("chr", "", regex=False)

    logger.debug(spec.log_fmt.format(n=len(df)))
    _validate(df, spec, out_cols)
    return df


def _gwas_schema(out_cols: dict[str, str]) -> ColumnSpec:
    """Strict GWAS contract over the caller's position/p-value columns."""
    return spec(
        Family.GWAS, Tier.LOAD, pos_col=out_cols["pos_col"], p_col=out_cols["p_col"]
    )


# =============================================================================
# GWAS Loaders
# =============================================================================


# No comment= here on purpose. PLINK 2's --glm header line starts with "#CHROM",
# so a "#" comment prefix makes pandas discard the header and promote the first
# variant to column labels. PLINK writes no comment lines of its own.
_PLINK_SPEC = LoaderSpec(
    log_fmt="Loaded PLINK file with {n} variants",
    read={"sep": r"\s+"},
    col_candidates={
        "pos_col": ("BP", "POS", "bp", "pos"),
        "p_col": ("P", "P_BOLT_LMM", "p", "PVAL", "pval", "P_LINREG"),
        "rs_col": ("SNP", "ID", "rsid", "RSID", "MarkerName", "variant_id"),
        "chr": ("CHR", "chr", "CHROM", "chrom", "#CHROM"),
    },
    schema=_gwas_schema,
)


def load_plink_assoc(
    filepath: Union[str, Path],
    pos_col: str = "ps",
    p_col: str = "p_wald",
    rs_col: str = "rs",
) -> pd.DataFrame:
    """Load PLINK association results (.assoc, .assoc.linear, .assoc.logistic, .qassoc).

    Automatically detects PLINK format variant and maps columns to standard names.
    See the module docstring for the shared GWAS loader arguments and return value.

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
    log_fmt="Loaded REGENIE file with {n} variants",
    read={"sep": r"\s+", "comment": "#"},
    col_map={"GENPOS": "pos_col", "ID": "rs_col", "CHROM": "chr"},
    transform=_regenie_pvalue,
    schema=_gwas_schema,
)


def load_regenie(
    filepath: Union[str, Path],
    pos_col: str = "ps",
    p_col: str = "p_wald",
    rs_col: str = "rs",
) -> pd.DataFrame:
    """Load REGENIE association results (.regenie).

    Prefers the computed LOG10P column over a raw P column.
    See the module docstring for the shared GWAS loader arguments and return value.

    Example:
        >>> gwas_df = load_regenie("results.regenie")
    """
    return _load_tabular(
        filepath, _REGENIE_SPEC, pos_col=pos_col, p_col=p_col, rs_col=rs_col
    )


_BOLT_LMM_SPEC = LoaderSpec(
    log_fmt="Loaded BOLT-LMM file with {n} variants",
    read={"sep": "\t"},
    col_map={"BP": "pos_col", "SNP": "rs_col", "CHR": "chr"},
    p_candidates=("P_BOLT_LMM", "P_BOLT_LMM_INF"),
    schema=_gwas_schema,
)


def load_bolt_lmm(
    filepath: Union[str, Path],
    pos_col: str = "ps",
    p_col: str = "p_wald",
    rs_col: str = "rs",
) -> pd.DataFrame:
    """Load BOLT-LMM association results (.stats).

    See the module docstring for the shared GWAS loader arguments and return value.

    Example:
        >>> gwas_df = load_bolt_lmm("results.stats")
    """
    return _load_tabular(
        filepath, _BOLT_LMM_SPEC, pos_col=pos_col, p_col=p_col, rs_col=rs_col
    )


_GEMMA_SPEC = LoaderSpec(
    log_fmt="Loaded GEMMA file with {n} variants",
    read={"sep": "\t"},
    col_map={"ps": "pos_col", "rs": "rs_col", "chr": "chr"},
    p_candidates=("p_wald", "p_lrt", "p_score"),
    schema=_gwas_schema,
)


def load_gemma(
    filepath: Union[str, Path],
    pos_col: str = "ps",
    p_col: str = "p_wald",
    rs_col: str = "rs",
) -> pd.DataFrame:
    """Load GEMMA association results (.assoc.txt).

    See the module docstring for the shared GWAS loader arguments and return value.

    Example:
        >>> gwas_df = load_gemma("output.assoc.txt")
    """
    return _load_tabular(
        filepath, _GEMMA_SPEC, pos_col=pos_col, p_col=p_col, rs_col=rs_col
    )


_SAIGE_SPEC = LoaderSpec(
    log_fmt="Loaded SAIGE file with {n} variants",
    read={"sep": "\t"},
    col_map={"POS": "pos_col", "MarkerID": "rs_col", "CHR": "chr"},
    p_candidates=("p.value.NA", "p.value"),
    schema=_gwas_schema,
)


def load_saige(
    filepath: Union[str, Path],
    pos_col: str = "ps",
    p_col: str = "p_wald",
    rs_col: str = "rs",
) -> pd.DataFrame:
    """Load SAIGE association results.

    See the module docstring for the shared GWAS loader arguments and return value.

    Example:
        >>> gwas_df = load_saige("results.txt")
    """
    return _load_tabular(
        filepath, _SAIGE_SPEC, pos_col=pos_col, p_col=p_col, rs_col=rs_col
    )


_GWAS_CATALOG_SPEC = LoaderSpec(
    log_fmt="Loaded GWAS Catalog file with {n} variants",
    read={"sep": "\t"},
    col_map={
        "base_pair_location": "pos_col",
        "variant_id": "rs_col",
        "chromosome": "chr",
        "p_value": "p_col",
    },
    schema=_gwas_schema,
)


def load_gwas_catalog(
    filepath: Union[str, Path],
    pos_col: str = "ps",
    p_col: str = "p_wald",
    rs_col: str = "rs",
) -> pd.DataFrame:
    """Load GWAS Catalog summary statistics format.

    See the module docstring for the shared GWAS loader arguments and return value.
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
    log_fmt="Loaded SuSiE file with {n} variants",
    read={"sep": "\t"},
    col_candidates={
        "pos": ("pos", "position", "BP", "bp", "POS"),
        "pip": ("pip", "PIP", "posterior_prob", "prob"),
        "cs_col": ("cs", "CS", "credible_set", "cs_index", "L"),
        "rs": ("snp", "SNP", "variant_id", "rsid"),
    },
    transform=_susie_cs,
    schema=lambda out_cols: spec(Family.FINEMAPPING, Tier.LOAD),
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
    log_fmt="Loaded FINEMAP file with {n} variants",
    read={"sep": r"\s+"},
    col_map={"position": "pos", "prob": "pip", "rsid": "rs", "chromosome": "chr"},
    transform=_finemap_cs,
    schema=lambda out_cols: spec(Family.FINEMAPPING, Tier.LOAD),
)


def load_finemap(
    filepath: Union[str, Path],
    cs_col: str = "cs",
) -> pd.DataFrame:
    """Load FINEMAP results (.snp output file).

    FINEMAP reports no credible set, so a 95% set is assigned from the
    cumulative PIP. See the module docstring for the shared fine-mapping
    arguments and return value.

    Example:
        >>> fm_df = load_finemap("results.snp")
    """
    return _load_tabular(filepath, _FINEMAP_SPEC, cs_col=cs_col)


_CAVIAR_SPEC = LoaderSpec(
    log_fmt="Loaded CAVIAR file with {n} variants",
    # CAVIAR .set files are headerless: SNP_ID Causal_Post_Prob.
    read={"sep": r"\s+", "header": None, "names": ["rs", "pip"]},
    transform=_finemap_cs,
    schema=lambda out_cols: spec(Family.FINEMAPPING, Tier.LOAD),
    # CAVIAR carries no position; callers merge a SNP annotation to add it.
    schema_requires=("pos",),
)


def load_caviar(
    filepath: Union[str, Path],
    cs_col: str = "cs",
) -> pd.DataFrame:
    """Load CAVIAR results (.set output file).

    A 95% credible set is assigned from the cumulative PIP. CAVIAR carries no
    position, so merge a SNP annotation file to add ``pos`` before plotting.
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
    schema=lambda out_cols: spec(Family.FINEMAPPING, Tier.LOAD),
)


def load_polyfun(
    filepath: Union[str, Path],
    cs_col: str = "cs",
) -> pd.DataFrame:
    """Load PolyFun/SuSiE fine-mapping results.

    See the module docstring for the shared fine-mapping arguments and return value.
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
    check(result, spec(Family.GENES, Tier.LOAD))
    return result


# Standard BED column names, in order, up to BED12.
_BED_COLUMNS = (
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
)


def _bed_names(df: pd.DataFrame, out_cols: dict[str, str]) -> pd.DataFrame:
    """Name a headerless BED frame, generically past BED12."""
    if not all(isinstance(c, int) for c in df.columns):
        return df
    n_cols = len(df.columns)
    names = list(_BED_COLUMNS[:n_cols])
    names += [f"col{i}" for i in range(len(_BED_COLUMNS), n_cols)]
    df.columns = names
    return df


_BED_SPEC = LoaderSpec(
    log_fmt="Loaded {n} features from BED",
    read={"sep": "\t"},
    col_map={
        "chrom": "chr",
        "chromStart": "start",
        "chromEnd": "end",
        "name": "gene_name",
    },
    transform=_bed_names,
    clean_chrom=True,
    schema=lambda out_cols: spec(Family.GENES, Tier.LOAD),
)


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
    spec = replace(
        _BED_SPEC, read={**_BED_SPEC.read, "header": 0 if has_header else None}
    )
    return _load_tabular(filepath, spec)


def _ensembl_strand(df: pd.DataFrame, out_cols: dict[str, str]) -> pd.DataFrame:
    """Convert Ensembl 1/-1 strand encoding to +/- symbols."""
    if "strand" in df.columns:
        df["strand"] = df["strand"].map({1: "+", -1: "-", "+": "+", "-": "-"})
    return df


_ENSEMBL_SPEC = LoaderSpec(
    log_fmt="Loaded {n} genes from Ensembl export",
    read={"sep": "\t"},
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
    schema=lambda out_cols: spec(Family.GENES, Tier.LOAD),
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


# Filename substrings that identify a GWAS format.
_FORMAT_HINTS: tuple[tuple[str, str], ...] = (
    (".assoc", "plink"),
    (".qassoc", "plink"),
    (".glm.", "plink"),
    (".regenie", "regenie"),
    (".stats", "bolt"),
    ("gemma", "gemma"),
    (".assoc.txt", "gemma"),
    ("saige", "saige"),
)

_GWAS_LOADERS: dict[str, Callable[..., pd.DataFrame]] = {
    "plink": load_plink_assoc,
    "regenie": load_regenie,
    "bolt": load_bolt_lmm,
    "gemma": load_gemma,
    "saige": load_saige,
    "catalog": load_gwas_catalog,
}


def _detect_format(filepath: Path) -> str:
    """Identify a GWAS format from the filename, defaulting to PLINK."""
    name = filepath.name.lower()
    for hint, format in sorted(_FORMAT_HINTS, key=lambda h: len(h[0]), reverse=True):
        if hint in name:
            return format
    logger.warning(
        f"Could not detect GWAS file format for '{filepath}'. "
        "Defaulting to 'plink'. Specify format= parameter explicitly."
    )
    return "plink"


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

    Raises:
        ValueError: If ``format`` names an unknown format.

    Example:
        >>> # Auto-detect format
        >>> gwas_df = load_gwas("results.assoc.linear")
        >>>
        >>> # Explicit format
        >>> gwas_df = load_gwas("results.txt", format="regenie")
    """
    filepath = Path(filepath)
    format = format or _detect_format(filepath)

    if format not in _GWAS_LOADERS:
        raise ValueError(
            f"Unknown format '{format}'. Options: {list(_GWAS_LOADERS.keys())}"
        )

    return _GWAS_LOADERS[format](
        filepath, pos_col=pos_col, p_col=p_col, rs_col=rs_col, **kwargs
    )
