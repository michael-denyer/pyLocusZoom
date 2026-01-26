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

from pathlib import Path
from typing import Optional, Union

import pandas as pd

from .logging import logger

# =============================================================================
# GWAS Loaders
# =============================================================================


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
    df = pd.read_csv(filepath, sep=r"\s+", comment="#")

    # Standardize column names (PLINK uses various conventions)
    col_map = {}

    # Position columns
    for col in ["BP", "POS", "bp", "pos"]:
        if col in df.columns:
            col_map[col] = pos_col
            break

    # P-value columns
    for col in ["P", "P_BOLT_LMM", "p", "PVAL", "pval", "P_LINREG"]:
        if col in df.columns:
            col_map[col] = p_col
            break

    # SNP ID columns
    for col in ["SNP", "ID", "rsid", "RSID", "MarkerName", "variant_id"]:
        if col in df.columns:
            col_map[col] = rs_col
            break

    # Chromosome column (keep as "chr" for reference)
    for col in ["CHR", "chr", "CHROM", "chrom", "#CHROM"]:
        if col in df.columns:
            col_map[col] = "chr"
            break

    df = df.rename(columns=col_map)
    logger.debug(f"Loaded PLINK file with {len(df)} variants")
    return df


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
    df = pd.read_csv(filepath, sep=r"\s+", comment="#")

    col_map = {
        "GENPOS": pos_col,
        "ID": rs_col,
        "CHROM": "chr",
    }

    # REGENIE uses LOG10P, need to convert
    if "LOG10P" in df.columns:
        df[p_col] = 10 ** (-df["LOG10P"])
    elif "P" in df.columns:
        col_map["P"] = p_col

    df = df.rename(columns=col_map)
    logger.debug(f"Loaded REGENIE file with {len(df)} variants")
    return df


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
    df = pd.read_csv(filepath, sep="\t")

    col_map = {
        "BP": pos_col,
        "SNP": rs_col,
        "CHR": "chr",
        "P_BOLT_LMM_INF": p_col,  # Infinitesimal model (default)
    }

    # Prefer P_BOLT_LMM if available (full model)
    if "P_BOLT_LMM" in df.columns:
        col_map["P_BOLT_LMM"] = p_col
        if "P_BOLT_LMM_INF" in col_map:
            del col_map["P_BOLT_LMM_INF"]

    df = df.rename(columns=col_map)
    logger.debug(f"Loaded BOLT-LMM file with {len(df)} variants")
    return df


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
    df = pd.read_csv(filepath, sep="\t")

    col_map = {
        "ps": pos_col,
        "rs": rs_col,
        "chr": "chr",
        "p_wald": p_col,
        "p_lrt": p_col,  # Alternative if p_wald not present
        "p_score": p_col,  # Alternative
    }

    # Only map first matching p-value column
    p_cols = ["p_wald", "p_lrt", "p_score"]
    found_p = False
    for p in p_cols:
        if p in df.columns and not found_p:
            col_map[p] = p_col
            found_p = True
        elif p in col_map:
            del col_map[p]

    df = df.rename(columns=col_map)
    logger.debug(f"Loaded GEMMA file with {len(df)} variants")
    return df


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
    df = pd.read_csv(filepath, sep="\t")

    col_map = {
        "POS": pos_col,
        "MarkerID": rs_col,
        "CHR": "chr",
        "p.value": p_col,
        "p.value.NA": p_col,  # SPA-adjusted
    }

    df = df.rename(columns=col_map)
    logger.debug(f"Loaded SAIGE file with {len(df)} variants")
    return df


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
    df = pd.read_csv(filepath, sep="\t")

    col_map = {
        "base_pair_location": pos_col,
        "variant_id": rs_col,
        "chromosome": "chr",
        "p_value": p_col,
    }

    df = df.rename(columns=col_map)
    logger.debug(f"Loaded GWAS Catalog file with {len(df)} variants")
    return df


# =============================================================================
# eQTL Loaders
# =============================================================================


def load_gtex_eqtl(
    filepath: Union[str, Path],
    gene: Optional[str] = None,
) -> pd.DataFrame:
    """Load GTEx eQTL significant pairs format.

    Args:
        filepath: Path to GTEx eQTL file (e.g., signif_variant_gene_pairs.txt.gz).
        gene: Optional gene to filter to (ENSG ID or gene symbol).

    Returns:
        DataFrame with columns: pos, p_value, gene, effect.

    Example:
        >>> eqtl_df = load_gtex_eqtl("GTEx_Analysis.signif_pairs.txt.gz", gene="BRCA1")
    """
    # GTEx files are often gzipped
    df = pd.read_csv(filepath, sep="\t")

    # Map GTEx columns to standard format
    col_map = {}

    # Variant position (GTEx uses variant_id like chr1_12345_A_G_b38)
    if "variant_id" in df.columns:
        # Extract position from variant_id
        df["pos"] = df["variant_id"].str.split("_").str[1].astype(int)
    elif "pos" not in df.columns:
        for col in ["tss_distance", "POS"]:
            if col in df.columns:
                col_map[col] = "pos"
                break

    # P-value
    for col in ["pval_nominal", "p_value", "pvalue", "P"]:
        if col in df.columns:
            col_map[col] = "p_value"
            break

    # Gene
    for col in ["gene_id", "gene_name", "phenotype_id"]:
        if col in df.columns:
            col_map[col] = "gene"
            break

    # Effect size (slope)
    for col in ["slope", "beta", "effect_size"]:
        if col in df.columns:
            col_map[col] = "effect"
            break

    df = df.rename(columns=col_map)

    # Filter to gene if specified
    if gene is not None and "gene" in df.columns:
        # Match either ENSG ID or gene symbol
        mask = df["gene"].str.contains(gene, case=False, na=False)
        df = df[mask]

    logger.debug(f"Loaded GTEx eQTL file with {len(df)} associations")
    return df


def load_eqtl_catalogue(
    filepath: Union[str, Path],
    gene: Optional[str] = None,
) -> pd.DataFrame:
    """Load eQTL Catalogue format.

    Args:
        filepath: Path to eQTL Catalogue file.
        gene: Optional gene to filter to.

    Returns:
        DataFrame with columns: pos, p_value, gene, effect.
    """
    df = pd.read_csv(filepath, sep="\t")

    col_map = {
        "position": "pos",
        "pvalue": "p_value",
        "gene_id": "gene",
        "beta": "effect",
        "chromosome": "chr",
    }

    df = df.rename(columns=col_map)

    if gene is not None and "gene" in df.columns:
        mask = df["gene"].str.contains(gene, case=False, na=False)
        df = df[mask]

    logger.debug(f"Loaded eQTL Catalogue file with {len(df)} associations")
    return df


def load_matrixeqtl(
    filepath: Union[str, Path],
    gene: Optional[str] = None,
) -> pd.DataFrame:
    """Load MatrixEQTL output format.

    Args:
        filepath: Path to MatrixEQTL output file.
        gene: Optional gene to filter to.

    Returns:
        DataFrame with columns: pos, p_value, gene, effect.

    Note:
        MatrixEQTL output doesn't include position by default.
        You may need to merge with a SNP annotation file.
    """
    df = pd.read_csv(filepath, sep="\t")

    col_map = {
        "SNP": "rs",
        "gene": "gene",
        "p-value": "p_value",
        "pvalue": "p_value",
        "beta": "effect",
        "t-stat": "t_stat",
    }

    df = df.rename(columns=col_map)

    if gene is not None and "gene" in df.columns:
        df = df[df["gene"] == gene]

    logger.debug(f"Loaded MatrixEQTL file with {len(df)} associations")
    return df


# =============================================================================
# Fine-mapping Loaders
# =============================================================================


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
    df = pd.read_csv(filepath, sep="\t")

    col_map = {}

    # Position
    for col in ["pos", "position", "BP", "bp", "POS"]:
        if col in df.columns:
            col_map[col] = "pos"
            break

    # PIP (posterior inclusion probability)
    for col in ["pip", "PIP", "posterior_prob", "prob"]:
        if col in df.columns:
            col_map[col] = "pip"
            break

    # Credible set
    for col in ["cs", "CS", "credible_set", "cs_index", "L"]:
        if col in df.columns:
            col_map[col] = cs_col
            break

    # SNP ID
    for col in ["snp", "SNP", "variant_id", "rsid"]:
        if col in df.columns:
            col_map[col] = "rs"
            break

    df = df.rename(columns=col_map)

    # SuSiE uses -1 or NA for variants not in a credible set; standardize to 0
    if cs_col in df.columns:
        df[cs_col] = df[cs_col].fillna(0).astype(int)
        df.loc[df[cs_col] < 0, cs_col] = 0

    logger.debug(f"Loaded SuSiE file with {len(df)} variants")
    return df


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
    df = pd.read_csv(filepath, sep=r"\s+")

    col_map = {
        "position": "pos",
        "prob": "pip",
        "rsid": "rs",
        "chromosome": "chr",
    }

    df = df.rename(columns=col_map)

    # FINEMAP doesn't directly output credible sets
    # Assign based on cumulative PIP threshold (95% default)
    if cs_col not in df.columns and "pip" in df.columns:
        df = df.sort_values("pip", ascending=False)
        df["cumsum_pip"] = df["pip"].cumsum()
        df[cs_col] = (df["cumsum_pip"] <= 0.95).astype(int)
        df = df.drop(columns=["cumsum_pip"])

    logger.debug(f"Loaded FINEMAP file with {len(df)} variants")
    return df


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
    return df


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
    df = pd.read_csv(filepath, sep=r"\s+")

    col_map = {
        "BP": "pos",
        "PIP": "pip",
        "SNP": "rs",
        "CHR": "chr",
        "CREDIBLE_SET": cs_col,
    }

    df = df.rename(columns=col_map)

    if cs_col in df.columns:
        df[cs_col] = df[cs_col].fillna(0).astype(int)

    logger.debug(f"Loaded PolyFun file with {len(df)} variants")
    return df


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
        col_names = ["chr", "start", "end"]
        if n_cols >= 4:
            col_names.append("gene_name")
        if n_cols >= 5:
            col_names.append("score")
        if n_cols >= 6:
            col_names.append("strand")
        df.columns = col_names[:n_cols]

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
    return df


def load_ensembl_genes(
    filepath: Union[str, Path],
) -> pd.DataFrame:
    """Load Ensembl BioMart gene export.

    Args:
        filepath: Path to BioMart export file (TSV).

    Returns:
        DataFrame with columns: chr, start, end, gene_name, strand.
    """
    df = pd.read_csv(filepath, sep="\t")

    col_map = {
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
    }

    df = df.rename(columns=col_map)

    # Convert strand (Ensembl uses 1/-1)
    if "strand" in df.columns:
        df["strand"] = df["strand"].map({1: "+", -1: "-", "+": "+", "-": "-"})

    logger.debug(f"Loaded {len(df)} genes from Ensembl export")
    return df


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
            format = "plink"  # Default fallback

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

    return loaders[format](filepath, pos_col=pos_col, p_col=p_col, rs_col=rs_col)
