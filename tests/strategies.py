"""Hypothesis strategies for GWAS-domain test data generation.

Provides reusable strategies for generating valid DataFrames that match
pyLocusZoom's expected input formats.
"""

import pandas as pd
from hypothesis import strategies as st

# =============================================================================
# Primitive Strategies
# =============================================================================


@st.composite
def chromosomes(draw, species: str = "canine") -> int:
    """Generate valid chromosome numbers.

    Args:
        species: Species name for chromosome count.
            - canine: 1-38
            - human: 1-22
            - feline: 1-18

    Returns:
        Integer chromosome number.
    """
    max_chrom = {"canine": 38, "human": 22, "feline": 18}.get(species, 38)
    return draw(st.integers(min_value=1, max_value=max_chrom))


@st.composite
def positions(draw, min_val: int = 1, max_val: int = 300_000_000) -> int:
    """Generate valid genomic positions.

    Args:
        min_val: Minimum position (inclusive).
        max_val: Maximum position (inclusive).

    Returns:
        Integer genomic position.
    """
    return draw(st.integers(min_value=min_val, max_value=max_val))


@st.composite
def pvalues(draw, allow_zero: bool = False, allow_one: bool = True) -> float:
    """Generate valid p-values.

    Args:
        allow_zero: If True, can return 0.0.
        allow_one: If True, can return 1.0.

    Returns:
        Float p-value in (0, 1] or [0, 1] depending on flags.
    """
    min_val = 0.0 if allow_zero else 1e-300
    max_val = 1.0 if allow_one else 1.0 - 1e-10
    return draw(st.floats(min_value=min_val, max_value=max_val, allow_nan=False))


@st.composite
def ld_values(draw, allow_nan: bool = True) -> float:
    """Generate valid LD R² values.

    Args:
        allow_nan: If True, can return NaN for missing LD.

    Returns:
        Float R² value in [0, 1] or NaN.
    """
    if allow_nan and draw(st.booleans()):
        return float("nan")
    return draw(st.floats(min_value=0.0, max_value=1.0, allow_nan=False))


# =============================================================================
# DataFrame Strategies
# =============================================================================


@st.composite
def gwas_dataframes(
    draw,
    min_snps: int = 1,
    max_snps: int = 100,
    chrom: int | None = None,
    start: int | None = None,
    end: int | None = None,
) -> pd.DataFrame:
    """Generate valid GWAS DataFrames.

    Args:
        min_snps: Minimum number of SNPs.
        max_snps: Maximum number of SNPs.
        chrom: Fixed chromosome (random if None).
        start: Region start (random if None).
        end: Region end (random if None).

    Returns:
        DataFrame with rs, chr, ps, p_wald columns.
    """
    n_snps = draw(st.integers(min_value=min_snps, max_value=max_snps))

    # Determine chromosome
    if chrom is None:
        chrom = draw(chromosomes())

    # Determine region
    if start is None or end is None:
        region_start = draw(st.integers(min_value=1_000_000, max_value=100_000_000))
        region_size = draw(st.integers(min_value=100_000, max_value=10_000_000))
        start = region_start
        end = region_start + region_size

    # Generate positions within region
    pos_list = sorted(
        draw(
            st.lists(
                st.integers(min_value=start, max_value=end),
                min_size=n_snps,
                max_size=n_snps,
                unique=True,
            )
        )
    )

    # Generate p-values (biased toward small values for realism)
    p_vals = [draw(pvalues()) for _ in range(n_snps)]

    return pd.DataFrame(
        {
            "rs": [f"rs{i}" for i in range(n_snps)],
            "chr": [chrom] * n_snps,
            "pos": pos_list,
            "p_value": p_vals,
        }
    )


@st.composite
def gwas_dataframes_multichrom(
    draw,
    min_snps_per_chrom: int = 5,
    max_snps_per_chrom: int = 50,
    species: str = "canine",
) -> pd.DataFrame:
    """Generate GWAS DataFrames spanning multiple chromosomes.

    Args:
        min_snps_per_chrom: Minimum SNPs per chromosome.
        max_snps_per_chrom: Maximum SNPs per chromosome.
        species: Species for chromosome count.

    Returns:
        DataFrame with SNPs across multiple chromosomes.
    """
    max_chrom = {"canine": 38, "human": 22, "feline": 18}.get(species, 38)
    n_chroms = draw(st.integers(min_value=2, max_value=min(5, max_chrom)))
    selected_chroms = draw(
        st.lists(
            st.integers(min_value=1, max_value=max_chrom),
            min_size=n_chroms,
            max_size=n_chroms,
            unique=True,
        )
    )

    dfs = []
    for chrom in selected_chroms:
        chrom_df = draw(
            gwas_dataframes(
                min_snps=min_snps_per_chrom,
                max_snps=max_snps_per_chrom,
                chrom=chrom,
            )
        )
        dfs.append(chrom_df)

    result = pd.concat(dfs, ignore_index=True)
    # Make rs IDs unique across chromosomes
    result["rs"] = [f"rs{i}" for i in range(len(result))]
    return result


@st.composite
def gene_dataframes(
    draw,
    min_genes: int = 0,
    max_genes: int = 20,
    chrom: int | None = None,
    start: int | None = None,
    end: int | None = None,
) -> pd.DataFrame:
    """Generate valid gene annotation DataFrames.

    Args:
        min_genes: Minimum number of genes.
        max_genes: Maximum number of genes.
        chrom: Fixed chromosome (random if None).
        start: Region start (random if None).
        end: Region end (random if None).

    Returns:
        DataFrame with gene_name, chr, start, end, strand columns.
    """
    n_genes = draw(st.integers(min_value=min_genes, max_value=max_genes))

    if n_genes == 0:
        return pd.DataFrame(columns=["gene_name", "chr", "start", "end", "strand"])

    if chrom is None:
        chrom = draw(chromosomes())

    if start is None or end is None:
        region_start = draw(st.integers(min_value=1_000_000, max_value=100_000_000))
        region_size = draw(st.integers(min_value=500_000, max_value=10_000_000))
        start = region_start
        end = region_start + region_size

    genes = []
    for i in range(n_genes):
        gene_start = draw(st.integers(min_value=start, max_value=end - 1000))
        gene_size = draw(
            st.integers(min_value=1000, max_value=min(100_000, end - gene_start))
        )
        gene_end = gene_start + gene_size
        strand = draw(st.sampled_from(["+", "-"]))
        genes.append(
            {
                "gene_name": f"GENE_{i}",
                "chr": chrom,
                "start": gene_start,
                "end": gene_end,
                "strand": strand,
            }
        )

    return pd.DataFrame(genes)


@st.composite
def pvalues_invalid(draw) -> float:
    """Generate invalid p-values for error handling tests.

    Returns:
        Float that is invalid as a p-value (negative, >1, or NaN).
    """
    return draw(
        st.one_of(
            st.floats(max_value=-0.001, allow_nan=False),  # Negative
            st.floats(min_value=1.001, allow_nan=False),  # > 1
            st.just(float("nan")),  # NaN
        )
    )
