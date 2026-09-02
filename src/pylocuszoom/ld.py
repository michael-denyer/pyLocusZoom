"""LD (Linkage Disequilibrium) calculation using PLINK.

Calculates R² values between a lead SNP and all other SNPs in a region
using PLINK 1.9's --r2 command.
"""

import os
import re
import shutil
import subprocess
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Literal, Optional, Union

import pandas as pd

from .exceptions import EmptyLDOutputError, PlinkError, ValidationError
from .logging import logger


def _log_rmtree_error(func, path, exc_info):
    """Log shutil.rmtree cleanup failures instead of silently ignoring them."""
    logger.debug(f"Failed to clean up temp file {path}: {exc_info[1]}")


LDMetric = Literal["r2", "dprime"]

_METRIC_FLAGS: dict[str, tuple[str, ...]] = {
    "r2": ("--r2", "square"),
    "dprime": ("--r", "dprime", "square"),
}


def _metric_flags(metric: str) -> tuple[str, ...]:
    """Return the PLINK flags for an LD metric, rejecting any other spelling."""
    try:
        return _METRIC_FLAGS[metric]
    except KeyError:
        raise ValidationError(
            f"Unknown LD metric {metric!r}; expected one of {sorted(_METRIC_FLAGS)}"
        ) from None


def _add_species_flags(cmd: list[str], species: str | None) -> None:
    """Add species-specific flags to PLINK command.

    Args:
        cmd: Command list to append flags to (modified in-place).
        species: Species name ("canine", "feline", or None for no species-specific flags).
    """
    if species == "canine":
        cmd.append("--dog")
    elif species == "feline":
        cmd.extend(["--chr-set", "18"])


def build_pairwise_ld_command(
    plink_path: str,
    bfile_path: str,
    output_path: str,
    snp_list_file: Optional[str] = None,
    chrom: Optional[int] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
    species: Optional[str] = "canine",
    metric: LDMetric = "r2",
) -> list:
    """Build PLINK command for pairwise LD matrix computation.

    Generates command for computing an N x N LD matrix using PLINK's
    --r2 square (or --r dprime square) command.

    Args:
        plink_path: Path to PLINK executable.
        bfile_path: Input binary fileset prefix (.bed/.bim/.fam).
        output_path: Output prefix (creates .ld and .snplist files).
        snp_list_file: Path to file with SNP IDs to extract (one per line).
        chrom: Chromosome number for region-based extraction.
        start: Start position (bp) for region-based extraction.
        end: End position (bp) for region-based extraction.
        species: Species flag ('canine', 'feline', or None for human).
        metric: LD metric ('r2' or 'dprime').

    Returns:
        List of command arguments for subprocess.
    """
    cmd = [plink_path]
    _add_species_flags(cmd, species)

    # Input and output
    cmd.extend(["--bfile", bfile_path])
    cmd.extend(["--out", output_path])

    cmd.extend(_metric_flags(metric))

    # Track SNP order in output
    cmd.append("--write-snplist")

    # SNP extraction mode
    if snp_list_file:
        cmd.extend(["--extract", snp_list_file])

    # Region-based extraction
    if chrom is not None:
        cmd.extend(["--chr", str(chrom)])
    if start is not None:
        cmd.extend(["--from-bp", str(start)])
    if end is not None:
        cmd.extend(["--to-bp", str(end)])

    return cmd


def find_plink() -> Optional[str]:
    """Find PLINK executable on PATH.

    Checks for plink1.9 first, then plink.

    Returns:
        Path to PLINK executable, or None if not found.
    """
    for name in ["plink1.9", "plink"]:
        path = shutil.which(name)
        if path:
            return path
    return None


def validate_plink_files(bfile_path: Union[str, Path]) -> Path:
    """Validate that PLINK binary fileset exists.

    Checks for .bed, .bim, and .fam files.

    Args:
        bfile_path: Path prefix for PLINK files (without extension).

    Returns:
        Path object if files exist.

    Raises:
        ValidationError: If any PLINK files are missing.
    """
    path = Path(bfile_path)
    missing = []
    # Use string concatenation rather than with_suffix() — PLINK prefixes
    # frequently contain dots (e.g. "ukbb.v3"), which with_suffix would
    # truncate, causing validation to check the wrong file.
    for ext in [".bed", ".bim", ".fam"]:
        if not Path(str(path) + ext).exists():
            missing.append(ext)

    if missing:
        raise ValidationError(
            f"PLINK files missing for {path}: {missing}. "
            f"Expected: {path}.bed, {path}.bim, {path}.fam"
        )
    return path


def _resolve_plink(plink_path: Optional[str]) -> str:
    """Return the PLINK executable to run.

    Args:
        plink_path: Caller-supplied path, or None to search PATH.

    Returns:
        Path to the PLINK executable.

    Raises:
        FileNotFoundError: If PLINK is neither supplied nor on PATH.
    """
    resolved = plink_path or find_plink()
    if resolved is None:
        raise FileNotFoundError(
            "PLINK not found. Install PLINK 1.9 or specify plink_path."
        )
    logger.debug(f"Using PLINK at {resolved}")
    return resolved


@contextmanager
def _plink_workdir(working_dir: Optional[str], prefix: str) -> Iterator[str]:
    """Yield a directory for PLINK output, removing it if we created it.

    Args:
        working_dir: Caller-supplied directory, or None for a temp directory.
        prefix: Temp-directory name prefix, used only when one is created.

    Yields:
        Path to a directory that exists for the duration of the block.
    """
    created = working_dir is None
    if working_dir is None:
        working_dir = tempfile.mkdtemp(prefix=prefix)
    try:
        os.makedirs(working_dir, exist_ok=True)
        yield working_dir
    finally:
        if created and os.path.exists(working_dir):
            shutil.rmtree(working_dir, onerror=_log_rmtree_error)


def _run_plink(cmd: list[str], working_dir: str, what: str, hint: str) -> None:
    """Run a PLINK command, raising PlinkError on a timeout or a non-zero exit.

    Args:
        cmd: Command arguments, as built by one of the build_*_command helpers.
        working_dir: Directory to run in.
        what: Name of the calculation, used in error messages.
        hint: What the caller can try after a timeout.

    Raises:
        PlinkError: If PLINK times out or exits non-zero.
    """
    logger.debug(f"Running PLINK command: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd,
            cwd=working_dir,
            capture_output=True,
            text=True,
            timeout=300,
        )
    except subprocess.TimeoutExpired:
        raise PlinkError(f"PLINK {what} timed out after 300s. {hint}")

    if result.returncode != 0:
        raise PlinkError(
            f"PLINK {what} failed (exit code {result.returncode}): "
            f"{result.stderr.strip()}"
        )


def build_ld_command(
    plink_path: str,
    bfile_path: str,
    lead_snp: str,
    output_path: str,
    window_kb: int = 500,
    ld_window_r2: float = 0.0,
    species: Optional[str] = "canine",
    threads: Optional[int] = None,
) -> list:
    """Build PLINK command for LD calculation.

    Args:
        plink_path: Path to PLINK executable.
        bfile_path: Input binary fileset prefix (.bed/.bim/.fam).
        lead_snp: SNP ID to calculate LD against.
        output_path: Output prefix (creates .ld file).
        window_kb: Window size in kilobases.
        ld_window_r2: Minimum R² to report (0.0 reports all).
        species: Species flag for PLINK ('canine', 'feline', or None for human).
        threads: Number of threads (auto-detect if None).

    Returns:
        List of command arguments for subprocess.
    """
    cmd = [plink_path]
    _add_species_flags(cmd, species)

    # Input and output
    cmd.extend(["--bfile", bfile_path])
    cmd.extend(["--out", output_path])

    # LD calculation flags
    cmd.append("--r2")
    cmd.extend(["--ld-snp", lead_snp])
    cmd.extend(["--ld-window-kb", str(window_kb)])
    cmd.extend(["--ld-window", "99999"])  # Remove default 10 SNP limit
    cmd.extend(["--ld-window-r2", str(ld_window_r2)])

    # Threads
    if threads is None:
        threads = os.cpu_count() or 1
    cmd.extend(["--threads", str(threads)])

    return cmd


def parse_pairwise_ld_output(
    ld_file: str, snplist_file: str
) -> tuple[pd.DataFrame, list[str]]:
    """Parse PLINK pairwise LD matrix output files.

    PLINK --r2 square outputs:
    - .ld file: N x N matrix of R2/D' values (whitespace-separated, no headers)
    - .snplist file: SNP IDs in order (one per line)

    Args:
        ld_file: Path to .ld output file (square matrix).
        snplist_file: Path to .snplist output file (SNP IDs).

    Returns:
        Tuple of (DataFrame with R2/D' values, list of SNP IDs).
        DataFrame has SNP IDs as both index and columns.

    Raises:
        PlinkError: If output files are missing after a successful PLINK run.
    """
    if not os.path.exists(ld_file) or not os.path.exists(snplist_file):
        missing = [f for f in (ld_file, snplist_file) if not os.path.exists(f)]
        raise PlinkError(
            f"PLINK reported success but output files missing: {missing}. "
            f"This may indicate a filesystem issue or PLINK bug."
        )

    with open(snplist_file) as f:
        snp_ids = [line.strip() for line in f if line.strip()]

    if not snp_ids:
        raise PlinkError(
            f"PLINK reported success but snplist file is empty: {snplist_file}. "
            f"No SNPs were retained after filtering."
        )

    # Read LD matrix (whitespace-separated, no headers)
    # Values can be numbers or 'nan'
    matrix = pd.read_csv(
        ld_file,
        sep=r"\s+",
        header=None,
        names=snp_ids,
        index_col=False,
    )

    # Set SNP IDs as row index
    matrix.index = snp_ids

    return matrix, snp_ids


def parse_ld_output(ld_file: str, lead_snp: str) -> pd.DataFrame:
    """Parse PLINK .ld output file.

    Args:
        ld_file: Path to .ld output file.
        lead_snp: SNP ID of the lead variant.

    Returns:
        DataFrame with columns: SNP, R2.

    Raises:
        PlinkError: If output file is missing after a successful PLINK run.
    """
    if not os.path.exists(ld_file):
        raise PlinkError(
            f"PLINK reported success but output file not found: {ld_file}. "
            f"This may indicate a filesystem issue or PLINK bug."
        )

    # PLINK outputs whitespace-separated: CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
    ld_df = pd.read_csv(ld_file, sep=r"\s+")

    if ld_df.empty:
        # PLINK exited 0 but produced no LD pairs. Most common causes: lead
        # SNP monomorphic, filtered out by --maf, or no variants in window.
        # Silently returning an empty DataFrame would leave the caller with
        # an uncoloured plot and no clue why.
        raise EmptyLDOutputError(
            f"PLINK produced an empty LD output ({ld_file}) for lead SNP "
            f"{lead_snp!r}. The lead SNP may be monomorphic, filtered by "
            f"PLINK quality thresholds, or absent from the reference panel. "
            f"Verify {lead_snp!r} exists in the PLINK .bim file."
        )

    # We want SNP_B (the other SNPs) and their R2 with lead SNP (SNP_A)
    result = ld_df[["SNP_B", "R2"]].rename(columns={"SNP_B": "SNP"})

    # PLINK's --ld-snp output includes the lead paired with itself (R2=1.0).
    # Drop that self-pair so the explicit lead row below is the sole lead entry.
    # Without this the lead key is duplicated, and the validate="many_to_one"
    # LD merge in LocusZoomPlotter.plot() raises MergeError — silently disabling
    # LD colouring for every real PLINK run.
    result = result[result["SNP"] != lead_snp]

    # Add the lead SNP itself with R2=1.0
    lead_row = pd.DataFrame({"SNP": [lead_snp], "R2": [1.0]})
    result = pd.concat([result, lead_row], ignore_index=True)

    return result


def calculate_ld(
    bfile_path: str,
    lead_snp: str,
    window_kb: int = 500,
    plink_path: Optional[str] = None,
    working_dir: Optional[str] = None,
    species: Optional[str] = "canine",
    threads: Optional[int] = None,
) -> pd.DataFrame:
    """Calculate LD (R²) between a lead SNP and all SNPs in a region.

    Runs PLINK --r2 to compute pairwise LD values, then returns a DataFrame
    that can be merged with GWAS results for regional plot coloring.

    Args:
        bfile_path: Path to PLINK binary fileset (.bed/.bim/.fam prefix).
        lead_snp: SNP ID of the lead variant to calculate LD against.
        window_kb: Window size in kilobases around lead SNP.
        plink_path: Path to PLINK executable. Auto-detects if None.
        working_dir: Directory for PLINK output files. Uses temp dir if None.
        species: Species flag ('canine', 'feline', or None for human).
        threads: Number of threads for PLINK.

    Returns:
        DataFrame with columns: SNP (rsid), R2 (LD with lead SNP).

    Raises:
        FileNotFoundError: If PLINK executable not found.
        ValidationError: If PLINK binary files (.bed/.bim/.fam) are missing.
        PlinkError: If PLINK subprocess fails or times out.

    Example:
        >>> ld_df = calculate_ld(
        ...     bfile_path="/path/to/genotypes",
        ...     lead_snp="rs12345",
        ...     window_kb=500,
        ... )
        >>> # Merge with GWAS results for plotting
        >>> gwas_with_ld = gwas_df.merge(ld_df, left_on="rs", right_on="SNP")
    """
    plink_path = _resolve_plink(plink_path)
    validate_plink_files(bfile_path)

    with _plink_workdir(working_dir, "pylocuszoom_ld_") as workdir:
        # Sanitize SNP ID for safe use in file paths (e.g., chr1:12345 contains ':')
        safe_snp_id = re.sub(r"[^\w\-.]", "_", lead_snp)
        if safe_snp_id != lead_snp:
            logger.debug(
                f"Sanitized SNP ID for file path: {lead_snp!r} -> {safe_snp_id!r}"
            )
        output_prefix = os.path.join(workdir, f"ld_{safe_snp_id}")

        cmd = build_ld_command(
            plink_path=plink_path,
            bfile_path=bfile_path,
            lead_snp=lead_snp,
            output_path=output_prefix,
            window_kb=window_kb,
            species=species,
            threads=threads,
        )
        _run_plink(
            cmd,
            workdir,
            "LD calculation",
            "Consider reducing window_kb or checking PLINK installation.",
        )

        return parse_ld_output(f"{output_prefix}.ld", lead_snp)


def calculate_pairwise_ld(
    bfile_path: str,
    snp_list: list[str] | None = None,
    chrom: int | None = None,
    start: int | None = None,
    end: int | None = None,
    plink_path: str | None = None,
    working_dir: str | None = None,
    species: str | None = "canine",
    metric: LDMetric = "r2",
) -> tuple[pd.DataFrame, list[str]]:
    """Calculate pairwise LD matrix for a set of variants.

    Runs PLINK --r2 square to compute an N x N LD matrix, suitable for
    LD heatmap visualization.

    Args:
        bfile_path: Path to PLINK binary fileset (.bed/.bim/.fam prefix).
        snp_list: List of SNP IDs to compute pairwise LD between.
        chrom: Chromosome number for region-based extraction.
        start: Start position (bp) for region-based extraction.
        end: End position (bp) for region-based extraction.
        plink_path: Path to PLINK executable. Auto-detects if None.
        working_dir: Directory for PLINK output files. Uses temp dir if None.
        species: Species flag ('canine', 'feline', or None for human).
        metric: LD metric ('r2' or 'dprime').

    Returns:
        Tuple of (LD matrix DataFrame, list of SNP IDs).
        DataFrame has SNP IDs as both index and columns.

    Raises:
        FileNotFoundError: If PLINK executable not found.
        ValidationError: If PLINK binary files (.bed/.bim/.fam) are missing.
        ValidationError: If requested SNPs are not found in reference panel.
        PlinkError: If PLINK subprocess fails or times out.

    Example:
        >>> matrix, snp_ids = calculate_pairwise_ld(
        ...     bfile_path="/path/to/genotypes",
        ...     snp_list=["rs1", "rs2", "rs3"],
        ... )
        >>> # matrix is 3x3 DataFrame with LD values
        >>> matrix.loc["rs1", "rs2"]  # LD between rs1 and rs2
    """
    _metric_flags(metric)
    plink_path = _resolve_plink(plink_path)
    validate_plink_files(bfile_path)

    with _plink_workdir(working_dir, "pylocuszoom_pairwise_ld_") as workdir:
        output_prefix = os.path.join(workdir, "pairwise_ld")

        snp_list_file = None
        if snp_list:
            snp_list_file = os.path.join(workdir, "snp_list.txt")
            with open(snp_list_file, "w") as f:
                for snp in snp_list:
                    f.write(f"{snp}\n")

        cmd = build_pairwise_ld_command(
            plink_path=plink_path,
            bfile_path=bfile_path,
            output_path=output_prefix,
            snp_list_file=snp_list_file,
            chrom=chrom,
            start=start,
            end=end,
            species=species,
            metric=metric,
        )
        _run_plink(
            cmd,
            workdir,
            "pairwise LD calculation",
            "Consider reducing the region size or SNP count.",
        )

        matrix, found_snps = parse_pairwise_ld_output(
            f"{output_prefix}.ld", f"{output_prefix}.snplist"
        )

        # Validate all requested SNPs were found
        if snp_list:
            missing_snps = set(snp_list) - set(found_snps)
            if missing_snps:
                raise ValidationError(
                    f"SNPs not found in reference panel: {', '.join(sorted(missing_snps))}"
                )

        return matrix, found_snps
