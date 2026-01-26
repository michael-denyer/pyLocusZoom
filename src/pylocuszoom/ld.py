"""LD (Linkage Disequilibrium) calculation using PLINK.

Calculates R² values between a lead SNP and all other SNPs in a region
using PLINK 1.9's --r2 command.
"""

import os
import shutil
import subprocess
import tempfile
from typing import Optional

import pandas as pd

from .logging import logger
from .utils import validate_plink_files


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


def build_ld_command(
    plink_path: str,
    bfile_path: str,
    lead_snp: str,
    output_path: str,
    window_kb: int = 500,
    ld_window_r2: float = 0.0,
    species: str = "canine",
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

    # Species flag (maps to PLINK's --dog flag)
    if species == "canine":
        cmd.append("--dog")
    elif species == "feline":
        # PLINK doesn't have --cat, use --chr-set for 18 autosomes + X
        cmd.extend(["--chr-set", "18"])

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


def parse_ld_output(ld_file: str, lead_snp: str) -> pd.DataFrame:
    """Parse PLINK .ld output file.

    Args:
        ld_file: Path to .ld output file.
        lead_snp: SNP ID of the lead variant.

    Returns:
        DataFrame with columns: SNP, R2.
    """
    if not os.path.exists(ld_file):
        return pd.DataFrame(columns=["SNP", "R2"])

    # PLINK outputs whitespace-separated: CHR_A BP_A SNP_A CHR_B BP_B SNP_B R2
    ld_df = pd.read_csv(ld_file, sep=r"\s+")

    if ld_df.empty:
        return pd.DataFrame(columns=["SNP", "R2"])

    # We want SNP_B (the other SNPs) and their R2 with lead SNP (SNP_A)
    result = ld_df[["SNP_B", "R2"]].rename(columns={"SNP_B": "SNP"})

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
    species: str = "canine",
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
        Returns empty DataFrame if PLINK fails or no LD values found.

    Raises:
        FileNotFoundError: If PLINK executable not found.
        ValidationError: If PLINK binary files (.bed/.bim/.fam) are missing.

    Example:
        >>> ld_df = calculate_ld(
        ...     bfile_path="/path/to/genotypes",
        ...     lead_snp="rs12345",
        ...     window_kb=500,
        ... )
        >>> # Merge with GWAS results for plotting
        >>> gwas_with_ld = gwas_df.merge(ld_df, left_on="rs", right_on="SNP")
    """
    # Find PLINK first (tests mock this to return None)
    if plink_path is None:
        plink_path = find_plink()
    if plink_path is None:
        raise FileNotFoundError(
            "PLINK not found. Install PLINK 1.9 or specify plink_path."
        )

    logger.debug(f"Using PLINK at {plink_path}")

    # Validate PLINK files exist
    validate_plink_files(bfile_path)

    # Use temp directory if working_dir not specified
    cleanup_working_dir = False
    if working_dir is None:
        working_dir = tempfile.mkdtemp(prefix="snp_scope_ld_")
        cleanup_working_dir = True

    try:
        os.makedirs(working_dir, exist_ok=True)
        output_prefix = os.path.join(working_dir, f"ld_{lead_snp}")

        # Build and run PLINK command
        cmd = build_ld_command(
            plink_path=plink_path,
            bfile_path=bfile_path,
            lead_snp=lead_snp,
            output_path=output_prefix,
            window_kb=window_kb,
            species=species,
            threads=threads,
        )

        logger.debug(f"Running PLINK command: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,
            cwd=working_dir,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            logger.warning(f"PLINK LD calculation failed: {result.stderr[:200]}")
            return pd.DataFrame(columns=["SNP", "R2"])

        # Parse output
        ld_file = f"{output_prefix}.ld"
        return parse_ld_output(ld_file, lead_snp)

    finally:
        # Clean up temp directory
        if cleanup_working_dir and os.path.exists(working_dir):
            shutil.rmtree(working_dir, ignore_errors=True)
