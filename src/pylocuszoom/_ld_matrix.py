"""Shared intake for standalone and regional LD matrices."""

from typing import List, Optional, Union

import numpy as np
import pandas as pd

from .exceptions import ValidationError


def prepare_ld_matrix(
    matrix: Union[pd.DataFrame, np.ndarray], snp_ids: Optional[List[str]]
) -> tuple[np.ndarray, List[str]]:
    """Resolve a square matrix and exactly one SNP id per row and column.

    Missing ids come from a DataFrame's index, or array row numbers.
    Raises ValidationError when the shape and ids do not agree.
    """
    data = matrix.to_numpy() if isinstance(matrix, pd.DataFrame) else np.asarray(matrix)
    if data.ndim != 2 or data.shape[0] != data.shape[1]:
        raise ValidationError(f"ld_matrix must be square, got shape {data.shape}")
    if snp_ids is None:
        snp_ids = (
            list(matrix.index.astype(str))
            if isinstance(matrix, pd.DataFrame)
            else [str(i) for i in range(data.shape[0])]
        )
    if data.shape[0] != len(snp_ids):
        raise ValidationError(
            f"snp_ids length ({len(snp_ids)}) does not match matrix "
            f"dimension ({data.shape[0]})"
        )
    return data, snp_ids
