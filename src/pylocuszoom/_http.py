# src/pylocuszoom/_http.py
"""Retrying HTTP transport shared by everything in the library that fetches.

Ensembl and UCSC want identical behaviour for a JSON GET: retry with doubling
backoff on connection errors and on 429/503, raise the client's own error class
on anything else. Only the error class, service name for messages, and headers
differ per caller.

``download_file`` streams a large file over the same retry policy, so the 50 MB
recombination tarball gets the attempts the 5 KB JSON payload always had.
"""

import os
import time
from pathlib import Path
from typing import Any

import requests
from tqdm import tqdm

from .exceptions import DataDownloadError
from .logging import logger

RETRYABLE_STATUS = (429, 503)


def request_json(
    url: str,
    params: dict,
    *,
    error_cls: type[Exception],
    service: str,
    headers: dict | None = None,
    timeout: float = 30,
    max_retries: int = 3,
    retry_delay: float = 1.0,
) -> Any:
    """GET a JSON payload, retrying on connection errors, 429 and 503.

    Always raises on failure; callers that want an empty result instead
    translate ``error_cls`` at their boundary.

    Args:
        url: Endpoint URL.
        params: Query parameters.
        error_cls: Exception class raised on failure.
        service: Service name used in error and log messages.
        headers: Optional request headers.
        timeout: Per-request timeout in seconds.
        max_retries: Attempts before giving up on a retryable error.
        retry_delay: Initial backoff in seconds; doubles on each retry.

    Returns:
        The decoded JSON payload.

    Raises:
        error_cls: If the request ultimately fails.
    """
    delay = retry_delay

    for attempt in range(max_retries):
        try:
            response = requests.get(
                url, params=params, headers=headers, timeout=timeout
            )
        except requests.RequestException as e:
            logger.warning(f"{service} API request failed (attempt {attempt + 1}): {e}")
            if attempt < max_retries - 1:
                time.sleep(delay)
                delay *= 2
                continue
            raise error_cls(
                f"{service} API request failed after {max_retries} attempts: {e}"
            )

        if response.ok:
            try:
                return response.json()
            except (ValueError, requests.exceptions.JSONDecodeError) as e:
                logger.warning(f"{service} API returned invalid JSON: {e}")
                raise error_cls(f"{service} API returned invalid JSON: {e}")

        # Retryable errors (429 rate limit, 503 service unavailable)
        if response.status_code in RETRYABLE_STATUS and attempt < max_retries - 1:
            logger.warning(
                f"{service} API returned {response.status_code} "
                f"(attempt {attempt + 1}), retrying..."
            )
            time.sleep(delay)
            delay *= 2
            continue

        error_msg = f"{service} API error {response.status_code}: {response.text[:200]}"
        logger.warning(error_msg)
        raise error_cls(error_msg)

    # All retries exhausted (e.g., repeated 429/503 responses)
    raise error_cls(
        f"{service} API request failed after {max_retries} attempts (rate limited)"
    )


def download_file(
    url: str,
    dest_path: Path,
    desc: str = "Downloading",
    *,
    timeout: float = 60,
    max_retries: int = 3,
    retry_delay: float = 1.0,
) -> None:
    """Stream a file to disk with a progress bar, retrying like request_json.

    Streams into a ``.part`` sibling and renames it onto ``dest_path`` only once
    the stream completes, so an interrupted download never leaves a truncated
    file where a later ``exists()`` check would trust it. A retry restarts the
    stream from the beginning; the partial file is removed either way.

    Connection failures retry, as they do for a JSON request. An HTTP status
    retries only when it is 429 or 503; a 404 is not going to become a file on
    the second attempt.

    Args:
        url: URL to download from.
        dest_path: Destination file path.
        desc: Description for the progress bar.
        timeout: Per-request timeout in seconds.
        max_retries: Attempts before giving up on a retryable error.
        retry_delay: Initial backoff in seconds; doubles on each retry.

    Raises:
        DataDownloadError: If the download ultimately fails.
    """
    partial_path = dest_path.with_name(dest_path.name + ".part")
    delay = retry_delay

    for attempt in range(max_retries):
        try:
            _stream_to(url, partial_path, desc, timeout)
        except requests.RequestException as e:
            partial_path.unlink(missing_ok=True)
            retryable = (
                not isinstance(e, requests.HTTPError)
                or _status_of(e) in RETRYABLE_STATUS
            )
            if retryable and attempt < max_retries - 1:
                logger.warning(f"Download of {url} failed (attempt {attempt + 1}): {e}")
                time.sleep(delay)
                delay *= 2
                continue
            raise DataDownloadError(f"Failed to download {url}: {e}") from e
        except BaseException:
            partial_path.unlink(missing_ok=True)
            raise

        os.replace(partial_path, dest_path)
        return


def _status_of(error: requests.HTTPError) -> int | None:
    """Read the status code off an HTTPError, or None if it carries no response."""
    response = getattr(error, "response", None)
    return getattr(response, "status_code", None)


def _stream_to(url: str, partial_path: Path, desc: str, timeout: float) -> None:
    """Stream one response body into partial_path, showing a progress bar."""
    response = requests.get(url, stream=True, timeout=timeout)
    response.raise_for_status()
    total_size = int(response.headers.get("content-length", 0))
    with (
        open(partial_path, "wb") as f,
        tqdm(
            total=total_size,
            unit="B",
            unit_scale=True,
            unit_divisor=1024,
            desc=desc,
            disable=total_size == 0,  # Disable if size unknown
        ) as pbar,
    ):
        for chunk in response.iter_content(chunk_size=8192):
            if chunk:
                f.write(chunk)
                pbar.update(len(chunk))
