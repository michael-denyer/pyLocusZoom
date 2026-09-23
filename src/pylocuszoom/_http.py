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
import tempfile
import time
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any, TypeVar

import requests
from tqdm import tqdm

from .exceptions import DataDownloadError
from .logging import logger

RETRYABLE_STATUS = (429, 503)

T = TypeVar("T")


def _retryable(error: requests.RequestException) -> bool:
    """Retry connection failures, and HTTP errors only on 429 or 503."""
    if isinstance(error, requests.HTTPError):
        return _status_of(error) in RETRYABLE_STATUS
    return True


def _status_of(error: requests.HTTPError) -> int | None:
    """Read the status code off an HTTPError, or None if it carries no response."""
    response = getattr(error, "response", None)
    return getattr(response, "status_code", None)


def _with_retries(
    attempt: Callable[[], T], *, what: str, max_retries: int, retry_delay: float
) -> T:
    """Run ``attempt``, retrying retryable request errors with doubling backoff.

    The one retry loop: a JSON GET and a streamed download differ only in
    what one attempt does.

    Raises:
        requests.RequestException: The last attempt's error, or the first
            error that is not worth retrying.
    """
    delay = retry_delay
    attempt_number = 1
    while True:
        try:
            return attempt()
        except requests.RequestException as e:
            if attempt_number >= max_retries or not _retryable(e):
                raise
            logger.warning(f"{what} failed (attempt {attempt_number}): {e}")
            time.sleep(delay)
            delay *= 2
            attempt_number += 1


@contextmanager
def staged_path(dest: Path) -> Iterator[Path]:
    """Yield a private sibling of ``dest``, published over it only on success.

    The one temp-file-then-replace writer. Concurrent writers never share an
    in-progress file, a reader never sees a partial one, and the sibling is
    removed whether or not the block succeeds.
    """
    with tempfile.NamedTemporaryFile(
        dir=dest.parent, prefix=f".{dest.name}.", suffix=".part", delete=False
    ) as partial:
        partial_path = Path(partial.name)
    try:
        yield partial_path
        os.replace(partial_path, dest)
    finally:
        partial_path.unlink(missing_ok=True)


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

    def get() -> requests.Response:
        response = requests.get(url, params=params, headers=headers, timeout=timeout)
        if not response.ok:
            raise requests.HTTPError(
                f"{service} API error {response.status_code}: {response.text[:200]}",
                response=response,
            )
        return response

    try:
        response = _with_retries(
            get,
            what=f"{service} API request",
            max_retries=max_retries,
            retry_delay=retry_delay,
        )
    except requests.HTTPError as e:
        raise error_cls(str(e)) from e
    except requests.RequestException as e:
        raise error_cls(
            f"{service} API request failed after {max_retries} attempts: {e}"
        ) from e

    try:
        return response.json()
    except ValueError as e:  # requests' JSONDecodeError is a ValueError
        raise error_cls(f"{service} API returned invalid JSON: {e}") from e


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

    Streams into a private ``.part`` sibling and replaces ``dest_path`` only once
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
    with staged_path(dest_path) as partial_path:
        try:
            _with_retries(
                lambda: _stream_to(url, partial_path, desc, timeout),
                what=f"Download of {url}",
                max_retries=max_retries,
                retry_delay=retry_delay,
            )
        except requests.RequestException as e:
            raise DataDownloadError(f"Failed to download {url}: {e}") from e


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
