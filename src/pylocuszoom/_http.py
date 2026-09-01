# src/pylocuszoom/_http.py
"""Retrying JSON GET shared by the reference-annotation clients.

Ensembl and UCSC want identical transport behaviour: retry with doubling
backoff on connection errors and on 429/503, raise the client's own error
class on anything else. Only the error class, service name for messages, and
headers differ per caller.
"""

import time
from typing import Any

import requests

from .logging import logger


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
        if response.status_code in (429, 503) and attempt < max_retries - 1:
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
