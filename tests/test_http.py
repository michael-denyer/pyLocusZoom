# tests/test_http.py
"""Tests for the shared retrying HTTP transport."""

from unittest.mock import MagicMock, patch

import pytest
import requests
from hypothesis import given
from hypothesis import strategies as st

from pylocuszoom._http import _with_retries, download_file
from pylocuszoom.exceptions import DataDownloadError


class TestDownloadFile:
    """A download must never leave a truncated file at the destination."""

    @staticmethod
    def _streaming_response(chunks):
        response = MagicMock()
        response.headers = {"content-length": str(sum(len(c) for c in chunks))}
        response.iter_content.return_value = iter(chunks)
        return response

    def test_complete_download_lands_at_dest(self, tmp_path):
        dest = tmp_path / "file.gz"
        response = self._streaming_response([b"abcd", b"efgh"])

        with patch("pylocuszoom._http.requests.get", return_value=response):
            download_file("https://example.invalid/f", dest)

        assert dest.read_bytes() == b"abcdefgh"
        assert [p.name for p in tmp_path.iterdir()] == ["file.gz"]

    def test_interrupted_download_leaves_nothing_behind(self, tmp_path):
        dest = tmp_path / "file.gz"

        def interrupted_response(*_args, **_kwargs):
            def chunks():
                yield b"abcd"
                raise requests.ConnectionError("connection reset")

            response = self._streaming_response([])
            response.iter_content.return_value = chunks()
            return response

        with (
            patch("pylocuszoom._http.time.sleep"),
            patch("pylocuszoom._http.requests.get", side_effect=interrupted_response),
            pytest.raises(DataDownloadError, match="example.invalid/f") as exc_info,
        ):
            download_file("https://example.invalid/f", dest)

        assert isinstance(exc_info.value.__cause__, requests.ConnectionError)
        assert not dest.exists()
        assert list(tmp_path.iterdir()) == []

    def test_connection_failure_retries(self, tmp_path):
        """A dropped connection retries, unlike the download this replaced."""
        dest = tmp_path / "file.gz"
        failed = self._streaming_response([])
        failed.iter_content.side_effect = requests.ConnectionError("reset")

        with (
            patch("pylocuszoom._http.time.sleep"),
            patch(
                "pylocuszoom._http.requests.get",
                side_effect=[failed, self._streaming_response([b"abcd"])],
            ) as mock_get,
        ):
            download_file("https://example.invalid/f", dest)

        assert mock_get.call_count == 2
        assert dest.read_bytes() == b"abcd"

    def test_http_error_does_not_retry(self, tmp_path):
        """A 404 will not become a file on the second attempt."""
        dest = tmp_path / "file.gz"
        response = self._streaming_response([])
        error = requests.HTTPError("404 Client Error")
        error.response = MagicMock(status_code=404)
        response.raise_for_status.side_effect = error

        with (
            patch("pylocuszoom._http.time.sleep"),
            patch("pylocuszoom._http.requests.get", return_value=response) as mock_get,
            pytest.raises(DataDownloadError),
        ):
            download_file("https://example.invalid/f", dest)

        assert mock_get.call_count == 1

    def test_http_error_raises_download_error(self, tmp_path):
        dest = tmp_path / "file.gz"
        response = self._streaming_response([])
        original = requests.HTTPError("404 Client Error")
        response.raise_for_status.side_effect = original

        with (
            patch("pylocuszoom._http.requests.get", return_value=response),
            pytest.raises(DataDownloadError, match="example.invalid/f") as exc_info,
        ):
            download_file("https://example.invalid/f", dest)

        assert exc_info.value.__cause__ is original
        assert list(tmp_path.iterdir()) == []


def test_concurrent_downloads_publish_only_their_own_complete_response(
    tmp_path, monkeypatch
):
    from concurrent.futures import ThreadPoolExecutor
    from threading import Event

    first_written, second_written, first_published = Event(), Event(), Event()
    dest = tmp_path / "shared.gz"

    def stream(url, partial, desc, timeout):
        if url == "first":
            partial.write_bytes(b"first")
            first_written.set()
            assert second_written.wait(5)
        else:
            assert first_written.wait(5)
            partial.write_bytes(b"second")
            second_written.set()
            assert first_published.wait(5)

    monkeypatch.setattr("pylocuszoom._http._stream_to", stream)
    with ThreadPoolExecutor(max_workers=2) as executor:
        first = executor.submit(download_file, "first", dest)
        second = executor.submit(download_file, "second", dest)
        try:
            first.result(timeout=10)
            first_result = dest.read_bytes()
        finally:
            first_published.set()
        second.result(timeout=10)
    assert first_result == b"first"
    assert dest.read_bytes() == b"second"
    assert list(tmp_path.iterdir()) == [dest]


class TestRequestJson:
    def test_invalid_json_is_the_callers_error_with_the_cause_kept(self):
        from pylocuszoom._http import request_json
        from pylocuszoom.exceptions import EnsemblAPIError

        response = MagicMock(ok=True)
        response.json.side_effect = ValueError("Expecting value")

        with (
            patch("pylocuszoom._http.requests.get", return_value=response),
            pytest.raises(EnsemblAPIError, match="invalid JSON") as exc_info,
        ):
            request_json(
                "https://example.invalid", {}, error_cls=EnsemblAPIError, service="X"
            )

        assert isinstance(exc_info.value.__cause__, ValueError)

    def test_a_503_retries_and_then_succeeds(self):
        from pylocuszoom._http import request_json
        from pylocuszoom.exceptions import EnsemblAPIError

        busy = MagicMock(ok=False, status_code=503, text="busy")
        ok = MagicMock(ok=True)
        ok.json.return_value = {"answer": 42}

        with (
            patch("pylocuszoom._http.time.sleep") as sleep,
            patch("pylocuszoom._http.requests.get", side_effect=[busy, busy, ok]),
        ):
            payload = request_json(
                "https://example.invalid", {}, error_cls=EnsemblAPIError, service="X"
            )

        assert payload == {"answer": 42}
        assert [call.args[0] for call in sleep.call_args_list] == [1.0, 2.0]


def _http_error(status):
    response = MagicMock(status_code=status)
    return requests.HTTPError(str(status), response=response)


_retryable_errors = st.sampled_from(
    [requests.ConnectionError("reset"), requests.Timeout("slow")]
) | st.sampled_from([429, 503]).map(_http_error)
_fatal_errors = st.sampled_from([400, 401, 403, 404, 500, 502]).map(_http_error)


class TestRetryProperties:
    """The one retry loop behind every JSON GET and download."""

    @staticmethod
    def _run(outcomes, max_retries, retry_delay=1.0):
        """Run the loop over a script of errors then success; return attempts, sleeps."""
        script = iter(outcomes)
        attempts = []

        def attempt():
            attempts.append(1)
            outcome = next(script, "ok")
            if isinstance(outcome, Exception):
                raise outcome
            return outcome

        with patch("pylocuszoom._http.time.sleep") as sleep:
            try:
                result = _with_retries(
                    attempt,
                    what="GET",
                    max_retries=max_retries,
                    retry_delay=retry_delay,
                )
            except requests.RequestException as e:
                result = e
        return result, len(attempts), [call.args[0] for call in sleep.call_args_list]

    @given(st.lists(_retryable_errors, max_size=6), st.integers(1, 6))
    def test_retryable_errors_retry_up_to_the_limit(self, errors, max_retries):
        result, attempts, sleeps = self._run(errors, max_retries)

        assert attempts == min(len(errors) + 1, max_retries)
        if len(errors) < max_retries:
            assert result == "ok"
        else:
            assert result is errors[max_retries - 1]
        assert sleeps == [2.0**k for k in range(attempts - 1)]

    @given(st.lists(_retryable_errors, max_size=4), _fatal_errors, st.integers(1, 6))
    def test_a_fatal_error_is_raised_without_another_attempt(
        self, retryable, fatal, max_retries
    ):
        result, attempts, _ = self._run([*retryable, fatal], max_retries)

        if len(retryable) < max_retries:
            assert result is fatal
            assert attempts == len(retryable) + 1
        else:
            assert attempts == max_retries
