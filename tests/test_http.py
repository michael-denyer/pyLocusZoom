# tests/test_http.py
"""Tests for the shared retrying HTTP transport."""

from unittest.mock import MagicMock, patch

import pytest
import requests

from pylocuszoom._http import download_file
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
