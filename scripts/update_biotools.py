"""Send docs/biotools.json to bio.tools as the pyLocusZoom registry entry.

The entry is validated by bio.tools first, so a rejected payload leaves the
live entry as it was. BIOTOOLS_TOKEN is the token that
https://bio.tools/api/rest-auth/login/ returns for the account that owns the
entry.

Usage: BIOTOOLS_TOKEN=<token> python3 scripts/update_biotools.py
"""

import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

ENTRY = Path(__file__).resolve().parents[1] / "docs" / "biotools.json"
API = "https://bio.tools/api/tool"


def put(url: str, payload: bytes, token: str) -> str:
    """PUT the entry to a bio.tools endpoint and return the response body.

    Raises:
        SystemExit: If bio.tools answers with an error status, with its body.
    """
    request = urllib.request.Request(
        url,
        data=payload,
        method="PUT",
        headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
            "Authorization": f"Token {token}",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return response.read().decode()
    except urllib.error.HTTPError as e:
        raise SystemExit(f"{url} answered {e.code}: {e.read().decode()}") from e


def main() -> None:
    token = os.environ.get("BIOTOOLS_TOKEN")
    if not token:
        sys.exit("BIOTOOLS_TOKEN is not set")
    entry = json.loads(ENTRY.read_text())
    payload = json.dumps(entry).encode()
    tool = f"{API}/{entry['biotoolsID']}"
    put(f"{tool}/validate/", payload, token)
    updated = json.loads(put(f"{tool}/", payload, token))
    print(f"bio.tools entry {updated['biotoolsID']} is at {updated['version']}")


if __name__ == "__main__":
    main()
