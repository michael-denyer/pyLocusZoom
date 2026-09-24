#!/usr/bin/env python3
"""Execute a notebook headless top to bottom; optionally save its outputs."""

import argparse
import sys
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("notebook", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        help="write the executed notebook here (pass the input path to update it)",
    )
    options = parser.parse_args()
    notebook = nbformat.read(options.notebook, as_version=4)
    client = NotebookClient(
        notebook,
        timeout=300,
        kernel_name="python3",
        resources={"metadata": {"path": str(options.notebook.parent)}},
    )
    try:
        client.execute()
    except CellExecutionError as error:
        print(error, file=sys.stderr)
        return 1
    if options.output:
        nbformat.write(notebook, options.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
