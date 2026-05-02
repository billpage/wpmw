"""
WPMW shared utilities.

Currently provides only the ``output_path`` helper, which routes file output
through a runtime-configurable directory.

Convention
----------
- ``WPMW_OUTPUT`` environment variable controls the output directory.
- If unset, files go to ``./output``.
- The directory is created on first call.
- Scripts must use this helper for *all* generated files (PNGs, MP4s, CSVs,
  patches, etc.). Hard-coded paths like ``/home/claude/`` or
  ``/mnt/user-data/outputs/`` are forbidden in committed code.

Example
-------
>>> from wpmw_utils import output_path
>>> fig.savefig(output_path("my_figure.png"), dpi=150, bbox_inches="tight")
"""

from __future__ import annotations

import os


def output_path(filename: str) -> str:
    """Return an absolute path under the WPMW output directory.

    Parameters
    ----------
    filename : str
        File name (or relative subpath) to place inside the output directory.

    Returns
    -------
    str
        Full path. The output directory is created if it does not exist.
    """
    base = os.environ.get("WPMW_OUTPUT", "./output")
    os.makedirs(base, exist_ok=True)
    return os.path.join(base, filename)
