"""
WPMW shared utilities.

Provides two path helpers for routing file output:

- ``output_path`` — runtime output (Claude container, Kaggle, local scratch).
- ``docs_path``   — figures/data destined for ``docs/`` in the repository,
  to be committed alongside the generating script.

Convention
----------
- ``WPMW_OUTPUT`` environment variable controls the output directory for
  ``output_path``.  If unset, files go to ``./output``.
- ``WPMW_DOCS`` environment variable must be set to the repository root for
  ``docs_path`` to write anything.  If unset, ``docs_path`` returns ``None``
  and callers should skip the write silently.
- Both helpers create their target directory on first call.
- Scripts must use these helpers for *all* generated files (PNGs, MP4s, CSVs,
  patches, etc.).  Hard-coded paths like ``/home/claude/`` or
  ``/mnt/user-data/outputs/`` are forbidden in committed code.

Environment variable summary
-----------------------------
+------------------+------------------------------+---------------------------+
| Context          | WPMW_OUTPUT                  | WPMW_DOCS                 |
+==================+==============================+===========================+
| Claude container | /mnt/user-data/outputs       | (unset)                   |
+------------------+------------------------------+---------------------------+
| Kaggle           | /kaggle/working              | (unset)                   |
+------------------+------------------------------+---------------------------+
| Local (post-patch| (unset → ./output)           | ~/wpmw-output (worktree)  |
+------------------+------------------------------+---------------------------+

Examples
--------
>>> from wpmw_utils import output_path, docs_path
>>> fig.savefig(output_path("phase_portrait.png"), dpi=150, bbox_inches="tight")
>>> dp = docs_path("phase_portrait.png")
>>> if dp:
...     fig.savefig(dp, dpi=150, bbox_inches="tight")
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


def docs_path(filename: str, subdir: str = "figures") -> str | None:
    """Return an absolute path under ``<subdir>/`` in the output branch worktree.

    Only active when the ``WPMW_DOCS`` environment variable is set to the
    root of a git worktree checked out to the orphan ``output`` branch.
    Returns ``None`` when unset so that callers can skip the write gracefully
    in cloud environments.

    Once committed and pushed, files are publicly accessible at::

        https://raw.githubusercontent.com/billpage/wpmw/output/<subdir>/<filename>

    and can be embedded in any markdown doc on ``main`` as::

        ![Alt text](https://raw.githubusercontent.com/billpage/wpmw/output/figures/<filename>)

    Parameters
    ----------
    filename : str
        File name to place inside ``<subdir>/``.
    subdir : str, optional
        Subdirectory under the worktree root to write into.  Defaults to
        ``"figures"``.

    Returns
    -------
    str or None
        Full path when ``WPMW_DOCS`` is set; ``None`` otherwise.
        The target directory is created if it does not exist.

    Example
    -------
    Set up the output branch worktree once::

        cd ~/wpmw
        git worktree add ../wpmw-output output

    Then set ``WPMW_DOCS`` when running a script locally::

        export WPMW_DOCS=~/wpmw-output
        python src/my_script.py

    Then commit the new figure from the worktree::

        cd ~/wpmw-output
        git add figures/phase_portrait.png
        git commit -m "add phase portrait output"
        git push origin output
    """
    worktree_root = os.environ.get("WPMW_DOCS")
    if not worktree_root:
        return None
    dest = os.path.join(worktree_root, subdir)
    os.makedirs(dest, exist_ok=True)
    return os.path.join(dest, filename)
