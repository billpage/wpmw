"""Run a wpmw script on Kaggle as a private batch kernel.

A kernel here is one directory holding ``run.py`` and
``kernel-metadata.json``.  ``run.py`` clones the public repository, checks
out a pinned commit, optionally applies a patch (changes not yet on main),
sets ``WPMW_OUTPUT=/kaggle/working``, runs ``src/<script>`` with the given
arguments and tees the output to ``run.log``.  Everything the script writes
through ``output_path()`` comes back with ``fetch``.

Requires the ``kaggle`` client (``pip install kaggle``; 2.x accepts an API
token from kaggle.com -> Settings -> API -> Generate New Token, saved as
``~/.kaggle/access_token`` with mode 600 or exported as KAGGLE_API_TOKEN).

Command line::

    python -m wpmwlib.kaggle_batch push   SLUG SCRIPT [--commit REF] [--patch FILE |
                                          --worktree-diff] [--gpu] -- ARGS...
    python -m wpmwlib.kaggle_batch status SLUG
    python -m wpmwlib.kaggle_batch fetch  SLUG [DEST]

Example, a sweep pinned to the current origin/main::

    PYTHONPATH=src python3 -m wpmwlib.kaggle_batch push wpmw-resonance-d \\
        demo_sea_resonance_clock.py --commit origin/main -- \\
        --parts D --nu 32 64 --seeds 8 --tag _sweep
    PYTHONPATH=src python3 -m wpmwlib.kaggle_batch status wpmw-resonance-d
    PYTHONPATH=src python3 -m wpmwlib.kaggle_batch fetch wpmw-resonance-d

``--commit`` accepts anything ``git rev-parse`` resolves; the full SHA is
written into the kernel, and pushing refuses a commit that no remote branch
contains, since Kaggle can only clone what is on GitHub.  ``--worktree-diff``
embeds ``git diff <commit>`` of the local working tree, so uncommitted
changes can be tested before they are applied upstream (a new file appears in
that diff only once tracked: ``git add -N <file>``).  Slugs are lower
case letters, digits and hyphens, 5 to 50 characters (Kaggle derives the
slug from the title, so the title is the slug with spaces).
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys

from wpmwlib.wpmw_utils import output_path

REPO_URL = "https://github.com/billpage/wpmw.git"
OUT_VAR = "WPMW_OUTPUT"
SLUG_RE = re.compile(r"^[a-z0-9][a-z0-9-]{3,48}[a-z0-9]$")

RUN_TEMPLATE = '''"""Kaggle batch kernel for {repo} (written by wpmwlib.kaggle_batch)."""
import os, subprocess, sys, time

REPO_URL = {repo!r}
COMMIT = {commit!r}
SCRIPT = {script!r}
ARGS = {args!r}
PATCH = {patch!r}

os.environ[{out_var!r}] = "/kaggle/working"
os.environ.pop("WPMW_DOCS", None)
repo = "/tmp/repo"                      # outside /kaggle/working: results only there
subprocess.run(["git", "clone", "-q", REPO_URL, repo], check=True)
subprocess.run(["git", "checkout", "-q", COMMIT], cwd=repo, check=True)
if PATCH:
    with open("/tmp/local.patch", "w") as fh:
        fh.write(PATCH)
    subprocess.run(["git", "apply", "/tmp/local.patch"], cwd=repo, check=True)
    subprocess.run(["git", "diff", "--stat"], cwd=repo, check=True)
subprocess.run(["git", "log", "--oneline", "-1"], cwd=repo, check=True)
print("cpu_count", os.cpu_count(), flush=True)
src = os.path.join(repo, "src")
env = dict(os.environ, PYTHONPATH=src)
cmd = [sys.executable, "-u", os.path.join(src, SCRIPT)] + ARGS
print(">>>", " ".join(cmd), flush=True)
t0 = time.time()
with open("/kaggle/working/run.log", "w") as log:
    proc = subprocess.Popen(cmd, cwd=src, env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, bufsize=1)
    for line in proc.stdout:
        print(line, end="", flush=True)
        log.write(line)
        log.flush()
    proc.wait()
print(f"exit {{proc.returncode}} after {{time.time() - t0:.0f}} s", flush=True)
sys.exit(proc.returncode)
'''


def _git(*args: str, raw: bool = False) -> str:
    out = subprocess.run(["git", *args], check=True, capture_output=True,
                         text=True).stdout
    return out if raw else out.strip()   # never strip a diff: blank context lines


def _kaggle(*args: str, capture: bool = True) -> str:
    res = subprocess.run(["kaggle", *args], capture_output=capture, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"kaggle {' '.join(args)} failed:\n{res.stdout}{res.stderr}")
    return res.stdout if capture else ""


def username() -> str:
    """The Kaggle account the configured credentials belong to."""
    m = re.search(r"username:\s*(\S+)", _kaggle("config", "view"))
    if not m or m.group(1) == "None":
        raise RuntimeError("no Kaggle credentials: save an API token to "
                           "~/.kaggle/access_token (mode 600)")
    return m.group(1)


def resolve_commit(ref: str) -> str:
    """Full SHA of ``ref``; refuses a commit no remote branch contains."""
    sha = _git("rev-parse", "--verify", f"{ref}^{{commit}}")
    if not _git("branch", "-r", "--contains", sha):
        raise RuntimeError(f"{ref} ({sha[:7]}) is on no remote branch, so Kaggle "
                           "cannot clone it; push it first or pin an upstream "
                           "commit and pass the change with --worktree-diff")
    return sha


def check_patch(patch: str, sha: str) -> None:
    """Refuse a patch that does not apply to commit ``sha`` (checked against
    that commit's tree in a scratch index, so the working tree is untouched)."""
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        env = dict(os.environ, GIT_INDEX_FILE=os.path.join(tmp, "index"))
        subprocess.run(["git", "read-tree", sha], env=env, check=True)
        res = subprocess.run(["git", "apply", "--check", "--cached", "-"], input=patch,
                             env=env, text=True, capture_output=True)
    if res.returncode != 0:
        raise RuntimeError(f"patch does not apply to {sha[:7]}:\n{res.stderr}")


def write_kernel(slug: str, script: str, args: list[str], commit: str, *,
                 user: str, patch: str = "", gpu: bool = False,
                 dest: str | None = None) -> str:
    """Write run.py and kernel-metadata.json; return the kernel directory."""
    if not SLUG_RE.match(slug):
        raise ValueError(f"slug {slug!r}: use 5-50 lower-case letters, digits, hyphens")
    path = dest or output_path(os.path.join("kaggle_kernels", slug))
    os.makedirs(path, exist_ok=True)
    with open(os.path.join(path, "run.py"), "w") as fh:
        fh.write(RUN_TEMPLATE.format(repo=REPO_URL, commit=commit, script=script,
                                     args=list(args), patch=patch, out_var=OUT_VAR))
    meta = {"id": f"{user}/{slug}", "title": slug.replace("-", " "),
            "code_file": "run.py", "language": "python", "kernel_type": "script",
            "is_private": True, "enable_gpu": bool(gpu), "enable_tpu": False,
            "enable_internet": True, "dataset_sources": [],
            "competition_sources": [], "kernel_sources": []}
    with open(os.path.join(path, "kernel-metadata.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    return path


def push(path: str) -> str:
    return _kaggle("kernels", "push", "-p", path).strip()


def status(ref: str) -> str:
    """'running', 'complete', 'error', ... for user/slug."""
    out = _kaggle("kernels", "status", ref)
    m = re.search(r'status "(?:KernelWorkerStatus\.)?(\w+)"', out)
    return m.group(1).lower() if m else out.strip()


def fetch(ref: str, dest: str) -> list[str]:
    """Download the kernel's output files and log into dest."""
    os.makedirs(dest, exist_ok=True)
    _kaggle("kernels", "output", ref, "-p", dest)
    return sorted(os.listdir(dest))


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    script_args: list[str] = []
    if "--" in argv:
        i = argv.index("--")
        argv, script_args = argv[:i], argv[i + 1:]
    ap = argparse.ArgumentParser(prog="python -m wpmwlib.kaggle_batch",
                                 description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("push", help="write and push a kernel")
    p.add_argument("slug")
    p.add_argument("script", help="file under src/, e.g. demo_sea_resonance_clock.py")
    p.add_argument("--commit", default="origin/main")
    g = p.add_mutually_exclusive_group()
    g.add_argument("--patch", help="patch file applied after checkout")
    g.add_argument("--worktree-diff", action="store_true",
                   help="embed `git diff <commit>` of the local working tree")
    p.add_argument("--gpu", action="store_true")
    p.add_argument("--dry-run", action="store_true", help="write the kernel, do not push")
    s = sub.add_parser("status", help="status of user/slug or slug")
    s.add_argument("slug")
    f = sub.add_parser("fetch", help="download outputs of user/slug or slug")
    f.add_argument("slug")
    f.add_argument("dest", nargs="?")
    a = ap.parse_args(argv)

    if a.cmd == "push":
        sha = resolve_commit(a.commit)
        patch = ""
        if a.patch:
            with open(a.patch) as fh:
                patch = fh.read()
        elif a.worktree_diff:
            patch = _git("diff", sha, raw=True)
        if patch:
            check_patch(patch, sha)
        user = username()
        path = write_kernel(a.slug, a.script, script_args, sha, user=user,
                            patch=patch, gpu=a.gpu)
        print(f"kernel {user}/{a.slug} at {path}: {a.script} {' '.join(script_args)}"
              f" @ {sha[:7]}" + (f" + patch ({len(patch.splitlines())} lines)"
                                 if patch else ""))
        if not a.dry_run:
            print(push(path))
        return 0
    ref = a.slug if "/" in a.slug else f"{username()}/{a.slug}"
    if a.cmd == "status":
        print(f"{ref}: {status(ref)}")
        return 0
    dest = a.dest or output_path(os.path.join("kaggle_out", ref.split("/")[1]))
    for name in fetch(ref, dest):
        print(os.path.join(dest, name))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (RuntimeError, ValueError, subprocess.CalledProcessError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
