#!/usr/bin/env bash
# Publish a narrated slide page that you rendered and checked locally.
#
#   bash src/narration/publish_narrated.sh ~/worlds_and_particles_narration/worlds_and_particles_narrated.html
#
# Copies the single-file page to narrated/ on the output branch (worktree
# ~/wpmw-output), commits and pushes it, then starts the Pages workflow so it
# is served at https://billpage.github.io/wpmw/narrated/<file>.
# Nothing is added to main: main keeps only the inputs (the slides page in
# docs/interactive/ with its notes, and src/narration/).
set -euo pipefail

PAGE="${1:?usage: bash publish_narrated.sh path/to/<name>_narrated.html}"
OUT="${WPMW_OUTPUT_WORKTREE:-$HOME/wpmw-output}"
REPO_SLUG="billpage/wpmw"

[ -f "$PAGE" ] || { echo "error: $PAGE not found" >&2; exit 1; }
case "$PAGE" in *_narrated.html) ;; *) echo "error: expected a *_narrated.html page from voiceover.py" >&2; exit 1;; esac
grep -q 'id="nv-audio"' "$PAGE" || { echo "error: $PAGE is not a narrated page" >&2; exit 1; }
if ! grep -q 'const MISSING = \[\];' "$PAGE"; then
  echo "error: some slides in $PAGE have no audio; run narrate.sh again to fill them in" >&2; exit 1
fi
if grep -q 'const EMBEDDED = false' "$PAGE"; then
  echo "error: this page links its audio (--link); publish the embedded version" >&2; exit 1
fi
[ -d "$OUT/.git" ] || [ -f "$OUT/.git" ] || {
  echo "error: no output worktree at $OUT; create it with: cd ~/wpmw && git worktree add ../wpmw-output output" >&2; exit 1; }

SIZE_MB=$(( $(stat -c %s "$PAGE") / 1000000 ))
if [ "$SIZE_MB" -ge 95 ]; then
  echo "error: $PAGE is ${SIZE_MB} MB; GitHub refuses files over 100 MB (rebuild with --reencode 24)" >&2; exit 1
elif [ "$SIZE_MB" -ge 50 ]; then
  echo "warning: $PAGE is ${SIZE_MB} MB; consider --reencode 32 to keep the output branch small"
fi

NAME="$(basename "$PAGE")"
cd "$OUT"
git pull --ff-only origin output
mkdir -p narrated
cp "$PAGE" "narrated/$NAME"
git add "narrated/$NAME"
if git diff --cached --quiet; then
  echo "narrated/$NAME is unchanged on output; nothing to commit"
else
  git commit -q -m "narrated: $NAME ($(date +%Y-%m-%d))"
  git push origin output
  git log --oneline -1
fi

if command -v gh > /dev/null && gh auth status > /dev/null 2>&1; then
  gh workflow run pages.yml --ref main -R "$REPO_SLUG"
  echo "Pages workflow started; it takes a minute or two."
else
  echo "Start the Pages deployment: https://github.com/$REPO_SLUG/actions/workflows/pages.yml -> Run workflow"
fi
echo
echo "Play online : https://billpage.github.io/wpmw/narrated/$NAME"
echo "Single file : the page's Save button, or https://github.com/$REPO_SLUG/blob/output/narrated/$NAME (Download raw file)"
