#!/usr/bin/env bash
# Render a narrated slide page locally with Microsoft Edge neural voices.
#
#   bash src/narration/narrate.sh docs/interactive/worlds_and_particles.html [voiceover.py options]
#
# Defaults: the WPMW pronunciation file, re-encoding to 32 kbit/s mono (a
# smaller page), output in ~/wpmw-narration/<name>/.  Any option given is
# passed on to voiceover.py, for example
#   --preview                 write the spoken text only
#   --reuse-audio             rebuild the page from the clips already made
#   --voice en-GB-RyanNeural  --rate "+8%"
#   --tts piper               offline voice (installed on first use), if Edge refuses
#   --pptx deck.pptx          also a narrated PowerPoint
# The first run creates a Python environment in src/narration/.venv.
set -euo pipefail

SRC="${1:?usage: bash narrate.sh slides_page.html [voiceover.py options]}"
shift
HERE="$(cd "$(dirname "$0")" && pwd)"
NAME="$(basename "${SRC%.*}")"
OUTDIR="${WPMW_NARRATION_DIR:-$HOME/wpmw-narration}/$NAME"

if [ ! -x "$HERE/.venv/bin/python" ]; then
  python3 -m venv "$HERE/.venv"
  "$HERE/.venv/bin/pip" install -q edge-tts python-pptx lxml
fi
case " $* " in
  *" --tts piper "*|*" --tts=piper "*)
    "$HERE/.venv/bin/python" -c "import piper" 2> /dev/null || "$HERE/.venv/bin/pip" install -q piper-tts ;;
esac
command -v ffmpeg > /dev/null || { echo "error: ffmpeg is needed for --reencode (sudo apt install ffmpeg)" >&2; exit 1; }

"$HERE/.venv/bin/python" -u "$HERE/voiceover.py" "$SRC" -d "$OUTDIR" \
  --pronunciations "$HERE/wpmw_terms.json" --reencode 32 "$@"

echo
echo "Open $OUTDIR/${NAME}_narrated.html in a browser to check it; publish with"
echo "  bash $HERE/publish_narrated.sh $OUTDIR/${NAME}_narrated.html"
