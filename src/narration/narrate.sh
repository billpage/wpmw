#!/usr/bin/env bash
# Render a narrated slide page locally with Microsoft Edge neural voices.
#
#   bash src/narration/narrate.sh docs/interactive/worlds_and_particles.html [voiceover.py options]
#
# Defaults: the Canadian voice en-CA-LiamNeural, a 30 s pause between Edge
# requests (so the service does not throttle a long deck), the WPMW
# pronunciation file, re-encoding to 32 kbit/s mono (a smaller page), output
# in ~/wpmw-narration/<name>/.  Any option given is passed on to voiceover.py
# and overrides these defaults, for example
#   --preview                 write the spoken text only
#   --reuse-audio             rebuild the page from the clips already made
#   --voice en-US-AndrewNeural  --rate "+8%"
#   --pause 1                 a short pause between Edge requests, for a quick run
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
# The Edge defaults; another --voice or --pause given later wins.  A different
# --tts keeps its own default voice (an Edge voice name means nothing to it).
EDGE_DEFAULTS=(--voice "${WPMW_VOICE:-en-CA-LiamNeural}" --pause "${WPMW_PAUSE:-30}")
case " $* " in
  *" --tts edge "*|*" --tts=edge "*) ;;
  *" --tts "*|*" --tts="*) EDGE_DEFAULTS=() ;;
esac
command -v ffmpeg > /dev/null || { echo "error: ffmpeg is needed for --reencode (sudo apt install ffmpeg)" >&2; exit 1; }

"$HERE/.venv/bin/python" -u "$HERE/voiceover.py" "$SRC" -d "$OUTDIR" \
  ${EDGE_DEFAULTS[@]+"${EDGE_DEFAULTS[@]}"} --pronunciations "$HERE/wpmw_terms.json" --reencode 32 "$@"

echo
echo "Open $OUTDIR/${NAME}_narrated.html in a browser to check it; publish with"
echo "  bash $HERE/publish_narrated.sh $OUTDIR/${NAME}_narrated.html"
