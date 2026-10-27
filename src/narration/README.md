# Narrated slide pages

Tools for turning a slides page in `docs/interactive/` into a narrated page:
the slides, an audio player per slide and the speaker notes as a transcript,
in one HTML file that plays online and also works when downloaded.

- [`voiceover.py`](voiceover.py) — reads the speaker notes, makes one MP3 per
  slide with a text-to-speech voice, and writes the narrated page (and,
  optionally, a PowerPoint with click-to-play audio). It also accepts a
  `.pptx` as the source. `python3 voiceover.py --help` lists every option.
- [`wpmw_terms.json`](wpmw_terms.json) — how to say this project's notation
  and words: `{"written": "spoken"}`, applied before synthesis.
- [`narrate.sh`](narrate.sh) — runs `voiceover.py` with the project's
  defaults (Edge voice, the terms above, 32 kbit/s audio).
- [`publish_narrated.sh`](publish_narrated.sh) — puts a page you have checked
  on the `output` branch and starts the Pages deployment.

**What is kept where.** `main` keeps only the inputs: the slides page, whose
`NOTES` are the transcript, and this folder. The rendered page, with its
audio embedded (about 15 to 20 MB for an hour of speech), goes to
`narrated/` on the `output` branch and is served at
`https://billpage.github.io/wpmw/narrated/<name>_narrated.html` by the
[`pages`](../../.github/workflows/pages.yml) workflow. The clips stay on your
machine (in `~/wpmw-narration/<name>/audio/`), so a page can be rebuilt
without synthesising again.

## Workflow

1. **Check the text first** (no audio, a few seconds):

   ```bash
   bash src/narration/narrate.sh docs/interactive/worlds_and_particles.html --preview
   ```

   Read `~/wpmw-narration/worlds_and_particles/worlds_and_particles_spoken.txt`.
   Any slide marked `NOT SPOKEN WELL` has characters the voice will stumble on.

2. **Render** (Edge voice, needs internet; about a minute per ten minutes of speech):

   ```bash
   bash src/narration/narrate.sh docs/interactive/worlds_and_particles.html
   ```

   Open `~/wpmw-narration/worlds_and_particles/worlds_and_particles_narrated.html`
   in a browser and listen.

   Microsoft's voice service sometimes stalls or refuses a request it accepts
   a minute later. The script asks for at most about 1200 characters at a time,
   abandons a request after 60 s (`--timeout`), and retries a failed one three
   times, 5, 20 and 60 s apart. If a slide still fails, the run ends with
   `MISSING AUDIO for slide(s) ...` and exit status 3, and the publish script
   refuses the page. Run the same command again: a slide whose spoken text,
   voice and rate are unchanged keeps its clip (`audio/slide_NN.txt` records
   what it was made from), so only the missing slides are synthesised.
   `--resynth` makes every clip again.

   If many first tries are refused and then go through on the retry (the
   service limits how often requests arrive), add `--pause 5`: the script then
   waits 5 s instead of 1 s between requests. With about 63 requests for an
   hour of speech that adds about four minutes and may save the retries.

   If even the first slide fails, Edge is refusing this machine for the time
   being (it does this after heavy use), and the run stops after two slides
   rather than retrying for hours. Check with
   `src/narration/.venv/bin/edge-tts --text test --write-media /tmp/t.mp3`;
   update with `src/narration/.venv/bin/pip install -U edge-tts`; or wait an
   hour or two. `--tts piper` uses an offline neural voice instead (Piper,
   installed and downloaded on first use, about 60 MB): steadier, a little less
   natural, and a different voice, so use it for a whole deck rather than for
   a few missing slides.

3. **Fix what sounds wrong** in `wpmw_terms.json` and render again. A word
   whose pronunciation depends on context (*winds*, *live*, *read*) needs the
   surrounding words in the key, because keys are matched as plain text,
   even inside other words: `"clock winds at": "clock wynds at"`, never
   `"wind": ...`. Only the slides whose spoken text changed are synthesised
   again; `--reuse-audio` rebuilds the page without contacting the voice
   service at all.

4. **Publish**:

   ```bash
   bash src/narration/publish_narrated.sh ~/wpmw-narration/worlds_and_particles/worlds_and_particles_narrated.html
   ```

   It commits the page to `output`, pushes, and starts the Pages workflow
   (with the GitHub CLI `gh` if it is logged in; otherwise it prints the link
   to the workflow's *Run workflow* button).

## The page

Left and right arrow keys change slide; Space or K plays and pauses; J and L
go back and forward 10 s; T shows or hides the transcript; C turns on
continuous play, which moves to the next slide at the end of each clip.
Drag the divider between slide and transcript to resize it, or focus it and
use the arrow keys; double-click resets it. The page remembers speed,
divider, transcript and continuous play in the browser. *Save* downloads the
page as one file that plays offline.

## Requirements

Python 3.9 or later, `ffmpeg`, and internet access for the Edge voice;
`narrate.sh` installs `edge-tts`, `python-pptx` and `lxml` into
`src/narration/.venv` on first use. A `.pptx` source also needs LibreOffice
and `pdftoppm` (poppler-utils). `voiceover.py --tts silent` makes silent clips
of the right length, for checking layout without a voice.
