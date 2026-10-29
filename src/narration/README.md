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
  defaults: the Canadian Edge voice `en-CA-LiamNeural`, a 30 s pause between
  requests to the voice service, the terms above, and 32 kbit/s audio. Any
  option you add overrides them (`--voice en-US-AndrewNeural`, `--pause 1`);
  so do the environment variables `WPMW_VOICE` and `WPMW_PAUSE`.
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

2. **Render** (Edge voice, needs internet; with the 30 s pause, about half an
   hour for an hour of speech, almost all of it waiting):

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
   `--resynth` makes every clip again; `--resynth 14` or `--resynth 5,14`
   makes only those.

   The service limits how often requests arrive, and refuses many first tries
   when they come quickly, so `narrate.sh` waits 30 s after each request
   (`--pause 30`; `voiceover.py` on its own waits 1 s). The pause is time
   between requests, not silence in the audio. An hour of speech is about 63
   requests, so the pause adds about half an hour; for a few changed slides it
   adds little, and `--pause 1` gives a quick run when the service is not
   throttling you.

   If even the first slide fails, Edge is refusing this machine for the time
   being (it does this after heavy use), and the run stops after two slides
   rather than retrying for hours. Check with
   `src/narration/.venv/bin/edge-tts --text test --write-media /tmp/t.mp3`;
   update with `src/narration/.venv/bin/pip install -U edge-tts`; or wait an
   hour or two. `--tts piper` uses an offline neural voice instead (Piper,
   installed and downloaded on first use, about 60 MB): steadier, a little less
   natural, and a different voice, so use it for a whole deck rather than for
   a few missing slides.

   **Clicks.** The voice service is asked for at most 1200 characters at a
   time, and the answers are appended to make each slide's clip. The log shows
   where they meet (`[4 pieces, bytes joins at 0:08.2, 0:28.0, 0:42.8]`), and
   so does `audio/slide_NN_joins.txt`, so a click can be checked against the
   joins. If clicks fall on the joins, `--resynth 14 --max-chars 0` sends
   slide 14 as one request (edge-tts still splits text over about 4096 bytes
   on its own), and `--join fade` decodes the pieces, fades 10 ms either side
   of each join and encodes the clip once. If they fall elsewhere, try
   `--reuse-audio --declick`: it runs ffmpeg's `adeclick` filter on the page's
   copy of each clip, takes seconds, and needs no voice service.
   Some voices (en-CA-LiamNeural more than en-US-AndrewNeural) put a soft
   low thump before an occasional word, like a "p" too close to a
   microphone. `--reuse-audio --highpass 80` removes what lies below 80 Hz
   from the page's copy of each clip; on a test signal it cut such thumps to
   about a tenth and lowered voiced speech by about 12 %. Use 60 if the voice
   sounds thin, 100 if thumps remain.

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
