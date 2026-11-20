# Narrated slide pages

A slide deck in `docs/interactive/` can be narrated: its speaker notes are read
by a text-to-speech voice into one HTML file with the slides, an audio player
per slide and the notes as a transcript, which plays online and also works when
downloaded. The tools are the [`narrate`](https://github.com/billpage/narrate)
package (they were in this folder until November 2026). This folder keeps only
the project's pronunciations:

- [`wpmw_terms.json`](wpmw_terms.json) — how to say this project's notation
  and words: `{"written": "spoken"}`, applied before synthesis.

The project's settings are in [`narration.toml`](../../narration.toml) at the
repository root: the Canadian Edge voice `en-CA-LiamNeural`, a 30 s pause
between requests to the voice service, the terms above, 32 kbit/s audio, and
where renders and the `output` worktree are. An option on the command line
overrides them (`--voice en-US-AndrewNeural`, `--pause 1`); `narrate config`
shows what is in force.

**Install** (once; again with `--force` and a new tag to move to it):

```bash
sudo apt install pipx ffmpeg
pipx install "narrate[pptx] @ git+https://github.com/billpage/narrate@v1.0.0"
```

**What is kept where.** `main` keeps the inputs: a deck's source in
`talks/<name>/` (`slides.html` and one `notes/slide_NN.md` per slide), the
slides page built from it in `docs/interactive/<name>.html`, and the terms
file. The rendered page, with its audio embedded (about 15 to 20 MB for an
hour of speech), goes to `narrated/` on the `output` branch and is served at
`https://billpage.github.io/wpmw/narrated/<name>_narrated.html` by the
[`pages`](../../.github/workflows/pages.yml) workflow. The clips stay on your
machine (in `~/wpmw-narration/<name>/audio/`), so a page can be rebuilt
without synthesising again.

## Workflow

The deck [worlds, particles and world-particles](../../talks/worlds_and_particles)
is the example.

1. **Edit the source and build the slides page.** Slides are in
   `talks/worlds_and_particles/slides.html`, the speaker notes in
   `talks/worlds_and_particles/notes/slide_NN.md` (plain text, read aloud as
   written). Build:

   ```bash
   narrate build-deck talks/worlds_and_particles --foot "WPMW · worlds and particles" \
       --out docs/interactive/worlds_and_particles.html
   ```

   This adds the footer `WPMW · worlds and particles · k / N` to every slide
   without `data-nofoot`, and narrate's page frame: the notes at the right of
   the slide (below it in a narrow window), N to hide and show them, arrow
   keys and the bar to move. Commit the source and the built page together.

2. **Check the text** (no audio, a few seconds):

   ```bash
   narrate render docs/interactive/worlds_and_particles.html --preview
   ```

   Read `~/wpmw-narration/worlds_and_particles/worlds_and_particles_spoken.txt`.
   Any slide marked `NOT SPOKEN WELL` has characters the voice will stumble on.
   Pronunciations come from narrate's common terms (13 physics symbols such as
   `ψ₀` and `ℏ`), then `wpmw_terms.json`, then a `pronunciations.json` beside
   the slides page if there is one; a later file wins on the same key.

3. **Render** (Edge voice, needs internet; with the 30 s pause, about half an
   hour for an hour of speech, almost all of it waiting):

   ```bash
   narrate render docs/interactive/worlds_and_particles.html
   ```

   Open `~/wpmw-narration/worlds_and_particles/worlds_and_particles_narrated.html`
   in a browser and listen.

   A clip is kept under a hash of its spoken text, voice, rate and settings
   (`audio/clips/`), so only slides whose text changed are synthesised again,
   even after slides are inserted or moved. `--resynth` makes every clip
   again; `--resynth 14` or `--resynth 5,14` makes only those;
   `--reuse-audio` rebuilds the page without contacting the voice service.

   Microsoft's voice service sometimes stalls or refuses a request it accepts
   a minute later. `voiceover` asks for at most about 1200 characters at a
   time, abandons a request after 60 s (`--timeout`), and retries a failed one
   three times, 5, 20 and 60 s apart. If a slide still fails, the run ends with
   `MISSING AUDIO for slide(s) ...` and exit status 3, and `narrate publish`
   refuses the page; run the same command again. The service also refuses
   many first tries when requests come quickly, hence the 30 s pause (time
   between requests, not silence in the audio); `--pause 1` gives a quick run
   when it is not throttling you.

   If even the first slide fails, Edge is refusing this machine for the time
   being (it does this after heavy use), and the run stops after two slides.
   Check with `~/.local/share/pipx/venvs/narrate/bin/edge-tts --text test
   --write-media /tmp/t.mp3` (`pipx list` shows the environment's folder if
   it is elsewhere); update with `pipx runpip narrate install -U
   edge-tts`; or wait an hour or two. `--tts piper` uses an offline neural
   voice instead (`pipx inject narrate piper-tts` once; the voice, about 60 MB,
   downloads on first use): steadier, a little less natural, and a different
   voice, so use it for a whole deck rather than for a few missing slides.

   **Clicks.** The answers to the requests are appended to make each slide's
   clip; the log shows where they meet (`[4 pieces, bytes joins at 0:08.2,
   0:28.0, 0:42.8]`), and so does `audio/slide_NN_joins.txt`. If clicks fall on
   the joins, `--resynth 14 --max-chars 0` sends slide 14 as one request, and
   `--join fade` fades 10 ms either side of each join and encodes the clip
   once. If they fall elsewhere, `--reuse-audio --declick` runs ffmpeg's
   `adeclick` on the page's copy of each clip. Some voices (en-CA-LiamNeural
   more than en-US-AndrewNeural) put a soft low thump before an occasional
   word; `--reuse-audio --highpass 80` removes what lies below 80 Hz (60 if the
   voice sounds thin, 100 if thumps remain).

4. **Fix what sounds wrong** in `wpmw_terms.json` and render again. A word
   whose pronunciation depends on context (*winds*, *live*, *read*) needs the
   surrounding words in the key, because keys are matched as plain text,
   even inside other words: `"clock winds at": "clock wynds at"`, never
   `"wind": ...`. To choose a spelling by ear, `narrate test-pronunciation`
   makes a page of players (see `narrate --help`).

5. **Publish**:

   ```bash
   narrate publish ~/wpmw-narration/worlds_and_particles/worlds_and_particles_narrated.html
   ```

   It refuses a page with missing audio or with linked (`--link`) audio,
   commits the page to `narrated/` on `output` in `~/wpmw-output`, pushes, and
   starts the Pages workflow (with the GitHub CLI `gh` if it is logged in;
   otherwise it prints the link to the workflow's *Run workflow* button).

## The page

Left and right arrow keys change slide; Space or K plays and pauses; J and L
go back and forward 10 s; T shows or hides the transcript; C turns on
continuous play, which moves to the next slide at the end of each clip.
Drag the divider between slide and transcript to resize it, or focus it and
use the arrow keys; double-click resets it. The page remembers speed,
divider, transcript and continuous play in the browser. *Save* downloads the
page as one file that plays offline.

`voiceover --help` lists every option of the narrator itself, which `narrate
render` passes on. A `.pptx` source also needs LibreOffice and `pdftoppm`
(poppler-utils); `--tts silent` makes silent clips of the right length, for
checking layout without a voice.
