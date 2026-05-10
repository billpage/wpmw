# Empirical test: does GitHub strip backslash escapes inside fenced math blocks?

**Status:** experimental, not for permanent inclusion. Push to a throwaway
branch only.

GitHub's markdown preprocessor strips backslashes from CommonMark-escapable
punctuation inside `$$...$$` math blocks (community discussion
[#17143](https://github.com/orgs/community/discussions/17143)). Whether the
same happens inside fenced math blocks (the alternative GitHub display-math
syntax) is **undocumented** — the spec implies fenced-code content shouldn't
go through inline processing, but GitHub's actual implementation has its
own quirks.

Each test below pairs the same pattern in both syntactic forms. View the
rendered file on GitHub and compare. If the two forms render identically,
the strip applies to both and our style-guide caveat stays. If the fenced
form renders correctly where the `$$` form is broken, fenced math is the
clean fix.

---

## Test 1 — sized brace `\bigl\{...\bigr\}` (binary: error or no error)

This is the test that gave the "Missing or unrecognized delimiter for
\bigl" MathJax error originally. Pure yes/no signal.

`$$` form (known broken):

$$\Omega = \bigl\{(x, y)\bigr\}$$

Fenced form (under test):

```math
\Omega = \bigl\{(x, y)\bigr\}
```

What to look for:
- If `\bigl\{` is **broken** in the fenced form too: a red MathJax error
  identical to the one above the fence.
- If `\bigl\{` is **exempt** in the fenced form: a cleanly rendered
  set-builder expression with curly braces around `(x, y)`.

---

## Test 2 — thick space `\;`

`$$` form (known broken — renders as literal semicolons):

$$a \; = \; b$$

Fenced form (under test):

```math
a \; = \; b
```

What to look for:
- If broken: visible semicolons in `a ; = ; b`.
- If exempt: `a = b` with a noticeable extra gap on each side of `=`.

---

## Test 3 — thin space `\,`

`$$` form (known broken — renders as literal commas):

$$x \, y \, z$$

Fenced form (under test):

```math
x \, y \, z
```

What to look for:
- If broken: visible commas in `x , y , z`.
- If exempt: three letters separated by narrow gaps (no commas).

---

## Test 4 — control: `\thinspace` (no backslash-escape characters)

The macro name `\thinspace` contains no ASCII punctuation, so CommonMark
has nothing to strip. Both forms should render identically — a sanity
check that the rendering pipeline is alive.

`$$` form:

$$p \thinspace q \thinspace r$$

Fenced form:

```math
p \thinspace q \thinspace r
```

If either of these shows a literal `\thinspace` in the rendered output,
something else is wrong (different MathJax config between the two
syntaxes, or `\thinspace` not in the loaded package set) and the other
results need a second look.

---

## How to record the result

After viewing on GitHub, one of these is true:

- [ ] **Tests 1–3 render identically broken in both forms.** Fenced math
      is *not* exempt; the strip applies inside fenced math too.
      Style-guide caveat stands; same escape rules everywhere.
- [ ] **Tests 1–3 render correctly in the fenced form, broken in the `$$`
      form.** Fenced math *is* exempt; recommend it as the clean
      alternative whenever backslash-escapes are needed.
- [ ] **Mixed results** — different tests behave differently. Worth a
      closer look at which patterns are affected and which aren't.

Reply with which box applies and I'll prepare a follow-up patch
adjusting the style guide and (if exempt) recommending the fenced form.
