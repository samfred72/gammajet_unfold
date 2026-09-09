# gammajet_unfold

Gamma-jet xJ unfolding analysis for sPHENIX. RooUnfold-based unfolding pipeline (`src/`),
driver macros (`macros/`), result/diagnostic plotting (`drawing/`), and an in-situ JES
calibration study (`insitu/`).

## Physics Ground Rules

- **xJγ is NOT bounded by construction**, unlike dijet xJ. Never truncate or drop
  high-xJ bins as "unphysical" — that assumption is valid for dijet, not gamma-jet.
- Low-count bins may still be excluded from chi2, but only with an explicit, stated
  threshold — see `drawing/draw_covariance_chi2.C:60`, which drops the last 3 xJ bins
  per pT bin for low statistics and documents why in the comment above it. Follow that
  precedent (explicit constant + comment) rather than silently filtering.
- Asymmetric systematic uncertainties are combined with a **per-source, per-bin sign
  split**, not a symmetric envelope: each two-point (high/low) source contributes its
  own signed value to the "up" quadrature sum if positive, "down" if negative, in each
  bin independently. See `asymmetricSystematics` and the total-uncertainty loop in
  `drawing/draw_systematics.C:103-390`. Symmetric sources instead contribute their full
  squared magnitude to both totals. Reuse this pattern for any new systematic source.
- State which of the above assumptions you're relying on before writing new
  uncertainty-combination or bin-selection code — don't import dijet-style conventions
  by default.

## Plotting Conventions

- Never combine statistical and systematic uncertainty into one band. Statistical stays
  on the data points' own error bars; systematic is drawn as a separate shape (currently
  a `TGraphAsymmErrors`-based box via `gSystBox`) — see `drawing/draw_final_result.C:34,
  132-158`. Follow this split for any new final-result plot.
- Output final results as PDF (written under `pdfs/`) — no HTML dashboards or
  interactive artifacts unless explicitly asked for.
- Sanitize histograms for NaN/Inf bin content before drawing. NaN poisons ROOT's axis
  auto-ranging and silently produces a blank canvas — this has happened before with
  closure plots.
- Set pad margins explicitly rather than relying on ROOT defaults, e.g.
  `p1->SetLeftMargin(.15)` / `SetBottomMargin(...)` as done throughout
  `drawing/draw_final_result.C:169-208`, so axis titles don't get clipped.
- The sPHENIX label, data/MC sample name, cut list, and jet radius are drawn together
  via `drawer::drawAll(samples, features, ...)` (`src/drawer.cc:30-42`), never as ad hoc
  `TLatex` calls in the macro itself. `samples` is a one-element vector naming the
  dataset (e.g. `"p+p Run24 Data"` or `"Pythia8 #gamma+jet MC"`); `features` lists the
  cuts/kinematics as separate strings, conventionally the pT bin range
  (`"%.0f GeV < p_{T}^{#gamma} < %.0f GeV"`) and jet radius (`"Jet R=%.1f"`, from
  `ana::JetRs[ir]`), optionally with the iteration count appended — see
  `drawing/draw_final_result.C:202-203` for the canonical call. Reuse this helper and
  format for any new plot rather than hand-rolling label text.

## Debugging Protocol

- Before proposing a root cause, read the actual code path first (grep for where a
  branch/node/histogram is filled and read) and report file:line evidence. Don't offer
  a sequence of speculative theories — verify the plumbing exists and executes before
  suspecting a subtle algorithmic bug.
- After any change to a `.cc`/`.h` in `src/`, rebuild via `src/make.sh` and report the
  actual compile output, not an assumption of success.

## SDCC Access

- The treemaking code (upstream of everything in this repo) lives on SDCC at
  `/sphenix/user/samfred/projects/gammajet/treemaking`.
- To reach it: `ssh sphnxuser04` (an interactive-node login, not `sphnx` — the shell
  alias `sphnx` in `~/.bashrc` expands to this, but Bash tool calls don't source
  `~/.bashrc`, so invoke `ssh sphnxuser04` directly rather than `ssh sphnx`). This
  matches the `Host sphnx*` block in `~/.ssh/config` (ProxyJump through `bnl`), so no
  further config is needed — just don't type the bare alias.

## Local RooUnfold Patch

- `/home/samson72/RooUnfold` (installed to `/home/samson72/root/lib/libRooUnfold.so`, the
  copy actually loaded at runtime - `/home/samson72/RooUnfold/libRooUnfold.so` is a build
  artifact, not itself on ROOT's library search path) has a local, uncommitted source
  patch in `src/RooUnfoldBayes.cxx`'s `getCovariance()`. Confirmed via a live backtrace
  (custom `SetErrorHandler` + `gSystem->StackTrace()`) during a gammajet_unfold session:
  when a response has fakes and `handleFakes=true` (this project's `unfold_utility.cc`
  always passes it), `setup()` increments `_nc` by one to add a synthetic truth-side
  "fakes" bookkeeping bin - but `getCovariance()`'s `_dosys`-gated "covariance due to
  unfolding matrix" term (only reached via `IncludeSystematics(...)`, which this project's
  `unfold_utility::unfoldOnce` does NOT call, so this bug is dormant unless you explicitly
  enable it) then reads `Eres(j, i)` for `i` up to the new (incremented) `_nc-1` - one
  column past the response matrix's own `Eresponse()`, which was never grown for the
  synthetic bin since there's no real response-matrix error to report there. ROOT's
  bounds-checked `TMatrixT::operator()` prints "Request column(N) outside matrix range of
  0 - N" and returns a fallback value instead of crashing, silently poisoning that
  covariance into NaN. Patched to treat the synthetic fakes bin's contribution as zero
  variance (`i < Eres.GetNcols() ? Eres(j,i) : 0.0`) instead of reading out of bounds -
  see the "LOCAL PATCH (gammajet_unfold...)" comment at that line in RooUnfoldBayes.cxx.
- A second, independent local patch: `RooUnfoldTH1Helpers.cxx`'s `h2meNorm<TH1,TH2>` -
  the function `RooUnfoldResponse::Eresponse()` calls to get the response matrix's error
  - was a byte-for-byte copy of `h2mNorm` (the CONTENT accessor) just above it: both
  called `GetBinContent()`, neither ever called `GetBinError()`. Confirmed against the
  `RooFitHist` specialization of the same function (`RooUnfoldFitHelpers.cxx`), which
  correctly uses `binError(...)` and is commented "sets Matrix to errors of bins" -
  proving the TH1/TH2 version was a copy-paste-and-forgot-to-change bug, not deliberate.
  Verified directly against a real response matrix: `Eresponse()` was returning the same
  values as `Mresponse()` (normalized migration probability, e.g. ~0.05) instead of any
  actual uncertainty, so `getCovariance()`'s response-matrix-statistics term
  (`IncludeSystematics(...)`-gated, same as above) was squaring migration probabilities
  and using that as the per-element variance - inflating that whole covariance term by
  (content/error)^2, a factor of ~2800x per element on the matrix checked. This is what
  made that term look absurdly large (e.g. ~219 analytic vs ~100 combined-toy at one
  bin) when first tried, well beyond what the toy bootstrap's own response-matrix-toy
  term ever found. Fixed by changing `GetBinContent()` to `GetBinError()` in that
  function - after the fix, the same bin's analytic error with systematics (~103) lines
  up with the toy bootstrap's combined estimate (~100) as expected.
- Neither patch is part of gammajet_unfold's own git history (RooUnfold is a separate
  checkout/build outside this repo) and there was already one other pre-existing local,
  uncommitted RooUnfold patch (a `verbose()>=1` gate on the "additional truth bin" print
  in `RooUnfoldBayes.cxx`'s `setup()`, presumably from an earlier session) - `git diff` in
  `/home/samson72/RooUnfold` is the actual record of all three. After any further edit
  there, rebuild with `cd /home/samson72/RooUnfold && make` (uses `ROOTSYS`/`root-config`,
  both already set up), then `cp libRooUnfold.so /home/samson72/root/lib/libRooUnfold.so`
  to actually deploy it - the build's own output `.so` is not on ROOT's load path by
  itself.

## Build & Run

- Rebuild the shared lib after source changes: `cd src && ./make.sh`.
- Driver scripts (`macros/runall_unfold.sh`, `runall_sys.sh`, `runall_puritymaker.sh`,
  `insitu/run_grid.sh`) are long-running — launch in background with output logged, and
  disable/skip the per-event progress bar for batch runs (it has previously been
  mistaken for the source of a slowdown that was actually just progress-bar overhead).
- **Do not launch production pipeline/driver scripts yourself.** Wire up the code and
  let the user run the pipeline/shell scripts once ready.
