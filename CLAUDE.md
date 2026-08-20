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

## Build & Run

- Rebuild the shared lib after source changes: `cd src && ./make.sh`.
- Driver scripts (`macros/runall_unfold.sh`, `runall_sys.sh`, `runall_puritymaker.sh`,
  `insitu/run_grid.sh`) are long-running — launch in background with output logged, and
  disable/skip the per-event progress bar for batch runs (it has previously been
  mistaken for the source of a slowdown that was actually just progress-bar overhead).
- **Do not launch production pipeline/driver scripts yourself.** Wire up the code and
  let the user run the pipeline/shell scripts once ready.
