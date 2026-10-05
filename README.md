# gammajet_unfold

Gamma-jet x_J unfolding analysis for sPHENIX: measures the unfolded x_J = p_T^jet /
p_T^gamma distribution in p+p Run24 Data, corrected for detector effects via
RooUnfold, with a full systematic-uncertainty budget. See `CLAUDE.md` for the physics
ground rules (xJ is not bounded, asymmetric-systematic combination convention, etc.)
and plotting conventions this repo follows - this file is about the *procedure*: what
to run, in what order, and why.

## Setup

Every macro and script finds the repository through `$GAMMAJET_UNFOLD`, and ROOT finds
`libgammajet_unfold.so` through `LD_LIBRARY_PATH`. Add to `~/.bashrc`, with the path of
your checkout:

```bash
export GAMMAJET_UNFOLD=/path/to/gammajet_unfold
export LD_LIBRARY_PATH=$GAMMAJET_UNFOLD/src${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
```

Then build with `src/make.sh`. ROOT and RooUnfold (with the local patches described in
`CLAUDE.md`) must be installed.

## Layout

| Directory | Contents |
|---|---|
| `src/` | The compiled analysis core (`libgammajet_unfold.so`): `unfolder.cc/.h` (response matrix + region A/B/C/D histogram production), `ana.h` (binning, cuts, physics constants - including the in-situ JES correction, see below), `unfold_utility.cc/.h` (purity correction), `insitu_utility.cc/.h` (shared helpers for `insitu/`), `treeuser.h`, `object`/`pho_object`/`jet_object`, `drawer.cc/.h`. |
| `macros/` | Driver scripts and the top-level `unfold.C`/`unfold_allsys.C`/`puritymaker.C` macros. `run_pipeline.sh` is the main entry point. |
| `drawing/` | Plotting/analysis macros that read `hists/` and produce `pdfs/`: purity-corrected spectra, systematics, the final result, closure/diagnostic checks. |
| `hists/` | `unfolder.cc`/`puritymaker.C` output (response matrices, region histograms, purity bootstrap) - one set per (trigger, sim, systag). |
| `pdfs/` | Final plots from `drawing/`. |
| `insitu/` | The in-situ jet-energy-scale (JES) calibration study - reads the `insitu_tree` ntuples `unfolder.cc` writes to `insitu/inputs/`, fits the Data-to-MC JES correction per jet radius, writes results to `insitu/output/` and `insitu/pdfs/`. **Its output feeds back into `src/ana.h` - see below.** |
| `insitu_closure/` | A closure test *of the in-situ fitting method itself* (not part of the physics result) - see its own section below. |
| `reweight/` | The Data/MC vz + cluster-pT reweighting used throughout `unfolder.cc`. |
| `latex/` | Note/paper source (gitignored). |

Raw input trees (`towerntup`, one file per trigger/sim, produced upstream by the
treemaking code) live in this repo's `trees/` directory (gitignored, not committed) -
see `src/treeuser.h` for the exact filename convention. The treemaking code itself
lives on SDCC (`ssh sphnxuser04`, see `CLAUDE.md`'s SDCC Access section) at
`/sphenix/user/samfred/projects/gammajet/treemaking`, now its own git repo
(`samfred72/gammajet_treemaking` on GitHub).

## The pipeline

`macros/run_pipeline.sh` is the reference driver - read its own header comment first
(`macros/run_pipeline.sh:1-33`), it has an explicit "which stage to rerun, by what you
changed" table. Summary of the four stages it runs, in order:

1. **Compile** (`src/make.sh`) - builds `libgammajet_unfold.so` from everything under
   `src/`. Needed whenever anything in `src/` changed, ana.h included.
2. **Unfold** (`unfolder.cc` via `runall_unfold.sh` / `runall_unfold_allsys.sh`) -
   reads the raw ntuples once per trigger, fills every systag's response matrix and
   region A/B/C/D histograms in one pass, and (as a side effect) writes the
   `insitu_tree` ntuples that `insitu/` consumes. **This is where
   `ana::jesNominal[ir]` (src/ana.h) is applied to Data's reconstructed jet p_T** -
   see `fill_matrix()`'s `jesCorrectionArr` in `src/unfolder.cc`.
3. **Purity** (`puritymaker.C` via `run_puritymaker.sh` / `runall_puritymaker.sh`) -
   reads stage 2's Data region histograms, writes the purity bootstrap
   `hists/purity_<systag>.root` that `ana::getPurity`/`getPurityC` read by index.
   Must be rerun any time stage 2 reran for Data, or `ana::ptBins` changed.
4. **Draw** (`drawing/*.C`) - purity-corrected spectra, systematics combination
   (`draw_systematics.C`), and the final result (`draw_final_result.C`): unfolded
   `(1/N)dN/dx_J` per p_T bin with statistical and systematic uncertainty shown
   separately (never combined into one band - see `CLAUDE.md`).

```
bash macros/run_pipeline.sh           # full: nominal + every systematic variation
bash macros/run_pipeline.sh nominal   # fast path for iterating on a code change
```

## In-situ JES calibration - an important loop-back, not a one-shot step

The in-situ study (`insitu/`) measures the Data-to-MC jet-energy-scale gap by
comparing the mean/shape of x_J between Data and MC around the photon-jet p_T balance
point, per jet radius, and writes the result as `ana::jesNominal[nJetR]` (plus its
stat+syst uncertainty, `jesTotalErrLow/High`) in `src/ana.h`. But **stage 2 above
already applies whatever `ana::jesNominal` currently holds** to Data's reconstructed
jet p_T before anything else happens - the response matrix, the purity histograms, the
region A/C events the in-situ study itself reads all depend on it.

That means the in-situ calibration and the main unfolding pipeline are coupled in a
loop, not a strict A-then-B order:

1. Run the full pipeline once (stages 1-4) with whatever `jesNominal` is currently in
   `src/ana.h` (a placeholder or a previous iteration's result) - this is what
   populates `insitu/inputs/*.root`.
2. Run the in-situ scan (from within `insitu/`, since its driver scripts/macros refer
   to each other by bare filename): `cd insitu && bash run_grid.sh` (or `cd insitu &&
   root -b -l -q 'grid_insitu.C("nominal")'` for just the nominal systag) - writes the
   purity-corrected best-fit p_a per radius to `insitu/output/` and `insitu/pdfs/`.
3. **Copy the new `jesNominal`/`jesTotalErrLow`/`jesTotalErrHigh` values into
   `src/ana.h`** (`cd insitu && root -b -l -q 'draw_jes_summary.C()'` plots/prints the
   per-radius summary across systags first, if you want a sanity check).
4. **Rerun the full pipeline (stages 1-4) again.** This is not optional and not just
   stage 2 - changing `ana::jesNominal` is an `ana.h` change, which per
   `run_pipeline.sh`'s own header comment means everything downstream (response
   matrix, purity, every drawing macro) is stale until stage 1-4 all rerun.
5. If you want to confirm the JES-fitting *method* itself is sound (not just rerun
   with new numbers), see `insitu_closure/` below - independent of this loop.

Skipping step 4 is the most likely way to end up with a "final" result that quietly
mixes an old, uncalibrated Data JES correction into the response matrix/purity/plots
while `insitu/`'s own plots show the new, correct calibration - i.e. two inconsistent
numbers for the same thing, split across `hists/` and `insitu/`.

## insitu_closure/ - validating the fitting method itself

`insitu_closure/` is a closure test of the in-situ *fitting procedure*
(mean-x_J-matching and shape-chi2 variants), not a step in the production pipeline
above and not something that feeds back into `ana::jesNominal`. It splits a stitched
Photon5+10+20 pythia MC sample into two random halves, injects a known, artificial
jet-energy-scale factor into one half ("data"), and checks whether the same
mean-x_J/shape-chi2 scan `insitu/grid_insitu.C` runs on real Data recovers that known
injected value from the other, unscaled half ("sim"). Run via:

```
cd insitu_closure
./run_closure.sh          # random injected scale in [0.9,1.1)
./run_closure.sh 0.95     # or an exact scale you choose
```

Useful when changing the in-situ scan's own fit logic (e.g. the low-xJ floor, the
chi2 error convention) to check the method still recovers a known answer, independent
of whatever the real Data/MC gap happens to be this production cycle.

## Systematics and the final result

`draw_systematics.C` re-unfolds Data's purity-corrected spectrum through each
systematic variation's own response + purity (one per `ana::systags` entry, see
`src/unfolder.h`'s constructor comment for what each one changes), including
`jes_high`/`jes_low` (Data's JES shifted by the in-situ study's own stat+syst
uncertainty on `jesNominal` - not a separate placeholder). Asymmetric (two-point
high/low) sources are combined with a per-source, per-bin sign split; symmetric
sources contribute their full magnitude to both totals - see `CLAUDE.md` and
`drawing/draw_systematics.C:103-390`. `draw_final_result.C` then plots the nominal
unfolded result with statistical and systematic uncertainty shown as separate bands.
