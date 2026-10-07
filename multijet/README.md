# multijet

Multijet balance analysis, moved here from SDCC `multiJet_legacy/` (Oct 2026). Tree
making stays in the SDCC `multiJet` repo (`DijetTreeMaker`).

The balance is x_j = pT,lead / |p⃗T,sub + p⃗T,subsub|, binned in leading-jet pT
(20, 25, 30, 35, 50 GeV), for all seven radii.

Selection: leading jet ≥ 20 GeV, jets 2 and 3 ≥ 7 GeV each, recoil |p⃗T,2 + p⃗T,3| ≥ 14 GeV,
|η| < 1.1 − R for the three jets, Δφ₁₂ > 3π/4, Δφ₁₃ > π/2. MC jets get the JER smearing from
`aux/jer_smear_templates.root`. In MC events with exactly two truth jets ≥ 7 GeV, jets 2 and 3 are
ranked and added unsmeared and the sum is smeared once, at the pT of the truth jet matched to the recoil
(the 7 GeV cuts use the values before that smearing). Every other event smears each jet on its own.
The leading reco jet is capped at 35 GeV in Jet12 and 50 GeV in Jet20.

## Inputs

The local multiJet trees in `../trees/`, the same files the in-situ combined mode reads
(see `../insitu/README.md` for where to copy them from):
- `multijet_Data.root`
- `multijet_<sim>_Jet{5,8,12,20,30}.root`, with sim = pythia or herwig (herwig has no Jet8;
  Jet12, Jet20 and Jet30 are required, the others are used if present)

MC samples are stitched on the leading truth jet with the MDC2 full-efficiency thresholds: each
sample covers its threshold up to the next used sample's, and the lowest used sample starts at 0.
Jet8 is not in the MDC2 threshold table, so Jet5 and Jet8 cannot be combined until its
thresholds are filled in (`truthThreshold` in `analysis.cc`); use `--no-jet8` with Jet5.

## Programs

| Program | What it does |
|---|---|
| `analysis.cc` (`./make.sh`, then `./analysis <0=pythia\|1=herwig> [tree dir] [--no-jet5] [--no-jet8] [--truth-smear] [--no-smear] [--tight]`) | Selects multijet events, fills the x_j histograms per radius, pT bin and JER variation (RECO/HIGH/LOW = `jet_pt_smear_{,high_,low_}reco`), and writes `multijet_analysis_<sim>.root`. That file also holds per-radius trees (`leadingPT`, `SLPT`, `SSLPT`, φ, η, `PT23`, `weight`) and `pdfs/unmatched_event_display_<sim>.pdf`. The Data trees keep events down to a JES of 0.9 (cuts × 0.9) for the in-situ scans; the histograms use the nominal cuts. `--tight` (leading > 25, recoil > 9 GeV) writes `multijet_analysis_<sim>_tight.root`. |
| `makeratio.C` | From `multijet_analysis_<sim>.root`, fits the Data/MC leading-pT ratio (expo, per JER variation) and takes the z-vertex ratio, then writes `aux/ratio<R>_<sim>.root`. `analysis.cc` uses these to reweight MC on its **next** pass. |
| `draw_xj.C(radius, sys, sim)` | Mean x_j vs leading pT, Data vs MC → `pdfs/hxj_means_r<R>_<sys>_<sim>.pdf`, and the distributions per pT bin → `pdfs/hxj_dists_…`. `draw_xj_all(sim)` (after `.L draw_xj.C`) writes every radius and JER variation into `pdfs/hxj_means_<sim>.pdf` and `pdfs/hxj_dists_<sim>.pdf`. |
| `draw_xj_samples.C(file, pdf)` | x_j per MC sample per pT bin, with the total MC and Data, every radius and JER variation in one PDF |
| `smearcheck.C` | Calibrated vs smeared pT checks on Pythia Jet12 → `checkhists.root` |

`aux/ratio*.root` are the reweighting files copied from `multiJet_legacy/aux` (May 26). They
were derived from the old skimmed trees, so after the first pass on the new trees run
`makeratio.C` and then `analysis` again. The reweighting is iterative: the fits are built
from the previous pass's spectra.

## Changes from the legacy copy

- **Input:** the local `../trees/multijet_*` files, with an optional directory argument.
  The legacy copy read hard-coded SDCC `multiJet_skimmed` / `multiJet_06152026` paths.
- **Branch names:** the current `DijetTreeMaker` names.
  - `jet_pt_smearRECO/HIGH/LOW_<R>` → `jet_pt_smear_reco/high_reco/low_reco_<R>`
  - `jet_pt_smear_<R>` → `jet_pt_smear_truth_<R>`
- **`draw_xj.C`:** reads `analysis.cc`'s current output. It was written for an older
  `hists/hists_<R>_pythia.root` histogram layout.
- **Build:** `-std=c++17` (local ROOT 6.28), via `make.sh`.

The selection, weights and histogram content are unchanged.

## Used by `../insitu` (combined mode)

`insitu/grid_insitu.C("<systag>","combined")` reads this analysis's output,
`multijet_analysis_pythia.root` (per-radius trees `ttree_data_r<R>` and
`ttree_Jet*_r<R>_{RECO,HIGH,LOW}` with their `weight`). So the multijet selection and MC
weighting are defined only here.

The order is:
1. `./analysis 0`
2. `makeratio.C`, then `./analysis 0` again, to update the reweighting
3. `../insitu/run_grid.sh --mode combined`

## Running on Alpine (Slurm)

`bash slurm/multijet/submit.sh [--passes N] [--sims "pythia herwig"] [--no-draw] [--no-jet8] [--tight]` submits
build → `analysis` (one array task per sim) → `makeratio.C` → `analysis` … → `draw_xj.C`
(one array task per radius) as a dependency chain. The default is 2 passes, i.e. one reweighting
update, as in the order above. Both reweighting inputs (leading pT and z vertex) are filled after the
full selection with the cross-section weights only, so the second pass is final. `--tight` runs one
pass with the current fits and no `makeratio`.
Logs, each pass's `multijet_analysis_<sim>.root`, and the `aux/` fits after each `makeratio`
go to `logs/slurm/multijet_<timestamp>/`.
