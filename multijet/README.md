# multijet

Multijet balance analysis, moved here from SDCC `multiJet_legacy/` (Oct 2026). Tree
making stays in the SDCC `multiJet` repo (`DijetTreeMaker`).

The balance is x_j = pT,lead / |p⃗T,sub + p⃗T,subsub|, binned in leading-jet pT
(20, 25, 30, 35, 40, 50, 60, 70 GeV), for all seven radii.

## Inputs

The local multiJet trees in `../trees/`, the same files the in-situ combined mode reads
(see `../insitu/README.md` for where to copy them from):
- `multijet_Data.root`
- `multijet_<sim>_Jet{8,12,20,30}.root`, with sim = pythia or herwig

## Programs

| Program | What it does |
|---|---|
| `analysis.cc` (`./make.sh`, then `./analysis <0=pythia\|1=herwig> [tree dir]`) | Selects multijet events, fills the x_j histograms per radius, pT bin and JER variation (RECO/HIGH/LOW = `jet_pt_smear_{,high_,low_}reco`), and writes `multijet_analysis_<sim>.root`. That file also holds per-radius trees (`leadingPT`, `SLPT`, `SSLPT`, φ, η, `PT23`, `weight`) and `pdfs/unmatched_event_display_<sim>.pdf`. |
| `makeratio.C` | From `multijet_analysis_<sim>.root`, fits the Data/MC leading-pT ratio (expo, per JER variation) and takes the z-vertex ratio, then writes `aux/ratio<R>_<sim>.root`. `analysis.cc` uses these to reweight MC on its **next** pass. |
| `draw_xj.C(radius, sys, sim)` | Mean x_j vs leading pT, Data vs MC → `pdfs/hxj_means_r<R>_<sys>_<sim>.pdf` |
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

`bash slurm/multijet/submit.sh [--passes N] [--sims "pythia herwig"] [--no-draw]` submits
build → `analysis` (one array task per sim) → `makeratio.C` → `analysis` … → `draw_xj.C`
(one array task per radius) as a dependency chain. The default is 2 passes, i.e. one reweighting
update, as in the order above.
Logs, each pass's `multijet_analysis_<sim>.root`, and the `aux/` fits after each `makeratio`
go to `logs/slurm/multijet_<timestamp>/`.
