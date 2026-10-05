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

## Known differences with `../insitu` (combined mode)

1. **Event selection.** `analysis.cc` additionally requires:
   - Data trigger bit 22;
   - sub and subsub pT < 30 GeV;
   - its own timing cut (|t_lead| < 6, |t_lead − Δt| < 3);
   - MC reweighting in leading pT and z-vertex.

   The in-situ combined mode re-applies only the `DijetTreeMaker` skim cuts plus |η| < 1.1 − R.
2. **MC normalization.**
   - `analysis.cc` uses its own per-sample weights, normalized to Jet30: Jet12 = 1.4903e6,
     Jet20 = 6.2623e4 (pythia).
   - `insitu/grid_insitu.C` uses `drawer.h`'s scale map: Jet12 = 3.997e6, Jet20 = 6.218e4,
     and it also uses Jet50.
   - They also differ in the top pT-hat slice: Jet30 runs to 100 GeV here, 50 GeV in `treeuser`.

These should be reconciled before the two are compared.
