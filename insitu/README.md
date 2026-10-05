# insitu

In-situ jet energy scale (JES) calibration: Data jets are scaled so that Data matches the
MC reference.

There are two modes. Both run through `grid_insitu.C` and `run_grid.sh`.

| Mode | Fits | Uses | Output | Feeds `ana.h`? |
|---|---|---|---|---|
| `gammajet` (default) | constant scale `pa` | γ+jet only | `output/grid_insitu_<systag>.root`, `pdfs/grid_insitu_<systag>.pdf` | **yes**, via `draw_jes_summary.C` |
| `combined` | linear `f(pT) = pa + pb·pT` | γ+jet + multijet balance | `output/grid_insitu_combined_<systag>.root`, `pdfs/grid_insitu_combined_<systag>.pdf` | **no**, cross-check only |

## Running

```bash
./run_grid.sh                                        # gammajet, all systags, then draw_jes_summary.C (rewrites src/ana.h)
./run_grid.sh --systag nominal                       # gammajet, one systag
./run_grid.sh --mode combined                        # gamma+jet + multijet, all systags; ana.h untouched
./run_grid.sh --mode combined --systag nominal
./run_grid_nominal.sh [--mode combined]              # rebuild + nominal unfold + purity + scan
root -b -l -q 'grid_insitu.C("nominal","combined",2)'  # one radius (index into ana::JetRs) for a quick test
```

Run the gammajet mode first if you want the combined plots to overlay its constant scale.

## Inputs

**γ+jet (both modes):**
- `inputs/<trigger>[_<sim>]_<systag>_insitu.root` are the `insitutree` files written by
  `unfolder.cc`.
- Purity comes from `hists/purity_<systag>.root`.

**Multijet (combined mode only):** local copies of the multiJet analysis trees go in `../trees/`.

| File | Source on SDCC/tg |
|---|---|
| `multijet_Data.root` | `/sphenix/tg/tg01/jets/samfred/multiJet_full_hadded/multijet_Data.root` (`multiJet/FunforAll/hadd_data.sh`) |
| `multijet_pythia_Jet{8,12,20,30,50}.root` | `/sphenix/tg/tg01/jets/samfred/multiJet_hadded/` (`multiJet/FunforAll/condor_hadd.job`) |

A missing MC sample is skipped with a warning. If there is no multijet Data or MC at all,
that radius is skipped with an error.

## Combined mode: what it does

- **Multijet balance:** B = pT,lead / |p⃗T,sub + p⃗T,subsub|, in leading-jet pT bins
  20, 25, 30, 35, 40, 50, 60 GeV. Each jet is divided by f at its own pT. A constant scale
  cancels in B, so multijet constrains the slope `pb` and γ+jet fixes `pa`. This is the same
  method as the old `~/sphnx/gammajet/macros/grid_insitu.C`.
- **Selection:** re-applied on the pT actually used: Data `jet_pt_calib`, MC
  `jet_pt_smear_truth` (`_high`/`_low` for JERhigh/JERlow; no other systag changes the
  multijet side).
  - at least 3 jets, all with |η| < 1.1 − R;
  - sub and subsub pT ≥ 7 GeV;
  - Δφ(lead, sub) ≥ 3π/4 and Δφ(lead, subsub) ≥ π/2;
  - 0.4 ≤ B < 2.65.

  These thresholds match `DijetTreeMaker`'s skim. The trees already have the skim's
  |vz| < 60 cm and the Data jet timing cuts applied.
- **MC reference:** Pythia8 Jet8–50. Each sample is cut to its pT-hat slice
  (`treeuser::truthSliceLow/High`, leading truth jet at the same R) and weighted by the
  `drawer.h` cross-section scale. The multiJet MC trees carry no weights and no
  generated-event counts, so this follows the repo's convention. It is exact only if the
  samples' event counts are comparable.
- **Fit:**
  - The χ² is the γ+jet χ² (region A, or purity-corrected) plus the multijet χ². Both use
    the same per-bin mean-ratio formula as the gammajet mode.
  - The grid is `pa` 0.80–1.00 in 200 steps × `pb` ±0.005 /GeV in 100 steps.
  - The 1σ region is Δχ² < 2.30. Its `pa`/`pb` extent and its f(pT) envelope are saved,
    with a warning if the region touches the grid edge.
- **Pages per radius:**
  - γ+jet region A + multijet;
  - purity-corrected + multijet;
  - multijet balance (Data vs MC, raw and corrected);
  - f(pT) bands with the gammajet-only constant overlaid;
  - χ² map.

## Other macros

These are constant-scale and gammajet only.

| Macro | Purpose |
|---|---|
| `grid_insitu_shapechi2.C` | xJ shape χ² instead of mean |
| `grid_insitu_unfolded.C`, `grid_insitu_unfolded_shapechi2.C` | compare after unfolding |
| `grid_insitu_jet12.C`, `grid_insitu_jet12_shapechi2.C` | Jet12_long dijet-MC reference |
| `draw_grid_chi2.C` | χ² comparison across methods |
| `draw_jes_summary.C` | collects the gammajet-mode `pa` for all systags and rewrites `src/ana.h` (`jesNominal`, `jesBySystag`) |
| `draw_jes_variations.C` | `pa` under each variation |
| `draw_insitu_*.C` | xJ comparison plots |
| `debug_shapechi2_*.C` | diagnostics |

`archive_expanded_ptbins/` holds the archived outputs from the temporary 9-bin pT study.
