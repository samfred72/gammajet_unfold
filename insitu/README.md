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

**Multijet (combined mode only):** `../multijet/multijet_analysis_pythia.root`, written by
`../multijet/analysis` (see `../multijet/README.md`). Run that first: it reads the multiJet
trees in `../trees/` and does the multijet event selection and MC weighting. This macro
only reads its per-radius trees:
- Data: `ttree_data_r<10R>`
- MC: `ttree_Jet{8,12,20,30}_r<10R>_{RECO,HIGH,LOW}`

A missing MC tree is skipped with a warning. If there are no multijet Data or MC events at
all, that radius is skipped with an error.

## Combined mode: what it does

- **Multijet balance:** B = pT,lead / |p⃗T,sub + p⃗T,subsub|, in `analysis.cc`'s
  leading-jet pT bins (20, 25, 30, 35, 40, 50, 60, 70 GeV), for 0.4 ≤ B < 2.65. Each jet is
  divided by f at its own pT. A constant scale cancels in B, so multijet constrains the
  slope `pb` and γ+jet fixes `pa`. This is the same method as the old
  `~/sphnx/gammajet/macros/grid_insitu.C`, which read the same kind of tree.
- **Selection and MC weights:** taken entirely from `analysis.cc`. They are not
  re-implemented here, so the two cannot drift apart.
- **MC jet pT:** `analysis.cc`'s reco-based smearing `jet_pt_smear_reco`. JERhigh/JERlow
  use its HIGH/LOW variants; no other systag changes the multijet side. This differs from
  the γ+jet side, which uses `jet_pt_smear_truth` as `unfolder.cc` does.
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
| `grid_insitu.C(systag, "gammajet", -1, "shape")` | xJ shape χ² instead of mean (`grid_insitu_shapechi2_<systag>.*`) |
| `grid_insitu_unfolded.C(systag, na, "mean"\|"shape")` | compare after unfolding |
| `grid_insitu_jet12.C(systag, "mean"\|"shape")` | Jet12_long dijet-MC reference |
| `draw_grid_chi2.C` | χ² comparison across methods |
| `draw_jes_summary.C` | collects the gammajet-mode `pa` for all systags and rewrites `src/ana.h` (`jesNominal`, `jesBySystag`) |
| `draw_jes_variations.C` | `pa` under each variation |
| `draw_insitu_*.C` | xJ comparison plots |
| `debug_shapechi2_*.C` | diagnostics |

`archive_expanded_ptbins/` holds the archived outputs from the temporary 9-bin pT study.
