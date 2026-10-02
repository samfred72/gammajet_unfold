# Why is the data purity lower than the Jet12 (Pythia) purity at 15–20 GeV?

At 15–20 GeV, data P_A = 0.61 at every R, against Jet12_long P_A = 0.74 (ABCD method) and
0.77 (true, from truth tagging). At 20–25 and 25–35 GeV data and Jet12 agree within
errors. See `../insitu_scan_edge/README.md`, check 3.

The quick diagnostics here are run interpreted (`root -l -b -q macro.C`); each prints a
table.

| macro | what it checks |
|---|---|
| `j12comp.C` | truth content of Jet12_long vs Jet12 vs Photon10 |
| `abcdcomp.C` | ABCD populations (committed `hclusterpt_abcd`), data vs Jet12, plus leakage fractions |
| `isocomp.C` | iso4 of **tight**-BDT clusters: data vs truth-tagged γ+jet MC / Jet12, Jet12 background |
| `isocomp_bkg.C` | the same for **non-tight** (BDT 0.2–0.6) clusters, which are background-dominated |

## Findings (15–20 GeV unless noted)

1. **The Jet12 purity is a genuine Pythia prediction, not an artifact.** 3.3% of
   Jet12_long events have a truth photon at 15–20 GeV. 96% of those are isolated, and
   their truth ⟨jet/γ⟩ = 0.73 matches Photon10 (0.73). The plain Jet12 sample is
   identical. (The generator config itself couldn't be checked; SDCC `/sphenix/user` was
   stale.)
2. **The whole gap sits in one ABCD ratio.** Region A/C is the same in data and Jet12
   (1.05 vs 1.03 at R = 0.2; 1.08 vs 1.06 at R = 0.4). B/D, the tight-to-non-tight ratio
   among **non-isolated** (background-dominated) clusters, is **0.51 in data vs 0.34 in
   Jet12** (0.32 for Jet12 background only). Since P ≈ 1 − (B/D)(C/A), that ratio
   accounts for the difference. It narrows at 20–25 GeV (0.46 vs 0.40) and reverses at
   25–35 GeV (0.37 vs 0.50, low statistics). This matches the pT dependence of the
   purity gap.
3. **It is not extra signal leakage into the non-isolated region.**
   - At 25–35 GeV, tight data clusters have 8.3% with iso4 > 4 GeV. MC signal and
     background shapes mixed at the data purity (~0.8) give 8.9%, so MC's
     signal-leakage model holds there.
   - At 15–20 GeV the same mix needs P ≈ 0.63 to reproduce data's 11.4%, consistent with
     the ABCD 0.61.
   - Tight data clusters have a higher median iso (0.68 GeV) than MC signal (0.38 GeV).
     That is explained by the ~37% background content, not a shifted signal peak.

   So the data sample at 15–20 GeV really does seem to contain more background.
4. **Data fakes are somewhat more isolated than Pythia's.** After removing the ~15%
   signal (P_C) from the non-tight data sample, 55% of its background has iso4 < 2 GeV,
   against 50% for Jet12 background (> 4 GeV: 24% vs 27%). That is consistent with a
   lower calorimeter response in data (the in-situ JES is ~0.93) and/or dead or masked
   towers in the isolation cone. It moves the purity by only a few percent, not enough on
   its own.
5. **The main driver is that data fakes pass the tight BDT about 1.5× more often than
   Pythia fakes at 15–20 GeV** (finding 2). Candidates, not yet separated:
   - shower-shape mismodeling of low-pT merged π0/η clusters. The BDT
     (`model_base_v3E_single_tmva`) is trained on MC, and the effect shrinks with pT as
     π0 showers become more photon-like in both;
   - a different fake composition in data, e.g. η/π0 ratio, leading-meson momentum
     fraction, or hadrons (n̄, K0L) showering in the EMCal;
   - more fakes per prompt photon than LO Pythia predicts.

   All of these make the data **really** less pure, which ABCD measures correctly as
   long as isolation and BDT are uncorrelated for background. In Jet12 background that
   correlation (A·D/B·C) is 0.92 at 15–20 GeV, i.e. about an 8% violation. At P ≈ 0.6
   that shifts the purity estimate by roughly 0.03.

## Next checks (not done)

- BDT score distribution of non-isolated (iso4 > 4 GeV) clusters, data vs Jet12
  background, at 15–20 and 25–35 GeV, plus the individual shower-shape inputs. This
  would confirm finding 5 directly.
- A·D/(B·C)-type correlation in a data sideband (e.g. a still-looser BDT band) to check
  the ABCD independence assumption in data.
- Purity vs |η| (the reviewer's region-C η observation). Dead-tower or acceptance holes
  in the isolation cone would show up there.

---

# Issue 10: purity vs |η^γ| (`purity_vs_eta.C`)

**Method.** The production chain is run in bins of |η^γ| (0–0.35, 0.35–0.7, 0.7–1.1) and
integrated:
- data ABCD counts of paired clusters, rebuilt from `trees/gammajet_Data.root` with
  `unfolder.cc`'s data branch;
- truth-matched Photon5+10+20 leakage templates, using `unfolder.cc`'s MC branch,
  `treeuser.h` truth windows, `mcWeight` and `drawer` `scalemap` weights, each in the
  same |η| bin;
- the purity from `puritymaker.C::combine_hists` itself (`#include`d).

R = 0.2 and 0.4. Output is `pdfs/purity_vs_eta.pdf` and `.root`, and
`purity_vs_eta.log`.

**Validation.** The integrated set reproduces the committed `hclusterpt_abcd` data counts
**exactly** in every bin, the leakage c to ≤ 0.004, and `purity_nominal.root` to ±0.001
(R = 0.4: 0.621 / 0.784 / 0.809 vs 0.620 / 0.783 / 0.809).

**Results, R = 0.4, 15–20 GeV** (R = 0.2 is the same within ±0.02):

| \|η^γ\| | A | C/A | B/D | leakage c (b) | P_A | P_C |
|---|---|---|---|---|---|---|
| 0–0.35 | 1585 | 0.62 | 0.59 | 0.16 (0.022) | 0.716 ± 0.027 | 0.18 |
| 0.35–0.7 | 1136 | 0.93 | 0.48 | 0.22 (0.018) | 0.637 ± 0.039 | 0.15 |
| 0.7–1.1 | 813 | 1.53 | 0.33 | 0.34 (0.012) | 0.580 ± 0.065 | 0.13 |
| all | 3534 | 0.93 | 0.50 | 0.22 (0.018) | 0.621 ± 0.025 | 0.15 |

At 20–25 GeV the per-|η| purities are 0.77 / 0.82 / 0.90 (errors 0.05–0.07). At 25–35
GeV the errors are 0.05–0.39, so no trend can be read.

**Findings**

1. **The region populations depend strongly on |η|, in data and in the signal MC.**
   - C/A rises from 0.62 to 1.53 and B/D falls from 0.59 to 0.33 going forward.
   - Isolated clusters become more common at large |η|, for signal (b: 0.022 → 0.012)
     and background alike. That is the reviewer's point: the R = 0.4 isolation cone runs
     out of the EMCal acceptance for |η| > 0.7, so less energy is collected and clusters
     look more isolated.
   - Real photons also fail the tight BDT twice as often forward (c: 0.16 → 0.34),
     consistent with shower shapes changing with incidence angle and at the calorimeter
     edge.
2. **The purity falls with |η| at 15–20 GeV** (0.72 → 0.64 → 0.58, about 2σ between the
   outer bins). At 20–25 GeV it trends the other way within errors.
3. **The |η|-integrated ABCD underestimates the purity by about 0.03 in every bin and at
   both radii.** The region-A-weighted average of the per-|η| purities exceeds the
   integrated value by:
   - R = 0.2: +0.036 / +0.027 / +0.031 (15–20 / 20–25 / 25–35 GeV);
   - R = 0.4: +0.038 / +0.030 / +0.029.

   The shift is the same size and sign everywhere, including the noisy bins, so it
   comes from the method, not statistics. The reason: |η| is a common cause of both
   "isolated" (acceptance) and "fails tight" (shower shape). Integrated over |η|, that
   creates an isolation–BDT correlation among background clusters, which ABCD reads as
   extra background in A. Within an |η| slice the independence assumption holds better.
   The shift is comparable to the quoted purity uncertainty (±0.025 at 15–20 GeV).
4. **Size of the effect on the result.** From the reviewer's backup P_A ± 1σ rows
   (σ ≈ 0.025), a +0.03 shift in P_A moves the 15–20 GeV xJ shape by about −1.5% at low
   xJ and +5–10% above xJ ≈ 1.2. For the in-situ fit it slightly reduces the
   subtraction: the corrected ⟨xJ⟩ rises by ≈ 0.003–0.004, i.e. p_a by ≈ +0.5%.

**Options.**
- (a) Do the purity correction in |η| bins: purity per (pT, |η|), subtract per |η|, then
  sum. This is the most direct fix, and the macro already produces the inputs.
- (b) Keep the integrated purity but take the per-|η| difference (+0.03) as an
  additional one-sided purity systematic.
- (c) Restrict photons to |η| < 0.7 as in PPG12, which removes the forward bin where the
  cone leaks out. It costs about 23% of region A at 15–20 GeV, and the particle-level
  definition (|η^γ| < 1.1) would have to change with it.

I'd recommend (a). This is a physics choice for you and the reviewer; nothing here
changes the production code.

**Decision (Sep 28): option (c), |η^γ| < 0.7, for the next production.** See the
top-level `../README.md` "Decisions" section.

## R = 0.4 xJ with and without |η^γ| < 0.7 (`draw_xj_eta07.C`)

This is reco-level data (the unfolding input), uncorrected region A and purity-corrected.
- Nominal (|η| < 1.1) uses the committed purity.
- |η| < 0.7 uses `combine_hists` on the |η| < 0.7 data ABCD counts and Photon-MC leakage,
  taken from `purity_vs_eta.C`'s output.

Output is `pdfs/draw_xj_eta07.pdf` (page 1: corrected, page 2: uncorrected) and
`draw_xj_eta07.log`. Not unfolded: that needs the MC response rebuilt with the cut.

| pT [GeV] | P_A 1.1 → 0.7 | ⟨xJ⟩ uncorr. 1.1 → 0.7 | ⟨xJ⟩ corrected 1.1 → 0.7 |
|---|---|---|---|
| 15–20 | 0.62 → 0.68 | 0.878 → 0.882 | 0.852 ± 0.011 → 0.858 ± 0.011 |
| 20–25 | 0.78 → 0.79 | 0.842 → 0.845 | 0.823 ± 0.019 → 0.830 ± 0.021 |
| 25–35 | 0.81 → 0.86 | 0.850 → 0.838 | 0.822 ± 0.039 → 0.814 ± 0.042 |

The shapes agree bin by bin within statistics in every pT bin, corrected and uncorrected,
and the corrected ⟨xJ⟩ moves by ≤ 0.008. The central-only purity at 15–20 GeV is higher
(0.68 vs 0.62). The xJ shape barely changes because the higher purity is consistent with
dropping the less-pure forward photons, and their signal xJ looks like the central one.

## In-situ p_a with and without |η^γ| < 0.7 (`insitu_eta07.C`, all R)

**Method.** Every scan input is rebuilt from the trees with the cut applied consistently:
- data A/C at the raw jet scale, with no xJ pre-cut (extended scan, 0.80–1.05);
- Photon5+10+20 region-A reference ⟨xJ⟩;
- |η| < 0.7 purities from `combine_hists` per R;
- `grid_insitu.C`'s purity-corrected χ², verbatim.

**Validation (nominal).** The rebuilt data A/C counts equal the committed in-situ
inputs exactly at every R. MC reference means agree to ≤ 0.002 (photon-smearing random
sequence). The nominal p_a reproduces `ana.h` at every R.

Output is `pdfs/insitu_eta07.pdf` and `.root` (χ² curves), and `insitu_eta07.log`.

| R | p_a, \|η^γ\| < 1.1 | p_a, \|η^γ\| < 0.7 | shift | χ² (3 pts) 1.1 → 0.7 |
|---|---|---|---|---|
| 0.2 | 0.886 (production: 0.900 edge) | **0.903** −0.006 +0.010 | +0.018 | 14.4 → 8.2 |
| 0.3 | 0.903 −0.005 +0.007 | 0.916 −0.011 +0.007 | +0.012 | 3.1 → 2.3 |
| 0.4 | 0.922 ± 0.007 | 0.926 −0.007 +0.008 | +0.004 | 0.2 → 0.6 |
| 0.5 | 0.919 −0.008 +0.006 | 0.926 −0.007 +0.009 | +0.007 | 0.7 → 1.1 |
| 0.6 | 0.921 −0.007 +0.006 | 0.929 ± 0.007 | +0.009 | 0.1 → 0.5 |
| 0.7 | 0.919 −0.007 +0.008 | 0.930 ± 0.008 | +0.011 | 3.8 → 0.9 |
| 0.8 | 0.942 −0.009 +0.007 | 0.952 ± 0.009 | +0.010 | 2.1 → 0.7 |

(The R = 0.2 nominal upper error is not meaningful: the χ² there is jagged near its
minimum as events step across the xJ floor.)

**Findings**

1. **The cut raises p_a at every R**, by +0.004 (R = 0.4) to +0.018 (R = 0.2); +0.007 to
   +0.012 for the other radii. The |η| < 0.7 sample is a subset of the nominal one, so the
   shifts are more significant than the separate ±0.007 errors suggest.
2. **The shift comes from the purity, not the MC reference.** The γ+jet reference ⟨xJ⟩
   moves by ≤ 0.003 with the cut. The data purity at 15–20 GeV rises by 0.05–0.07 at
   every R. A higher purity means less background subtraction and a higher corrected
   data ⟨xJ⟩, so p_a rises. This matches the ~0.03 underestimate of the |η|-integrated
   purity (+0.5% p_a estimated for R = 0.4; +0.4% observed).
3. **With the cut, R = 0.2 is back inside the production scan range** (0.903) and its fit
   improves (χ² 14.4 → 8.2). R = 0.3 moves to 0.916. The forward photons' purity and
   background treatment is part of why the small-R scans sat at the edge (the xJ-floor
   sensitivity, `../insitu_scan_edge/` check 4, is the other part).
4. **The R dependence flattens.** With the cut, p_a is 0.903 / 0.916 / 0.926 / 0.926 /
   0.929 / 0.930 / 0.952, against 0.886–0.942 nominal. The spread over R = 0.3–0.7 shrinks
   from 0.019 to 0.014.
