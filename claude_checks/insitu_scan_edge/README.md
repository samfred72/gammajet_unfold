# Issue 3: in-situ scans stuck at the 0.900 scan edge

Reviewer (ppg18_checks.pdf, slide 5): the Jet12 reference scan sits at 0.9000 for every R
(note Fig. 10). The threejet scans at R = 0.2 and 0.3 and the JERhigh scan at R = 0.2 are
also at the edge, and they feed `jesTotalErr` for R = 0.2/0.3.

## Check 1: signal vs background xJ in Jet12_long (`draw_jet12_sig_bkg_xj.C`)

**Question.** The purity-corrected in-situ mean (`insitu_utility::computeCorrectedMeans`)
and the unfolding input (`unfold_utility::purityCorrect`) assume

    A = P_A·S + (1−P_A)·B,    C = P_C·S + (1−P_C)·B

with one background xJ shape B common to regions A and C. How different is the
background xJ from the signal xJ, and does region C's background reproduce region A's?
This matters most at R = 0.2/0.3.

**Method.** The input is Pythia8 `gammajet_pythia_Jet12_long.root`.
- Each region-A/C reco cluster is labelled:
  - **truth-tagged** (signal) if it is within ΔR < 0.1 of the tree's truth photon
    (`unfolder::check_match`). That is the purity method's own signal definition
    (`hclusterpt_abcd_truthmatched`, `unfolder.cc:317`). xJ stays reco level throughout.
  - **background** otherwise.
- Kinematics and weights are the nominal MC branch of `unfolder.cc`: EMR-smeared
  photon pT, `jet_pt_smear_truth[ir]`, the `check_pair` cuts including the xJ floor, the
  vz/cluster-pT `mcWeight`, and the Jet12_long truth-jet window.
- **Bias.** The MC component means are mixed with the **data** purities
  (`hists/purity_nominal.root`, per R) and passed through the analysis coefficients
  (`purityCorrectCoeffs`). The signal in C is set equal to S_A, so only the B_A ≠ B_C
  mismatch enters.

Output is in `pdfs/draw_jet12_sig_bkg_xj.pdf` (pages 1–3: shapes per R; page 4: means)
and `draw_jet12_sig_bkg_xj.log`.

**⟨xJ⟩ at 15–20 GeV** (errors about ±0.006 on S_A and B_C, ±0.012 on B_A and S_C):

| R | S_A | B_A | S_C | B_C | B_A − S_A | B_A − B_C | bias on corrected ⟨xJ⟩ |
|---|-----|-----|-----|-----|-----------|-----------|-------------------|
| 0.2 | 0.740 | 0.770 | 0.738 | 0.766 | +0.030 | +0.004 | +0.4 ± 1.3% |
| 0.3 | 0.795 | 0.850 | 0.792 | 0.829 | +0.055 | +0.021 | +1.8 ± 1.3% |
| 0.4 | 0.837 | 0.906 | 0.842 | 0.892 | +0.069 | +0.014 | +1.2 ± 1.2% |

At 20–25 GeV the pattern is the same (B_A − S_A = +0.09 / +0.13 / +0.11, bias +1.6 / +2.5
/ +1.4%). At 25–35 GeV only 20–80 effective entries per component are left, so the
biases are 0 within ±3%.

**Findings**

1. **Background xJ is higher than signal xJ.** For a fake photon (a leading π0 or η from
   a jet), the recoil jet balances the whole parent jet, not just the neutral meson. The
   gap grows with R, from +0.03 at R = 0.2 to +0.07 at R = 0.4 (15–20 GeV). At R = 0.2
   background and signal are closest, so the background treatment matters least there.
2. **Region C's background reproduces region A's to within 0.005–0.02 in ⟨xJ⟩** at 15–20
   GeV, with B_C slightly lower. Most shape comparisons are compatible: KS probability
   0.13–0.95 at 15–20 GeV. The exception is R = 0.3 at 20–25 GeV (0.03). Signal in C matches signal in
   A (S_C ≈ S_A in every bin), so the leakage term is well behaved.
3. **The resulting bias pushes p_a in the wrong direction to explain the edge.** The
   B_A ≠ B_C mismatch raises the purity-corrected data ⟨xJ⟩ by +0.4–2.5%, which would pull
   the fitted p_a *up*. The scans are stuck at the *lower* edge: the χ² wants
   p_a < 0.900. The effect is also smallest at R = 0.2 (+0.4%). So the region-C
   background model is **not** why R = 0.2/0.3 hit 0.900.
4. **Jet12_long is not a background-only sample.** About 77% of its region-A pairs
   (after `mcWeight`) are truth-tagged photons, and about 18% of its region-C pairs. `grid_insitu_jet12.C` and `draw_insitu_xj_jet12.C` describe it as "mostly
   fake-photon-triggered dijet background". The truth photon in the tree
   (`CaloAna.cc::ProcessFillTruthPhotonParticle`, `b78f36d`) is the leading non-hadron-decay
   photon above 10 GeV, so parton-shower photons count. Either those, or prompt-photon
   processes in the Jet12_long generator config, make up this fraction. That changes how Fig. 10 should be read (it isn't a like-for-like
   "raw data A vs raw background-dominated MC A"). **To check:** the Pythia config
   used for Jet12_long on SDCC.

**Caveats.** This is MC-only: the fake-photon composition in data can differ from
Pythia. Statistics above 20 GeV are thin. The purities used are nominal data purities
per R.

## Check 2: uncorrected vs purity-corrected xJ, closure on Jet12 (`draw_jet12_purity_corrected_xj.C`)

**Method.** Same selection as check 1, R = 0.2/0.3/0.4. The whole production chain is run
on Jet12_long:
- **ABCD counts:** paired clusters per photon-pT bin from this sample.
- **Leakage fractions:** truth-matched Photon5+10+20 templates
  (`hclusterpt_abcd_truthmatched`, stitched by `drawer::get`).
- **Purity:** `puritymaker.C::combine_hists` itself (`#include`d, not copied), i.e. the
  leakage-corrected quadratic, the bootstrap median and P_C = c·S/C.
- **Correction:** `unfold_utility::purityCorrect(A, C, P_A, P_C)`.

The corrected region-A xJ is compared with the **truth-tagged reco photons in region A**:
region-A clusters within ΔR < 0.1 of the tree's truth photon, the purity method's own
signal definition (`unfolder.cc:317`). This is still reco-level xJ, i.e. the true-photon
subset that the subtraction should leave, not a particle-level distribution.

**Selection cross-check.** The ABCD counts filled here agree with the committed
`hists/Jet12_long_pythia_nominal_unfolding.root` `hclusterpt_abcd` to ≲ 1% at 15–20 GeV
and ≲ 4% above. The residual comes from the photon-smearing random sequence at the
bin edges.

Output is `pdfs/draw_jet12_purity_corrected_xj.pdf` (one page per R) and
`draw_jet12_purity_corrected_xj.log`. Means below are binned; errors are ×10⁻³.

| R | pT [GeV] | P_A (ABCD) | P_A (true) | P_C | ⟨A⟩ uncorr. | ⟨corrected⟩ | ⟨truth-tagged⟩ | corrected/tagged − 1 |
|---|---|---|---|---|---|---|---|---|
| 0.2 | 15–20 | 0.739 | 0.766 | 0.172 | 0.749(6) | 0.743(9) | 0.742(6) | +0.1 ± 1.4% |
| 0.2 | 20–25 | 0.783 | 0.795 | 0.271 | 0.692(13) | 0.681(20) | 0.673(14) | +1.2 ± 3.6% |
| 0.2 | 25–35 | 0.819 | 0.809 | 0.362 | 0.713(29) | 0.712(43) | 0.683(32) | +4.3 ± 7.8% |
| 0.3 | 15–20 | 0.742 | 0.772 | 0.171 | 0.810(6) | 0.803(9) | 0.797(6) | +0.7 ± 1.3% |
| 0.3 | 20–25 | 0.783 | 0.801 | 0.269 | 0.774(13) | 0.763(20) | 0.749(14) | +2.0 ± 3.3% |
| 0.3 | 25–35 | 0.849 | 0.799 | 0.404 | 0.775(26) | 0.754(39) | 0.754(29) | −0.0 ± 6.4% |
| 0.4 | 15–20 | 0.747 | 0.772 | 0.179 | 0.854(6) | 0.841(9) | 0.838(6) | +0.3 ± 1.3% |
| 0.4 | 20–25 | 0.779 | 0.795 | 0.271 | 0.827(14) | 0.815(21) | 0.806(15) | +1.2 ± 3.3% |
| 0.4 | 25–35 | 0.843 | 0.835 | 0.406 | 0.828(32) | 0.797(48) | 0.803(33) | −0.8 ± 7.2% |

**Findings**

1. **The purity chain closes on Jet12 at every R.** At 15–20 GeV the corrected ⟨xJ⟩
   matches the truth-tagged photons to 0.1–0.7% (±1.3%). The shapes agree bin by bin in
   the core (pages 1–3). At 20–35 GeV there is agreement within the (larger) errors. The
   ABCD purity is 0.74–0.75 against a true 0.77 at 15–20 GeV, a 3–4% underestimate. It
   partly offsets the B_A > B_C mismatch from check 1, which is why the closure is
   better than check 1's bias estimate alone suggests.
2. **The correction is small in this sample.** It lowers ⟨xJ⟩ by only 0.6–1.3% at 15–20
   GeV, because Jet12_long's region A is already ~77% signal (check 1, finding 4). It
   therefore does not stress the method the way data does (P_A ≈ 0.61 at 15–20 GeV). Data
   subtracts about twice as much background, and check 1's +0.4–2.5% bias estimate with data
   purities is the relevant number there.
3. Nothing here pushes the corrected mean down at R = 0.2/0.3, so the purity correction
   is not what drives the scans to the 0.900 edge.

## Check 3: data purity and corrected vs uncorrected xJ, next to Jet12 (`draw_data_vs_jet12_purity_xj.C`)

**Method.** The data side uses exactly the in-situ scan's ingredients (`grid_insitu.C`):
- data A and C from `insitu/inputs/Data_nominal_insitu.root` via
  `insitu_utility::cacheDataEvents`, at the **raw jet scale** (pa = 1);
- the `lowXjFloor` floor;
- the committed per-R data purities (`ana::getPurity`/`getPurityC`,
  `hists/purity_nominal.root`);
- `purityCorrectByPtBin` for the correction;
- the Pythia8 γ+jet region-A reference (`buildMCXjByPtBin`, Photon5+10+20, same weights).

The Jet12 numbers are read from check 2's output. Output is in
`pdfs/draw_data_vs_jet12_purity_xj.pdf` (page 1: purities; pages 2–4: xJ per R, with
corrected/uncorrected ratios for data and Jet12) and `draw_data_vs_jet12_purity_xj.log`.

**Purity, P_A (P_C)**

| R | 15–20 data | 15–20 Jet12 ABCD (true) | 20–25 data | 20–25 Jet12 | 25–35 data | 25–35 Jet12 |
|---|---|---|---|---|---|---|
| 0.2 | 0.61 (0.15) | 0.74 (0.17) [0.77] | 0.77 (0.29) | 0.78 (0.27) | 0.77 (0.24) | 0.82 (0.36) |
| 0.3 | 0.62 (0.15) | 0.74 (0.17) [0.77] | 0.78 (0.30) | 0.78 (0.27) | 0.77 (0.25) | 0.85 (0.40) |
| 0.4 | 0.62 (0.15) | 0.75 (0.18) [0.77] | 0.78 (0.30) | 0.78 (0.27) | 0.81 (0.27) | 0.84 (0.41) |

**Binned ⟨xJ⟩, corrected − uncorrected, and data/MC**

| R | pT | data uncorr. → corr. (shift) | Jet12 shift | γ+jet MC | data/MC uncorr. → corr. |
|---|---|---|---|---|---|
| 0.2 | 15–20 | 0.726 → 0.705 (−0.021) | −0.006 | 0.751 | 0.966 → 0.939 |
| 0.2 | 20–25 | 0.662 → 0.654 (−0.008) | −0.011 | 0.704 | 0.941 → 0.930 |
| 0.2 | 25–35 | 0.650 → 0.596 (−0.054 ± 0.04) | −0.002 | 0.714 | 0.910 → **0.834** |
| 0.3 | 15–20 | 0.775 → 0.749 (−0.026) | −0.007 | 0.805 | 0.962 → 0.930 |
| 0.3 | 20–25 | 0.734 → 0.721 (−0.012) | −0.011 | 0.769 | 0.953 → 0.938 |
| 0.3 | 25–35 | 0.742 → 0.711 (−0.031) | −0.021 | 0.783 | 0.948 → 0.908 |
| 0.4 | 15–20 | 0.823 → 0.798 (−0.024) | −0.013 | 0.852 | 0.965 → 0.937 |
| 0.4 | 20–25 | 0.777 → 0.759 (−0.018) | −0.012 | 0.823 | 0.944 → 0.922 |
| 0.4 | 25–35 | 0.789 → 0.764 (−0.026) | −0.031 | 0.833 | 0.947 → 0.917 |

**Findings**

1. **Data purity barely depends on R** (≤ 0.015 spread at a given pT). The purity input
   cannot create an R = 0.2-specific pull.
2. **At 15–20 GeV the data sample is much less pure than Jet12**: P_A = 0.61 vs 0.74
   (0.77 true). So the data correction is 2–4× larger than the one Jet12 validated: a
   shift of −0.021 to −0.026 in ⟨xJ⟩ against −0.006 to −0.013. Check 1's bias estimate
   with data purities (+0.4 to +1.8% at 15–20 GeV) applies here, and it would make the
   corrected data mean slightly too *high*. At 20–25 GeV data and Jet12 have the same
   purity and similar shifts. At 25–35 GeV data P_C (0.24–0.27) is lower than Jet12's
   (0.36–0.41), though both have large errors.
3. **After correction, the per-bin data/MC ratio is ~0.93 at R = 0.2 at 15–25 GeV, the
   same as at R = 0.3 and 0.4.** The one R = 0.2 outlier is 25–35 GeV (0.834 ± ~0.06):
   its uncorrected ratio is already the lowest (0.910), and the correction subtracts a
   further 0.054, twice the R = 0.3/0.4 shift in that bin.

   An error-weighted average of the three R = 0.2 ratios is ≈ 0.93, well above the
   scan's 0.900 edge. So the ratio of corrected means does *not* by itself ask for
   p_a < 0.90 at R = 0.2. The edge comes from how the scan responds to p_a, not from
   the purity-corrected inputs at p_a = 1.
4. **A likely mechanism, not yet tested: the xJ floor.** At R = 0.2 and 15–20 GeV the
   xJ distribution is largest in the floor bin (0.4–0.5; page 2). Lowering p_a scales
   jets up and pulls more events in across the floor, so ⟨xJ⟩ rises more slowly than
   1/p_a and the χ² keeps falling toward the scan edge. R = 0.3/0.4 peak well above the
   floor. Printing corrected-data ⟨xJ⟩ vs p_a (0.80–1.05) per R and pT bin would confirm
   or rule this out without touching the production scan.

## Check 4: the scan beyond its range (`scan_extended_range.C`)

**Method.** The production scan can't look below 0.900: `Data_*_insitu.root` only keeps
events with raw xJ ≥ 0.90 × floor (`unfolder.cc:287`, `floorScale = scanLow`). Here data
A/C are rebuilt from `trees/gammajet_Data.root` with no xJ pre-cut. The reference,
purities and χ² are identical to `grid_insitu.C`'s purity-corrected fit, and p_a is
scanned over 0.70–1.10.

**Cross-check.** With raw xJ ≥ 0.90 × floor, the rebuilt A/C counts equal the production
input exactly at every R (e.g. R = 0.2: A 4272 / C 3919). The in-range minima reproduce
`ana.h`: 0.900 (edge), 0.903, 0.922.

The **fixed-population** variant applies the floor once at p_a = 1, so ⟨xJ⟩ scales
exactly as 1/p_a. Comparing it with production isolates the effect of the floor moving
with p_a. Output is `pdfs/scan_extended_range.pdf` and `scan_extended_range.log`.

| R | in-range min | full-range min (χ²) | fixed-population min (χ²) | per-bin min: 15–20 / 20–25 / 25–35 GeV |
|---|---|---|---|---|
| 0.2 | 0.900 (edge, χ² 17.8) | **0.886** (15.2) | 0.927 (22.7) | 0.899 / 0.928 / 0.806 |
| 0.3 | 0.903 (2.7) | 0.903 (2.7) | 0.929 (1.1) | 0.899 / 0.929 / 0.894 |
| 0.4 | 0.922 (0.3) | 0.922 (0.3) | 0.934 (0.8) | 0.923 / 0.924 / 0.910 |

**Findings**

1. **The R = 0.2 minimum is real and just outside the range, at 0.886.** Widening
   `scanLow` would move the result there, not much lower. The fit is poor, though:
   χ² = 15 for three points. The pT bins disagree: 15–20 GeV wants 0.899, 20–25 GeV wants
   0.928 (the same as every other R), and 25–35 GeV wants 0.806. The 25–35 GeV bin is
   the one whose purity correction was anomalously large (check 3). At R = 0.3/0.4 the
   bins agree.
2. **The moving floor lowers the fit, and more so at small R.** Production minus fixed
   population: −0.041 (R = 0.2), −0.026 (R = 0.3), −0.012 (R = 0.4). Small-R xJ
   distributions sit on the floor (R = 0.2 at 15–20 GeV peaks in the floor bin, check 3).
   As p_a drops, events cross the floor from below, so ⟨xJ⟩ rises more slowly than 1/p_a
   and the fit is pushed down. Applying the floor after scaling is the right procedure,
   because data and MC are then compared in the same corrected-xJ window. But it makes
   the small-R fit depend on how well data and MC agree *at the floor*, which is where
   modelling is weakest (see below).
3. **The 15–20 GeV bin is what pins R = 0.2 and R = 0.3 near 0.90** (both want 0.899),
   while 20–25 GeV gives ~0.93 at every R. At R = 0.4 15–20 GeV agrees with the rest
   (0.923).

**Hypotheses still open, with suggested tests**

- **JER smearing is not radius-specific.** `CaloAna.cc` smears every radius's MC jets
  with the same `h_jerband_quaddiff` (`quaddiff_bi_nominal.root`), so R = 0.2 MC jets get
  the extra resolution derived for one radius. A resolution mismatch changes how many
  events sit near the floor. The JERhigh scan also pins at R = 0.2. Test: compare the
  data vs MC xJ width at R = 0.2, or rescan with the JER smear scaled.
- **Low-pT jet response or efficiency near threshold.** Near the floor (15–20 GeV,
  xJ 0.4–0.5 ⇒ 6–10 GeV jets) the data/MC jet response and efficiency at the 5 GeV
  threshold matter most. That is also where the data-only timing window loses the most
  jets (24% at 5–9 GeV). The timing loss raises data ⟨xJ⟩, which is the wrong direction
  for the pull, but it distorts the population at the floor the fit is sensitive to.
  Test: rescan with a higher floor, e.g. jet pT > 7–8 GeV.
- **A genuine small-R JES difference**, e.g. out-of-cone energy or jet-shape
  mismodelling that affects narrow jets most. That fits the monotonic `jesNominal` trend
  with R (0.90 → 0.94), but not the 20–25 GeV bin agreeing at 0.93 for every R.
- **The R = 0.2, 25–35 GeV purity correction** (P_C = 0.24, shift −0.054) drags the
  combined minimum from ~0.90 to 0.886. It is a low-statistics bin.

## Open for this issue

- **Revisit at the data reprocess (agreed Sep 28):** the MC jet smearing
  (`CaloAna.cc::smear_pt`, one `h_jerband_quaddiff` for every R at `b78f36d`) may need a
  radius-specific JER. Check whether the current SDCC code already differs, and whether
  the data R = 0.2 xJ width supports a different smear.

- Findings 1–3 rule out the ABCD background model as the cause of the scan edge at small
  R. Candidates still to examine:
  - a genuine data/MC JES gap below 0.90 at R = 0.2 (nominal R = 0.2 has χ² = 22, against
    ≤ 4 at other R);
  - the threejet and JERhigh variations themselves at R = 0.2/0.3;
  - behaviour of the xJ floor near the scan edge (the floor is tied to `scanLow`).

  The direct test is to widen `insitu_utility::scanLow` (e.g. to 0.80) and rerun
  `run_grid.sh`; the unfolder's `ispairedInsitu` floor follows `scanLow`, so the
  inputs must be regenerated as well.
- Resolve finding 4 before regenerating Fig. 10.
