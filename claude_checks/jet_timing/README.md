# Issue 1: jet timing window `0 < t_MBD - t_jet < 4 ns`

Reviewer (ppg18_checks.pdf, round 2, slide 6 / question 1): the window is applied to data
only (gammajet_treemaking `CaloAna.cc:433`, `b78f36d`). The recoil-jet Δt distribution is
truncated at 0 ns. A Gaussian extrapolation gives 8–20% loss if the edge entries are
in-time jets, and ~0 if they are the out-of-time background the cut is meant to remove.
The note should document the cut's efficiency for in-time jets.

## How the cut works (read from the source, `b78f36d`)

- `t_jet` is the energy-weighted mean time of the jet's constituent towers with
  E > 0.1 GeV, over EMCal, iHCal and oHCal, times 17.6 ns/sample
  (`CaloAna.cc:393-432`). A jet with no such tower gets `-999`, which fails the window.
- In data, a jet outside the window is **skipped**, and the highest-pT *in-window* jet
  becomes `jet_*[ir]` (`CaloAna.cc:433-437`). A lost recoil jet is therefore either
  replaced by a softer jet or the event loses its pairing. The trees never show which.
- The photon cluster has the same window (`CaloAna.cc:177`). MC has neither cut, and MC
  timing is not a usable reference: its Δt_cluster peaks at −0.7 ns with RMS 2.0 ns,
  against +1.9 ns / 0.9 ns in data.

Because the trees only hold jets that passed, the removed jets can't be counted from the
current trees. What *can* be tested is whether the jets piled against the edges are
in-time recoil jets or background.

## Checks on the existing trees (`draw_jet_timing.C`, R = 0.4, Sep-1 data tree)

Selection follows `unfolder.cc`: `pho_object`/`findabcdBin` regions, |v_z| < 60 cm,
15 < p_T^γ < 35 GeV, the `check_pair` cuts and the xJ floor, and jet pT =
`jet_pt_calib / jesNominal`. Output is in `pdfs/draw_jet_timing_r04.pdf` and `.root`; the
numbers below are from `draw_jet_timing_r04.log`.

1. **The edge jets are recoil jets.** Out-of-time or noise jets would be uncorrelated
   with the photon in φ, putting 12.5% of them above Δφ = 7π/8. With no Δφ cut (region
   A), the fraction above 7π/8 is:
   - 0.794 ± 0.016 for jets at 0 < Δt < 0.5 ns (N = 630)
   - 0.834 for the core, 1 < Δt < 3 ns (N = 3192)
   - 0.798 at the upper edge, 3.5 < Δt < 4 ns (N = 173)

   KS probability between the edge and core Δφ shapes is 0.37. Modelling the edge sample
   as core-like signal plus φ-flat background gives a background fraction of
   **5.6 ± 2.5%**. The entries at the 0 ns edge are ≳ 90% in-time recoil jets, not the
   background the cut targets (page 7).
2. **The Δt offset is set by calorimeter composition.** Mean in-window Δt_jet goes from
   **2.85 ns** for EM fraction < 0.3 (HCal-dominated) to **1.09 ns** for EM fraction >
   0.85 (EMCal-dominated) (pages 3, 5). The photon cluster, which is pure EMCal, sits at
   1.9–2.0 ns. The EMCal and HCal tower-time scales are offset by ~2 ns relative to each
   other, and a jet's time is a mixture. The window clips **both** edges:
   - EM-rich jets are cut at 0 ns.
   - HCal-rich jets are cut at 4 ns; that distribution is still at ~60% of its peak
     height there (page 3).
3. **The loss depends on jet pT.** In a Gaussian model fitted inside the window, the
   fraction outside [0, 4] ns falls from 26 ± 4% at 5–8 GeV to 12% at 14–24 GeV and 6%
   at 24–35 GeV (page 9). By EM fraction it is 23% for HCal-dominated jets and 7–11%
   otherwise (page 10). The model is crude: the EM-fraction slices are not Gaussian.
   Treat these as the size of the effect, not a correction.
4. **Implied xJ shape bias.** Reweighting region-A reco-level xJ by
   1/(1 − loss(p_T^jet)) changes the normalized shape by **+13 to +16% at the lowest xJ
   bin and −4 to −9% at high xJ**, similar in all three photon-pT bins (page 11). The
   direction is expected: low-xJ pairs have softer jets, which are lost more often.

   This is not a correction. It treats a lost jet as a lost event (really a softer
   in-window jet may be substituted), it is not unfolded or purity-corrected, and the
   in-situ JES is derived from the same windowed data, so part of the bias may be
   absorbed into p_a. But it is the same order as the listed systematics (JES ~7–15%
   RMS), and nothing in the current systematic budget covers it.
5. Consistency with the reviewer's numbers: 207 of 4257 region-A jets fall in the first
   0.25 ns of the window, against 10 clusters. The reviewer found 250/4246 vs 11,
   probably from a slightly different selection.

**Conclusion:** the "data-only by design, removes background" reading isn't supported.
The edge population is signal, and the cut removes in-time recoil jets with a jet-pT- and
composition-dependent efficiency that MC doesn't model. Documenting the cut is not enough
on its own.

An aside on the stated purpose: RHIC Run-24 p+p bunches are ~106 ns apart, so jets from
other crossings would sit near Δt ≈ ±106 ns. A window many ns wide removes them just as
well as 0–4 ns does.

## What is needed to close it

The loss has to be *measured*, which needs the jets that are cut today.

- **`treemaking_jet_notimecut.patch`** (against gammajet_treemaking `b78f36d`; checked
  with `git apply --check`, **not yet compiled** because SDCC `/sphenix/user` was
  returning stale file handles).
  - It adds data-only branches `jet_nt_{pt,pt_calib,eta,phi,emfrac,time}[7]`: the
    leading jet after the same pT/ΔR selection, with **no** timing window.
  - All existing branches and the current selection are unchanged, so the current
    result stays reproducible from the new trees.
  - After reprocessing the data trees (MC doesn't need it), this measures directly:
    - the efficiency: among region-A events whose no-window leading jet is
      back-to-back, the fraction inside the window, vs jet pT and EM fraction. This is
      the figure the reviewer asks for, measured rather than modelled.
    - how often the stored `jet_*` differs from `jet_nt_*`, i.e. how often the pairing
      jet is substituted.
    - the unfolded result with a wider window, or with no jet window, as a systematic
      or a new nominal.
- **Longer term:** equalize the EMCal/HCal tower-time offsets, or cut per calorimeter,
  before applying any narrow window.

The analysis-side code (the efficiency plot and a configurable window in `unfolder.cc`)
should be written once reprocessed trees exist. Nothing here has been run in production.

## Resolution (agreed Sep 28)

- **The Sep-1 data tree has the window applied.** `trees/gammajet_Data.root`, the same
  file as SDCC `gammajet_full_hadded/gammajet_Data.root` (465,273,712 bytes, Sep 1), has
  0 of 714,009 R=0.4 jets, 0 of 599,338 R=0.2 jets and 0 of 1,253,057 clusters outside
  0 < Δt < 4 ns, with sharp edges at both limits. Whatever the current SDCC treemaking
  code does, the result in the note was built with the cut. The findings above describe
  that sample.
- **For the data reprocess, keep "leading object that passes the timing window".** That
  is the correct choice when the window rejects only genuinely out-of-time objects (other
  crossings at about ±106 ns, noise): the next in-time object really is the collision's
  leading jet or photon. The problem today is the narrow 0–4 ns window, not the
  leading-passing choice. The narrow window rejects in-time jets (the EMCal/HCal offset
  plus resolution), and the next-jet swap then happens in data only and is not modelled
  in MC.
- **What to change in treemaking:**
  - Make the jet and cluster windows wide, several ns either side of the in-time peak.
  - Ideally, correct the EMCal/HCal time offset first: store the EMCal and HCal parts of
    the jet time separately.
  - Also store the leading jet with no timing requirement (the `jet_nt_*` branches in
    `treemaking_jet_notimecut.patch`, or equivalent in the current SDCC code). This
    measures how often "leading" and "leading-in-time" differ, and whether the discarded
    objects look out-of-time: flat in Δφ, far from the peak. It is also the efficiency
    figure for §3.2.
- **Downstream:** tighter windows can be applied in `unfolder.cc` as a systematic
  variation. Tightening downstream can only drop events; it can't re-select a softer
  object.

## Radius dependence (`draw_jet_timing_vs_R.C`, all 7 radii, Sep-1 data tree)

This uses the same selection as `draw_jet_timing.C`, evaluated per radius (the treemaking
window is applied per R, and each R picks its own leading in-window jet). Region A,
paired, 15–35 GeV. Output is `pdfs/draw_jet_timing_vs_R.pdf` and
`draw_jet_timing_vs_R.log`.

| R | N | ⟨Δt⟩ in window [ns] | RMS | in first 0.25 ns | ⟨EM frac.⟩ | loss 5–9 / 9–13 / 13–18 / 18–24 / 24–35 GeV (Gaussian model) |
|---|---|---|---|---|---|---|
| 0.2 | 4272 | 1.73 | 0.97 | 4.2% | 0.65 | 0.24 / 0.18 / 0.13 / 0.09 / 0.05 |
| 0.3 | 4521 | 1.69 | 0.95 | 4.4% | 0.64 | 0.24 / 0.17 / 0.12 / 0.10 / 0.06 |
| 0.4 | 4257 | 1.62 | 0.91 | 4.9% | 0.63 | 0.22 / 0.17 / 0.13 / 0.11 / 0.06 |
| 0.5 | 3816 | 1.55 | 0.90 | 5.8% | 0.63 | 0.36 / 0.20 / 0.14 / 0.12 / 0.07 |
| 0.6 | 3252 | 1.49 | 0.90 | 6.5% | 0.63 | (fit unstable) / 0.22 / 0.17 / 0.13 / 0.08 |
| 0.7 | 2602 | 1.45 | 0.90 | 7.5% | 0.63 | (fit unstable) / 0.30 / 0.24 / 0.14 / 0.10 |
| 0.8 | 1939 | 1.43 | 0.91 | 8.0% | 0.63 | (fit unstable) / 0.47 / 0.31 / 0.17 / 0.12 |

N falls with R because of the jet acceptance |η_jet| < 1.1 − R.

**Findings**

1. **There is a radius dependence: the jet-time distribution moves toward the 0 ns edge
   as R grows.** ⟨Δt⟩ drops from 1.73 ns (R = 0.2) to 1.43 ns (R = 0.8), and the pile-up
   in the first 0.25 ns doubles (4.2% → 8.0%). Small-R jets instead have the largest tail
   toward the 4 ns edge (page 1).
2. **It isn't a change in calorimeter mix.** ⟨EM fraction⟩ is flat (0.63–0.65). At
   *fixed* EM fraction, EM-dominated jets have the same Δt at every R (≈ 1.25 ns for EM
   fraction > 0.85). HCal-dominated jets shift strongly with R: 2.61 ns (R = 0.2) →
   1.97 ns (R = 0.8) for EM fraction < 0.3 (page 4). A larger cone adds more low-energy
   peripheral towers above the 0.1 GeV threshold, and their energy-weighted times pull
   the average. That looks like an amplitude-dependent (time-walk) effect in the HCal
   tower times, on top of the EMCal/HCal offset. The window is therefore not equally
   placed for every R.
3. **The implied loss is about the same for R = 0.2–0.4 and grows for R ≥ 0.5.** At
   13–18 GeV it is 0.12–0.13 for R ≤ 0.4, then 0.14 / 0.17 / 0.24 / 0.31 for
   R = 0.5 / 0.6 / 0.7 / 0.8. At 24–35 GeV it is 0.05–0.06 for R ≤ 0.4 and 0.07–0.12 above.
   At R ≥ 0.6 the 5–9 GeV Gaussian fits are unstable (too few low-pT jets for a
   reliable truncated fit) and should not be read as losses.
4. **For issue 3 (R = 0.2/0.3 scans at 0.900):** R = 0.2 and 0.3 have the same
   pT-dependent loss as R = 0.4, so the timing window does not single out small R. As
   at R = 0.4, the loss removes low-pT jets, which raises data ⟨xJ⟩, the opposite of what
   would pull p_a down. For the R ≥ 0.5 results, the timing window is a larger and
   R-dependent data-only effect.
5. **For the reprocess:** a single fixed window cannot be equally efficient at every R.
   A window wide enough for the largest R, or a jet time built with a higher tower-energy
   threshold (or energy-squared weighting) to suppress low-amplitude tower times, removes
   most of this R dependence. Store the per-calorimeter times so it can be studied
   downstream.

Same caveat as above: this is from in-window entries only, and the loss uses the Gaussian
model. The direct measurement needs the no-window jet from the reprocess.

## Tower-level timing study (TimingAna, Oct 1)

To find out what drives the jet time, a separate Fun4All sub-project on SDCC,
`/sphenix/user/samfred/projects/gammajet/timingana`, stores the leading R = 0.4 jet with
**no** timing cut and no vertex or η cuts.
- **Code:** `src/TimingAna.{cc,h}`, `macros/Fun4All_timing.C`, `macros/run_one.sh`. It
  builds `libtimingana` into the shared install.
- **Contents of the tree:**
  - MBD time and vz
  - jet-level time per calorimeter
  - every constituent tower's E, t, χ² and status
  - every raw EMCal tower under the jet's retowers
- **Sample:** one DST, `DST_JETCALO_run2pp_ana521_2025p007_v001-00047289-00000` (6,366
  events, 5,566 with a leading jet).
- **Plotting:** `draw_timingana.C` writes `pdfs/draw_timingana_r04.pdf`, with the numbers
  in `draw_timingana_r04.log`. It selects |vz| < 60, a finite MBD time,
  p_T^jet > 5 GeV and |η| < 0.7, which leaves 1,507 jets.

**Findings**

1. **Zero-suppressed EMCal towers carry t = 0 exactly, and the retower time includes
   them.**
   - ZS towers (status bit 5) aren't waveform-fit, but they are still flagged `isGood`.
   - ZS share of raw towers: 100% below 0.05 GeV, 87% at 0.05–0.1 GeV, 6% at
     0.1–0.2 GeV.
   - `RetowerCEMC` energy-weights them into the retower time. As a result, **54% of EMCal
     retowers above 0.1 GeV have t = 0 exactly** (77% at 0.1–0.3 GeV, ~0 above 0.3 GeV).
   - CaloAna's `jet_time` is built from retowers, so its EMCal part is pulled toward
     t = 0. Δt_EMCal is ~0.05 ns from retowers vs ~0.85 ns from raw non-ZS towers at
     11–13 GeV, and ≈ −0.4 vs 0.66 ns for EM fraction < 0.1.
2. **The non-ZS EMCal time is flat in jet pT and in EM fraction:** Δt ≈ 0.5–1.0 ns.
3. **The HCal time, dominated by the oHCal, rises with jet pT:**
   - from Δt ≈ 1 ns at 5–9 GeV to 3.3–3.8 ns at 13–18 GeV;
   - with falling EM fraction, it reaches 4.3 ns below 0.1.
   - The MBD-independent difference t_EMCal − t_HCal also grows, from 0.2 to 1.7 ns, so
     this is in the calorimeter, not MBD jitter. The jet-by-jet EMCal/HCal correlation
     is only 0.29.
4. **The cause is the oHCal tower time's dependence on tower energy:**
   - Δt_tower ≈ 1 ns at 0.1–0.4 GeV,
   - rising to ≈ 6–8 ns at 1–1.4 GeV,
   - then back down to ≈ 2.5–3 ns above 3 GeV.
   - Not monotonic, unlike simple time-walk, so it looks like a waveform-fit or
     calibration artifact.
   - Higher-pT and HCal-heavy jets have more oHCal towers in the 0.5–2 GeV range, which
     produces the pT and EM-fraction dependence of the jet time and pushes those jets
     toward the 4 ns edge.
   - The EMCal shows only mild walk: 0.1 → 1.3 ns over 0.1–4 GeV.

**Caveats:** one DST, mostly low-pT jets (only ~40 above 18 GeV). The iHCal has few towers
above 0.1 GeV.

**Implications for the reprocess:**
- Build the EMCal jet time from raw EMCal towers with ZS excluded, not from retowers.
- Store the EMCal and HCal times separately.
- Correct the oHCal tower time for its energy dependence, or use a wide window, before
  cutting.
- The oHCal structure needs more statistics (more DSTs) and a per-channel or gain check
  first.

## Dading Chen's tower-time calibration applied (Oct 2)

**Sample and setup:**
- 20 condor jobs (cluster 189721), segment 0 of 20 runs spread over his 284 calibrated runs (47352–52755).
- About 2M events, 462,517 leading R=0.4 jets. After |vz| < 60, a finite MBD time, pT > 5 GeV and |η| < 0.7, 120,792 remain.
- The calibration is the release sidecar `caloreco/CaloTowerTimeCalibration` with his run-by-run payloads.
- Towers enter with E > 0.5 GeV, EMCal ZS towers excluded, and the same towers in both versions.

**Outputs:** `draw_timingana_tcal.C` writes `pdfs/draw_timingana_tcal_r04.pdf` and `draw_timingana_tcal_r04.log`. `tcal_window_eff.C` writes `tcal_window_eff.log`.

**Findings**

1. **Offsets and slew are fixed.**
   - Median tower times go from −1.2 ns (EMCal) and −3 ns (oHCal) to ≈ 0 in every calorimeter.
   - EMCal time is flat in tower energy after correction.
   - **The median jet time is flat in EM fraction:** −3.4 → −1.5 ns standard, ≈ 0 corrected. It is also flat in pT up to 15 GeV (−0.6 ns at 16–20 GeV).
2. **The pT and EM-fraction dependence of the window efficiency mostly goes away.**
   - Fraction of jets with |t_MBD − t_jet − median| < 2 ns, by EM fraction (0–0.2 → 0.8–1): 0.47 → 0.81 standard, 0.73 → 0.82 corrected.
   - By pT (5 → 15 GeV): 0.75 → 0.62 standard, 0.80 → 0.74 corrected.
   - The half 16–84% width of t_MBD − t_jet goes from 1.78 to 1.61 ns (1.56 with outlier rejection).
3. **The one-sample-shifted oHCal towers are NOT removed.**
   - About 18–20% of oHCal towers at 1–1.5 GeV still have |t| > 9 ns. Overall: 10.5% standard, 8.5% corrected.
   - Group alignment would have removed a group-shift effect, so this is something else, probably the waveform fit at that amplitude.
   - Dropping towers with |t_corr| > 9 ns (2.5% of towers) before averaging adds only about 1% efficiency.
4. **The EMCal − HCal difference still grows slightly with pT:** 1.3 → 2.6 ns standard, −0.1 → 0.8 ns corrected. It is probably the oHCal ~1 GeV population above.
5. **iHCal has only 787 towers above 0.5 GeV** in 120k jets, and its corrected times look worse (RMS 3.1 → 5.8 ns). It is too sparse to judge, and it barely contributes to the jet time.

**Caveats:**
- These are MB + jet-triggered leading jets, not photon+jet.
- Window "efficiency" here is the in-window fraction of all jets, so any out-of-time background is in the denominator. It is not a signal efficiency.
- Even ±4 ns keeps only about 90% of jets: t0/MBD resolution, plus whatever background is present.

**Correction to the Sep 28 conclusion above:** the "edge jets are ≳ 90% signal because they are back-to-back" argument assumed background flat in Δφ. Out-of-time background can be back-to-back (the user confirmed this on Oct 1), so that argument does not hold.

## Without per-run constants (Oct 2)

**Test (`tcal_runindep.C`, writes `tcal_runindep.log`):**
- A run-independent correction δ(calorimeter, tower E) = median(t_corr − t_std) was derived on 10 runs and applied to the standard time on the other 10.
- It recovers only about a third of the gain. ±2 ns in-fraction across EM fraction: 0.47–0.80 standard, 0.55–0.81 with δ, 0.68–0.82 fully calibrated.
- So the per-run group alignment (phase 1) carries most of the improvement.

**Wide jet windows (`tcal_window_eff_wide.C`):**
- ±7 ns around the median of t_MBD − t_jet keeps 0.93–0.98 vs pT and 0.97–0.98 vs EM fraction.
- ±10 ns keeps 0.98–1.00.
- At these widths the standard and calibrated times give the same numbers.

**Recommendation (pending confirmation that the background is whole out-of-time photon+jet topologies):**
- Put the tight cut on the photon cluster time. It is EMCal-only and well behaved (Sep 28: RMS 0.9 ns).
- Use a wide jet window (about ±7–10 ns) on the standard time.
- Cross-check on Dading's 284 calibrated runs.
