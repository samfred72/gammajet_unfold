# claude_checks

Responses to the PPG18 independent-check deck (J. Nagle, "Checks of the PPG18 p+p γ-jet
xJγ analysis – round 2", Sep 10 2026, `ppg18_checks.pdf`). One subdirectory per issue.
Macros run interpreted (`root -l -b -q macro.C`, no ACLiC) and load `libgammajet_unfold.so`
explicitly.

| # | Issue | Directory | Status |
|---|-------|-----------|--------|
| 1 | Jet timing window 0 < t_MBD − t_jet < 4 ns (data only) | `jet_timing/` | **Resolved (note).** Edge jets are in-time recoil jets; for the reprocess, keep leading-passing selection with a wide (offset-corrected) window, store the no-window leading jet for the efficiency figure. R dependence: HCal-rich jet time shifts with R (time-walk); loss same for R=0.2–0.4, grows for R≥0.5. |
| 3 | In-situ scans at the 0.900 scan edge; Fig. 10; Appendix 11 | `insitu_scan_edge/` | Jet12_long needs a reprocess (decided Sep 28) before Fig. 10 / Jet12 scan are redone. Jet12 signal/bkg xJ: region-C background model OK to 0.01–0.03, bias pushes p_a up (wrong way for the edge); Jet12_long found to be ~77% truth-tagged in region A. Purity chain closes on Jet12 at R=0.2/0.3/0.4 (≤0.7% at 15–20 GeV). Data: purity R-independent; corrected data/MC ≈ 0.93 at R=0.2 (15–25 GeV) like other R — edge likely from scan response (xJ floor), untested. |
| 4 | Cluster-pT reweight follows the fake photons | `reweight/` (repo) | **Done (Sep 28):** `reweight/make_vz_pt_reweight.C` now uses the purity-corrected region-A Data spectrum vs truth-matched MC; new `vz_pt_reweight.root` written: fit to all 1-GeV bins: slope (−11±7)×10⁻³/GeV, weight 1.00 (15) → 0.80 (35 GeV), consistent with flat at 1.7σ; merged-bin fit gives (−0.2±6.0)×10⁻³ (old weight was 0.79→0.29). Plots split into pdfs/reweight_vz.pdf and pdfs/reweight_pt.pdf. Needs the pipeline rerun (user) to propagate. Note §2.1 + Fig. caption updated. |
| – | Data vs Pythia purity at 15–20 GeV (0.61 vs 0.74) | `purity/` | Diagnosed: gap is entirely B/D (data fakes pass tight BDT ~1.5× more); not signal leakage; data fakes slightly more isolated. BDT-score comparison pending. |
| 5 | Three-jet veto systematic / double count | code (`src/`, `insitu/`) | **Fixed (Sep 28):** each systag with its own in-situ scan now unfolds Data with its own p_a (`ana::jesBySystag`, `jesForSystag`); jes_high/low use only the nominal fit's stat error (`jesStatErrLow/High`); `draw_jes_summary.C` regenerates both; lib rebuilt clean. Note §5.4, §6 intro, §6.2 updated. Needs pipeline rerun; scan-edge values (R=0.2, threejet) must be fixed first. |
| 6 | Purity uncertainty inside the stat errors | – | **Ignored** (decided Sep 28). |
| 7 | Photon smearing 2 / 0 / 6% | – | **Done (Sep 28):** pipeline already uses the PPG12 EM-resolution prescription (`ana::emResolutionSigma`); note §6.4 rewritten with the parameterization and variations. |
| 8 | EM scale ±1.1% (code) vs ±1.5% (note) | – | **Done (Sep 28):** pipeline uses ±1.48% (PPG12 / Calorimeter Calibration WG, `ana::emscaleShift`); note §2, §3.2 threshold sentence and §6.3 updated. |
| 9 | 35–100 GeV buffer bin has P_A = 0 | – | **Ignored** (decided Sep 28). |
| 10 | Particle-level / fiducial definition; purity vs \|η\| | `purity/` | Purity vs \|η\| done: falls 0.72→0.58 with \|η\| (15–20 GeV); \|η\|-integrated ABCD underestimates P_A by ~0.03. Decision: \|η^γ\| < 0.7 (as PPG12, confirmed from its note). **Note changes deferred until the reprocess** (incl. the particle-level definition text). |
| 11 | Note-vs-code list | – | **Done in note (Sep 28):** trigger (no bit required; calojet stream), \|t_MBD\| < 20 ns (not applied; none beyond), random seed (default-seeded member generator), purity values, EM scale (issue 8), Fig. 10 (page → R=0.4, honest edge statement), Fig. 11 text (p_T^γ > 10 GeV, A/C background means). Particle-level definition deferred with issue 10. |
| 12 | Reproducibility (RooUnfold patches, stale scans, archiving) | – | **RooUnfold patches documented** in a new note appendix (Sep 28). Stale Jet12/shape-χ² scans: covered by the Jet12_long reprocess (issue 3). |

## Decisions for the next production

- **Photon acceptance |η^γ| < 0.7** (decided Sep 28, issue 10; **implemented Sep 28**: `ana::photonEtaCut = 0.7` used by `unfolder::check_pair` for reco and truth, jets still `|η_jet| < ana::etacut − R` = 1.1 − R; reweight macro follows it). Apply it to the data and
  MC reco selection, the purity (ABCD counts and Photon-MC leakage templates), the
  in-situ JES inputs, and the particle-level truth definition (currently
  |η^γ| < 1.1), and update the note's fiducial definition to match.
  - Why: at |η| > 0.7 the isolation cone leaves the EMCal acceptance. The purity then
    depends on |η| (0.72 → 0.58 at 15–20 GeV), and the |η|-integrated ABCD
    underestimates it by ~0.03.
  - Consequences: the reco xJ shape is unchanged within statistics; the in-situ p_a
    rises by +0.004–0.018 and R = 0.2 returns inside the scan range (0.903).
  - Cost: ~23% of region-A pairs; stat error +3% at 15–20 GeV, +10–14% above 20 GeV.
  - Details: `purity/README.md`.

- **Jet/cluster timing** (agreed Sep 28, issue 1). Keep "leading object passing the
  window", but use a wide window, ideally with the EMCal/HCal time offset corrected. Store
  the leading jet with no timing requirement for the efficiency figure. Details:
  `jet_timing/README.md` ("Resolution").
- **Revisit radius-specific JER smearing** at the reprocess (issue 3): one
  `h_jerband_quaddiff` is currently used for every R. Details:
  `insitu_scan_edge/README.md` ("Open for this issue").
- **Reprocess Jet12_long** (decided Sep 28, issue 3). Needed before the Jet12-reference
  scan and Fig. 10 are redone. Also check its truth-photon content: ~77% of its region-A
  pairs are truth-tagged photons (`insitu_scan_edge/README.md`).
- **In-situ scan range 0.80–1.00** (implemented Sep 28, issue 3): `insitu_utility::scanLow = 0.80`, `scanN = 2000` (1e-4 step); the unfolder's in-situ pairing floor follows `scanLow`.
- **Three-jet veto definition** (diagnosed Sep 29, `threejet/`; **implemented Oct 1–2**, trees not yet reprocessed). The veto had removed ~78% of selected Pythia events vs ~35% of data, because MC soft jets passed the threshold only via the OR over six smeared pTs.
  - **Treemaking (CaloAna):**
    - The third jet is now the highest-pT in-time jet outside the photon cone, other than the leading recoil jet.
    - Only one third jet is stored: `thirdjet_pt/eta/phi[7]`.
    - `thirdjet_pt` is calibrated pT in data and the nominal smeared pT in MC.
    - It is counted after the timing cut.
  - **Unfolder:** the event is vetoed if the third jet's pT > `ana::thirdJetPtCut` = 5 GeV, in the recoil jet's pT definition (data in-situ corrected).
  - **Scope:** any η, any ΔR. Reco only: there is no truth veto, so the observable is unchanged.
  - **Decided against:** jet-quality cuts, fake-rate tuning, storing all jets.
  Plot of the old behaviour: `threejet/pdfs/draw_thirdjet_dr.pdf`.
- **Full pipeline script** `macros/run_full_pipeline.sh`: unfold → purity → reweight → [in-situ scan → recompile → unfold → purity] until the p_a table converges → in-situ note plots → drawing.
- **Later (user, Sep 28):** regenerate the note's purity_check plots (Fig. 11 region A/C xJ, ABCD-validation appendix) with |η^γ| < 0.7 in `/home/samson72/sphnx/purity_check`, then `sync_figs.sh` and update the Fig. 11 numbers/caption.
