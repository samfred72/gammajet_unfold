# insitu_variations

Checks of the gamma+jet in-situ JES (Oct 7-8 2026), summarized in `slides/insitu_variations.pdf`.

- `threejet_dr.C`: photon - third-jet Delta R by ABCD region. At R = 0.2 the three-jet veto removes the
  non-isolated regions B/D about twice as often as A/C, partly through jets at Delta R = 0.2-0.4 (inside
  the iso_topo_04 cone, outside the treemaker's Delta R < R exclusion). `pdfs/threejet_dr.pdf`.
- `threejet_runaway.C`: the R = 0.2 threejet scan vs the assumed p_a (veto fractions, simple ABCD purity,
  the pipeline runs). With the 5 GeV veto the 25-35 GeV purity is low at any p_a and each scan lands below
  its start, so the in-situ loop walks to the 0.80 edge. `pdfs/threejet_runaway.pdf`.
  Outcome: threejet uses the nominal purity at R = 0.2 (`ana::getPurity*`, src/ana.cc).
- The jet-cut (7, 9 GeV), timing-cut and JER low-pT extrapolation numbers in the deck come from one-scan
  trials in separate worktrees; their logs are archived in `logs/trials/` (not in git). The JER study's
  macros and patch are in `../smear_test/` (`draw_gammajet_extrap*.C`, `unfolder_jerextrap.patch`).
