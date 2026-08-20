root -b -l -q "grid_insitu.C(\"nominal\",2)"
# Non-purity-corrected cross-check: Data Region A (raw) vs. Pythia8 Jet12_full Region A
# (raw) - requires insitu/Jet12_full_pythia_nominal_insitu.root, produced by running
# unfold.C("Jet12_full","pythia","nominal") first (see run_unfold.sh / runall_unfold.sh;
# Jet12_full reads from purity_check/gammajet_pythia_Jet12_full.root - see
# src/treeuser.h's trigger=="Jet12_full" special case).
root -b -l -q "grid_insitu_jet12.C(\"nominal\",2)"
# Unfolded extension of grid_insitu.C's purity-corrected branch: unfolds the
# purity-corrected Data x_J through the nominal photon+jet response matrix and finds the
# JES scale that matches the UNFOLDED mean to the TRUTH-level MC mean, instead of two
# reco-level means - requires hists/Data_nominal_unfolding.root and
# hists/Photon{5,10,20}_pythia_nominal_unfolding.root (see macros/run_unfold.sh /
# runall_unfold.sh) to already exist, since it reuses their stored response matrix.
root -b -l -q "grid_insitu_unfolded.C(\"nominal\")"
