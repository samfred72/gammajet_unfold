#/bin/bash
# Usage: runall_puritymaker.sh
# Regenerates hists/purity_<systag>.root (all seven jet radii, one ana::rnames[ir]
# subdirectory each) for nominal and every systematic variation currently produced by
# runall_unfold.sh - purity is "of paired photons" (see puritymaker.C's header comment),
# and pairing genuinely differs by jet radius, so it's no longer just an R=0.4 number
# reused everywhere (see src/ana.h's purityFilename/getPurity ir parameter).
# One ROOT process per systag now - puritymaker.C itself loops over every radius
# internally and writes them all into the one file/subdirectory structure, so this no
# longer needs a separate process per radius.
#
# Keep SYSTAGS in sync with ana::systags (src/ana.h) - the definitive systag list - bash
# can't read a C++ static vector<string> directly, so this is a duplicated, explicit
# list (same convention as insitu/run_grid.sh's SYSTAGS).
SYSTAGS=(nominal JERhigh JERlow emscale_high emscale_low jes_high jes_low threejet narrowBDT narrowISO)

for SYSTAG in "${SYSTAGS[@]}"; do
  bash run_puritymaker.sh $SYSTAG &
done
wait
echo "All Done!"
