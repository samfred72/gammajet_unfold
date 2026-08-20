#/bin/bash
# Usage: runall_puritymaker.sh
# Regenerates hists/purity_<systag>.root for nominal and every systematic variation
# currently produced by runall_unfold.sh.
# Each systag needs its own fresh ROOT process: combine_hists() creates fixed-name
# objects (e.g. "combined", "func", "bootstrap0"..) that would collide if puritymaker()
# were called repeatedly within one session.
for SYSTAG in nominal JERhigh JERlow emscale_high emscale_low jes_high jes_low threejet narrowBDT narrowISO; do
  bash run_puritymaker.sh $SYSTAG &
done
wait
echo "All Done!"
