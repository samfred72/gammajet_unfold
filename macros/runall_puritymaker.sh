#!/bin/bash
# Usage: runall_puritymaker.sh
# hists/purity_<systag>.root for every systag (all radii per file).
# Keep SYSTAGS in sync with ana::systags (src/ana.h).
SYSTAGS=(nominal JERhigh JERlow emscale_high emscale_low EMRhigh EMRlow jes_high jes_low threejet narrowBDT narrowISO narrowBDTbkg narrowISObkg wideISObkg timingwide)

for SYSTAG in "${SYSTAGS[@]}"; do
  bash run_puritymaker.sh $SYSTAG &
done
wait
echo "All Done!"
