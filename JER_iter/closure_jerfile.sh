#!/bin/bash
# Closure of the GAMMAJET_JER=<templates> re-smear on the v1jerz trees: Photon10 nominal unfold with the
# stored smeared branches (tree) and re-smeared from the ORIGINAL templates (file) must agree.
# Runs in its own work dir (JER_iter/work/closure_jerfile).
# usage: JER_iter/closure_jerfile.sh
set -eo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MAIN=$(dirname "$HERE")
WORK=$HERE/work/closure_jerfile
[ -d "$WORK" ] || bash "$HERE/setup_work.sh" closure_jerfile
export GAMMAJET_UNFOLD=$WORK GAMMAJET_NO_PROGRESS=1 LD_LIBRARY_PATH=$MAIN/src:$LD_LIBRARY_PATH
H=$WORK/hists/Photon10_pythia_nominal_unfolding.root
I=$WORK/insitu/inputs/Photon10_pythia_nominal_insitu.root
TMPL=$HERE/templates/jerband_smearing_templates.root
for mode in tree file; do
  if [[ $mode == tree ]]; then export GAMMAJET_JER=tree; else export GAMMAJET_JER=$TMPL; fi
  (cd "$MAIN/macros" && root -b -l -q 'unfold.C("Photon10","pythia","nominal")') > "$WORK/logs/unfold_$mode.log" 2>&1
  cp "$H" "$WORK/hists_$mode.root"
  cp "$I" "$WORK/insitu_$mode.root"
done
cd "$HERE" && root -b -l -q "compare_files.C(\"$WORK/hists_tree.root\",\"$WORK/hists_file.root\",\"$WORK/insitu_tree.root\",\"$WORK/insitu_file.root\")"
