#!/bin/bash
# Makes a JER-iteration work dir JER_iter/work/TAG: the GAMMAJET_UNFOLD tree that run_iter.sh points
# main's macros at. hists/ and reweight/ start as copies of main's (Jet12 purity inputs, the previous
# reweight); trees/ links to JER_iter/trees (v1jerz Photon trees + links to main's trees).
# usage: JER_iter/setup_work.sh TAG
set -eo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MAIN=$(dirname "$HERE")
TAG=${1:?tag}
WORK=$HERE/work/$TAG
[ -e "$WORK" ] && { echo "$WORK exists"; exit 1; }
mkdir -p "$WORK"/{hists,insitu/inputs,insitu/output,insitu/pdfs,reweight,pdfs,logs}
ln -s ../../trees "$WORK/trees"
cp "$MAIN"/hists/*.root "$WORK/hists/"
cp "$MAIN"/reweight/*.root "$WORK/reweight/"
echo "made $WORK"
