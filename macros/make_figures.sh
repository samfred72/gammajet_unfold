#!/bin/bash
# Redraw every figure from the existing ROOT outputs (no calculation), with note or paper labels.
# Usage: bash macros/make_figures.sh [--paper] [--status WORD]
#   default   note labels: "sPHENIX Internal", "p+p Run24 Data", "Pythia8 #gamma+jet MC"
#   --paper   publication labels: "sPHENIX", "p+p #sqrt{s}=200 GeV", "Pythia8 #gamma+jet"
#   --status  word after "sPHENIX" (e.g. Preliminary, Internal; "" for none)
# The PDFs overwrite pdfs/ and insitu/pdfs/ in place; copy them out before switching styles.
set -eo pipefail
: "${GAMMAJET_UNFOLD:?set GAMMAJET_UNFOLD to the gammajet_unfold checkout (see README.md)}"
unset GJ_PAPER_LABELS GJ_LABEL_STATUS
while [[ $# -gt 0 ]]; do
  case $1 in
    --paper) export GJ_PAPER_LABELS=1; shift ;;
    --status) export GJ_LABEL_STATUS="$2"; shift 2 ;;
    *) echo "Usage: $0 [--paper] [--status WORD]" >&2; exit 1 ;;
  esac
done
source "$GAMMAJET_UNFOLD/macros/pipeline_steps.sh"
ps_figures
