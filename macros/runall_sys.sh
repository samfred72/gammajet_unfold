#!/bin/bash
# One-pass, all-systematics version: replaces the old approach of calling
# runall_unfold.sh once per systag (ten separate full tree reads per trigger) with a
# single call that reads each trigger's tree once and fills every systag's histograms
# in that one pass - see unfold_allsys.C / runall_unfold_allsys.sh. Validated bin-for-bin
# identical to the old per-systag-loop output.
bash runall_unfold_allsys.sh
