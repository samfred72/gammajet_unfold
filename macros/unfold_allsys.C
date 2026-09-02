#include "/home/samson72/sphnx/gammajet_unfold/src/unfolder.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/object.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);
R__LOAD_LIBRARY(libRooUnfold.so);

// One-pass, all-systematics version of unfold.C: reads the raw tree once and fills
// every systag's histograms per event (see unfolder.h's constructor comment / the
// systag loop inside unfolder::fill_matrix()), instead of re-running the whole analysis
// once per systag (each re-reading and re-deserializing the same tree from disk).
// Produces exactly the same output files (hists/<trigger>_<sim>_<systag>_unfolding.root,
// insitu/<trigger>_<sim>_<systag>_insitu.root) as calling unfold.C once per systag would.
//
// systags defaults to ana::systags (src/ana.h) - the definitive nominal + 14-named-
// systematic-variation list every consumer of systags reads, so adding a new one there
// propagates here automatically; pass a shorter list (e.g. just {"nominal","JERhigh"})
// to validate against unfold.C's single-systag output for a subset without paying for
// the full list.
void unfold_allsys(string trigger="Photon5", string sim="pythia",
    vector<string> systags=ana::systags) {

  bool dodraw = false;

  unfolder uf(trigger, sim, systags);
  uf.set_dodraw(dodraw);
  uf.fill_matrix();
  if (!dodraw) {
    uf.unfold();
    uf.end();
  }
}
