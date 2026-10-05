#include "../src/unfolder.h"
#include "../src/object.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);
R__LOAD_LIBRARY(libRooUnfold.so);

void unfold(string trigger="Photon5", string sim="pythia", string systag="nominal") {

  bool dodraw = false;

  unfolder uf(trigger, sim, systag);
  uf.set_dodraw(dodraw);
  // GAMMAJET_NO_PROGRESS=1 (set by macros/run_full_pipeline.sh) turns off the per-event
  // progress line for batch runs.
  if (gSystem->Getenv("GAMMAJET_NO_PROGRESS")) uf.set_progress(false);
  uf.fill_matrix();
  if (!dodraw) {
    uf.unfold();
    uf.end();
  }
}
