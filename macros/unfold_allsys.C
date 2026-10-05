#include "../src/unfolder.h"
#include "../src/object.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);
R__LOAD_LIBRARY(libRooUnfold.so);

// unfold.C for many systags in one tree read (default ana::systags); same output files.
void unfold_allsys(string trigger="Photon5", string sim="pythia",
    vector<string> systags=ana::systags) {

  bool dodraw = false;

  unfolder uf(trigger, sim, systags);
  uf.set_dodraw(dodraw);
  // GAMMAJET_NO_PROGRESS=1 turns off the per-event progress line.
  if (gSystem->Getenv("GAMMAJET_NO_PROGRESS")) uf.set_progress(false);
  uf.fill_matrix();
  if (!dodraw) {
    uf.unfold();
    uf.end();
  }
}
