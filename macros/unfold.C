#include "/home/samson72/sphnx/gammajet_unfold/src/unfolder.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/object.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);
R__LOAD_LIBRARY(libRooUnfold.so);

void unfold(string trigger="Photon5", string sim="pythia", string systag="nominal") {

  bool dodraw = false;

  unfolder uf(trigger, sim, systag);
  uf.set_dodraw(dodraw);
  uf.fill_matrix();
  if (!dodraw) {
    uf.unfold();
    uf.end();
  }
}
