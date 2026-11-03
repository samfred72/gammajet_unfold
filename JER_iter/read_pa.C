// Prints "PA p_a errLow errHigh" per jet radius from insitu/output/grid_insitu_<systag>.root.
#include "../src/ana.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);

void read_pa(string systag = "nominal") {
  TFile f(ana::path("insitu/output/grid_insitu_" + systag + ".root"));
  for (int ir = 0; ir < ana::nJetR; ir++) {
    TTree * t = (TTree*)f.Get(Form("%s/results", ana::rnames[ir]));
    float pa, lo, hi;
    t->SetBranchAddress("pa_puritycorrected", &pa);
    t->SetBranchAddress("errLow_puritycorrected", &lo);
    t->SetBranchAddress("errHigh_puritycorrected", &hi);
    t->GetEntry(0);
    printf("PA %.4f %.4f %.4f\n", pa, lo, hi);
  }
}
