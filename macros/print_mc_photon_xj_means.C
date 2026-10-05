// Print the mean of the reco-level MC photon xJ (Photon5/10/20 combined) per pT bin.

#include "../src/drawer.h"
#include "../src/unfold_utility.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);

void print_mc_photon_xj_means() {
  const int ir = 2; // nominal R=0.4

  drawer d("pythia", "nominal");
  TH1D * flatMC = d.get(Form("hrecoxj%i", ir), 1);
  if (!flatMC) {
    cout << "ERROR: couldn't read combined MC photon hrecoxj" << ir << endl;
    return;
  }

  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    TH1D * h = unfold_utility::unflattenXj(flatMC, ipt, Form("hMCPhotonXj_pt%d", ipt));
    cout << "pT bin " << ipt << " [" << ana::ptBins[ipt] << "," << ana::ptBins[ipt+1]
         << "): mean x_{J} = " << h->GetMean() << endl;
    delete h;
  }
}
