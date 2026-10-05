// Print the mean of the raw region-A Data xJ per pT bin.

#include "../src/drawer.h"
#include "../src/unfold_utility.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);

void print_data_xj_means() {
  const int ir = 2; // nominal R=0.4

  drawer d("pythia", "nominal");
  TH1D * flatA = d.get(Form("hrecoxj%i_0", ir), 0);
  if (!flatA) {
    cout << "ERROR: couldn't read hrecoxj" << ir << "_0 from Data_nominal_unfolding.root" << endl;
    return;
  }

  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    TH1D * h = unfold_utility::unflattenXj(flatA, ipt, Form("hDataXj_pt%d", ipt));
    cout << "pT bin " << ipt << " [" << ana::ptBins[ipt] << "," << ana::ptBins[ipt+1]
         << "): mean x_{J} = " << h->GetMean() << endl;
    delete h;
  }
}
