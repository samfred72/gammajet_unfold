// Quick diagnostic: print TH1::GetMean() of the raw (region-A, not purity-corrected or
// unfolded) Data xJ histogram for each of the ana::nPtBins pT bins.
//
// hrecoxj%i_0 (ir = nominal jet radius index) is the flattened (pT,xJ) measured histogram
// for region A, written by unfolder.cc - same source drawing/draw_final_result.C:67 reads
// via drawer::get(..., 0). unfold_utility::unflattenXj pulls out one pT slice's
// ana::nUnfoldXjBins real xJ bins (see src/unfold_utility.h:15-19).

#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/unfold_utility.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);

void print_data_xj_means() {
  const int ir = 2; // nominal jet radius index, R=0.4 (drawing/draw_final_result.C:38)

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
