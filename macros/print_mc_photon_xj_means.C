// Quick diagnostic: print TH1::GetMean() of the reco-level MC photon xJ histogram (all
// photon samples - Photon5/10/20 - scale-combined, see drawer::combineMC) for each of the
// ana::nPtBins pT bins.
//
// hrecoxj%i (ir = nominal jet radius index, no _0/_2 ABCD-region suffix - that split is
// Data-only) is the same flattened (pT,xJ) reco template drawing/draw_final_result.C:64
// uses as respRecoTemplate. unfold_utility::unflattenXj pulls out one pT slice's
// ana::nUnfoldXjBins real xJ bins (see src/unfold_utility.h:15-19) - mirrors
// print_data_xj_means.C but for MC photon (type=1) instead of Data (type=0).

#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/unfold_utility.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);

void print_mc_photon_xj_means() {
  const int ir = 2; // nominal jet radius index, R=0.4 (drawing/draw_final_result.C:38)

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
