#include "/home/samson72/sphnx/gammajet_unfold/src/unfold_utility.h"
#include "RooUnfoldBayes.h"
using namespace std;

TH1D * unfold_utility::unflattenXj(TH1D * flat, int ipt, const char * name) {
  TH1D * h = new TH1D(name, ";x_{J#gamma};Counts", ana::nUnfoldXjBins, ana::unfoldXjBins);
  for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
    int flatbin = ipt*(ana::nUnfoldXjBins+2) + ixj + 1;
    h->SetBinContent(ixj+1, flat->GetBinContent(flatbin+1));
    h->SetBinError(ixj+1, flat->GetBinError(flatbin+1));
  }
  return h;
}

void unfold_utility::reflattenXj(TH1D * perPt, int ipt, TH1D * flatOut) {
  for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
    int flatbin = ipt*(ana::nUnfoldXjBins+2) + ixj + 1;
    flatOut->SetBinContent(flatbin+1, perPt->GetBinContent(ixj+1));
    flatOut->SetBinError(flatbin+1, perPt->GetBinError(ixj+1));
  }
}

TH1D * unfold_utility::unfoldOnce(TH1D * respRecoTemplate, TH1D * respTruthTemplate, TH2D * matrix, TH1D * flatMeasured, int niter, const char * name) {
  RooUnfoldResponse response(respRecoTemplate, respTruthTemplate, matrix);
  return unfoldOnce(&response, flatMeasured, niter, name);
}

TH1D * unfold_utility::unfoldOnce(RooUnfoldResponse * response, TH1D * flatMeasured, int niter, const char * name) {
  RooUnfoldBayes unfold(response, flatMeasured, niter, 0, 1);
  // -1 fully silences RooUnfoldBayes's per-call console spam ("Now unfolding...",
  // "Iteration : N", "Chi^2 of change ...", "Calculating covariances...", the priors
  // vector dump, etc.) - the default verbosity of 1 makes toy/iteration-scan loops here
  // print thousands of lines per macro run.
  unfold.SetVerbose(-1);
  return (TH1D*)((TH1D*)unfold.Hreco())->Clone(name);
}
