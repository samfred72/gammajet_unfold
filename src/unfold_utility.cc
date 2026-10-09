#include "unfold_utility.h"
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

TH1D * unfold_utility::unfoldOnce(TH1D * respRecoTemplate, TH1D * respTruthTemplate, TH2D * matrix, TH1D * flatMeasured, int niter, const char * name, bool includeSystematics) {
  RooUnfoldResponse response(respRecoTemplate, respTruthTemplate, matrix);
  TH1D * result = unfoldOnce(&response, flatMeasured, niter, name, includeSystematics);
  // RooUnfoldResponseT's destructor only calls ClearCache(), not Reset(), so its histogram clones
  // leak once per call (fatal in toy loops). Reset() frees them.
  response.Reset();
  return result;
}

TH1D * unfold_utility::unfoldOnce(RooUnfoldResponse * response, TH1D * flatMeasured, int niter, const char * name, bool includeSystematics) {
  RooUnfoldBayes unfold(response, flatMeasured, niter, 0, 1);
  // Silence RooUnfoldBayes's per-call printout.
  unfold.SetVerbose(-1);
  // Response-matrix-statistics covariance term (needs the local RooUnfold patch, see CLAUDE.md).
  // ~3 s per call, so toy loops that only read bin contents pass includeSystematics=false.
  if (includeSystematics) unfold.IncludeSystematics(RooUnfolding::kAll);
  // Hunfold() returns a new histogram the caller owns: clone into `name` and free the original.
  TH1D * raw = (TH1D*)unfold.Hreco();
  TH1D * result = (TH1D*)raw->Clone(name);
  delete raw;
  // RooUnfoldT keeps its own deep copy of the response and its destructor frees nothing, so free
  // that copy's four histograms here (the base type has no public Reset()).
  RooUnfoldResponseT<TH1,TH2> * internalResponse = unfold.response();
  if (internalResponse) {
    delete internalResponse->Hmeasured();
    delete internalResponse->Hfakes();
    delete internalResponse->Htruth();
    delete internalResponse->Hresponse();
    delete internalResponse;
  }
  return result;
}

bool unfold_utility::purityCorrectCoeffs(float pA, float pC, float NA, float NC, float & coeffA, float & coeffC) {
  if (NC <= 0) { coeffA = 1; coeffC = 0; return false; }
  float K = NA/NC;
  if (fabs(pA - pC) < minPurityDiff) {
    coeffA = 1;
    coeffC = (1-pA)*K;
    return false;
  }
  coeffA = pA*(1-pC)/(pA-pC);
  coeffC = pA*(1-pA)*K/(pA-pC);
  return true;
}

TH1D * unfold_utility::purityCorrect(TH1D * A, TH1D * C, float pA, float pAErrLow, float pAErrHigh,
    float pC, float pCErrLow, float pCErrHigh, const char * name, TGraphAsymmErrors ** graphOut, bool quiet) {
  float NA = A->Integral();
  float NC = C->Integral();
  if (NC <= 0) {
    if (!quiet) cout << "WARNING: " << name << " has zero region-C statistics - cannot cross-normalize, skipping." << endl;
    if (graphOut) *graphOut = nullptr;
    return nullptr;
  }
  float K = NA/NC;
  float coeffA, coeffC;
  bool useExact = purityCorrectCoeffs(pA, pC, NA, NC, coeffA, coeffC);
  if (!useExact && !quiet) {
    cout << "WARNING: " << name << " has |P_A-P_C| = " << fabs(pA-pC) << " < " << minPurityDiff
         << " - falling back to the single-purity (region-C-is-background) formula for this bin." << endl;
  }
  TH1D * h = (TH1D*)A->Clone(name);
  TGraphAsymmErrors * g = graphOut ? new TGraphAsymmErrors(A->GetNbinsX()) : nullptr;
  if (g) g->SetName(Form("%s_graph", name));
  for (int i = 1; i <= A->GetNbinsX(); i++) {
    float a  = A->GetBinContent(i);
    float ae = A->GetBinError(i);
    float c  = C->GetBinContent(i);
    float ce = C->GetBinError(i);

    float content = coeffA*a - coeffC*c;
    float errLow, errHigh;
    if (!useExact) {
      float dLow  = K*c*pAErrLow;
      float dHigh = K*c*pAErrHigh;
      errLow  = sqrt(ae*ae + pow(coeffC*ce,2) + pow(dLow,2));
      errHigh = sqrt(ae*ae + pow(coeffC*ce,2) + pow(dHigh,2));
    }
    else {
      auto signalVal = [&](float pa, float pc) {
        return pa*((1-pc)*a - (1-pa)*K*c) / (pa-pc);
      };
      float statErr = sqrt(pow(coeffA*ae,2) + pow(coeffC*ce,2));
      float devPAHigh = signalVal(pA+pAErrHigh, pC) - content;
      float devPALow  = signalVal(pA-pAErrLow,  pC) - content;
      float devPCHigh = signalVal(pA, pC+pCErrHigh) - content;
      float devPCLow  = signalVal(pA, pC-pCErrLow)  - content;
      float upSum2 = 0, downSum2 = 0;
      for (float dev : {devPAHigh, devPALow, devPCHigh, devPCLow}) {
        if (dev > 0) upSum2   += dev*dev;
        else         downSum2 += dev*dev;
      }
      errHigh = sqrt(statErr*statErr + upSum2);
      errLow  = sqrt(statErr*statErr + downSum2);
    }

    h->SetBinContent(i, content);
    h->SetBinError(i, std::max(errLow, errHigh));
    if (g) {
      double xc  = A->GetXaxis()->GetBinCenter(i);
      double xlo = xc - A->GetXaxis()->GetBinLowEdge(i);
      double xhi = A->GetXaxis()->GetBinUpEdge(i) - xc;
      g->SetPoint(i-1, xc, content);
      g->SetPointError(i-1, xlo, xhi, errLow, errHigh);
    }
  }
  if (graphOut) *graphOut = g;
  return h;
}

TH1D * unfold_utility::purityCorrectBkg(TH1D * A, TH1D * C, float pA, float pC, const char * name) {
  float NA = A->Integral();
  float NC = C->Integral();
  if (NC <= 0) return nullptr;
  float coeffA, coeffC;
  purityCorrectCoeffs(pA, pC, NA, NC, coeffA, coeffC);
  TH1D * h = (TH1D*)A->Clone(name);
  for (int i = 1; i <= A->GetNbinsX(); i++) {
    float a  = A->GetBinContent(i);
    float ae = A->GetBinError(i);
    float c  = C->GetBinContent(i);
    float ce = C->GetBinError(i);
    float bkg = (1-coeffA)*a + coeffC*c;
    float err = sqrt(pow((1-coeffA)*ae,2) + pow(coeffC*ce,2));
    h->SetBinContent(i, bkg);
    h->SetBinError(i, err);
  }
  return h;
}

TH1D * unfold_utility::densityForDisplay(TH1D * h, const char * name) {
  TH1D * hd = (TH1D*)h->Clone(name);
  hd->Scale(1., "width");
  hd->GetYaxis()->SetTitle("Counts / bin width");
  return hd;
}

TH1D * unfold_utility::buildFullyCorrected(TH1D * flatA, TH1D * flatC, const char * tag, string systag, int ir, bool quietAll) {
  TH1D * flatCorrected = (TH1D*)flatA->Clone(Form("hxjcorrected_flat_%s", tag));
  flatCorrected->Reset("ICES");
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    float ptlow  = ana::ptBins[ipt];
    float pthigh = ana::ptBins[ipt+1];
    float pA        = ana::getPurity(ptlow, pthigh, systag, ir);
    float pAErrLow  = ana::getPurityErrorLow(ptlow, pthigh, systag, ir);
    float pAErrHigh = ana::getPurityErrorHigh(ptlow, pthigh, systag, ir);
    float pC        = ana::getPurityC(ptlow, pthigh, systag, ir);
    float pCErrLow  = ana::getPurityCErrorLow(ptlow, pthigh, systag, ir);
    float pCErrHigh = ana::getPurityCErrorHigh(ptlow, pthigh, systag, ir);
    TH1D * A = unflattenXj(flatA, ipt, Form("htmpA_%s_pt%d", tag, ipt));
    TH1D * C = unflattenXj(flatC, ipt, Form("htmpC_%s_pt%d", tag, ipt));
    // Last pT bin (35-100 GeV migration buffer) has an empty region B.
    bool quiet = quietAll || (ipt == ana::nPtBins - 1);
    TH1D * hcorr = purityCorrect(A, C, pA, pAErrLow, pAErrHigh, pC, pCErrLow, pCErrHigh, Form("htmpcorr_%s_pt%d", tag, ipt), nullptr, quiet);
    // Buffer pT bins (not reported): negative purity-corrected contents set to 0, errors kept. Bayesian
    // unfolding assumes a non-negative input; a negative buffer bin (35-100 GeV, x_J 0.4-0.5: -4 +- 4)
    // made the unfolding blow up at isolated iterations beyond ~7 (refolding and response-toy chi2
    // spikes; claude_checks/unfold_iterations). Reported bins are left as measured.
    bool isBuffer = ipt < ana::firstUsedPtBin || ipt >= ana::firstUsedPtBin + ana::nPtBinsUsed;
    if (hcorr && isBuffer && clampBufferNegatives) {
      for (int b = 1; b <= hcorr->GetNbinsX(); b++) if (hcorr->GetBinContent(b) < 0) hcorr->SetBinContent(b, 0);
    }
    reflattenXj(hcorr ? hcorr : A, ipt, flatCorrected);
    delete A; delete C; if (hcorr) delete hcorr;
  }
  return flatCorrected;
}
