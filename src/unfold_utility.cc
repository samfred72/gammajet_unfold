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
  // RooUnfoldResponseT's destructor (RooUnfold/src/RooUnfoldResponse.cxx) only calls
  // ClearCache(), NOT Reset() - so the _mes/_fak/_tru/_res histogram clones its own
  // Setup() makes (including a full clone of `matrix`, the largest of the four) are
  // never freed when `response` goes out of scope below. Reset() is what actually
  // deletes them; call it explicitly here since the destructor won't. Confirmed by
  // reading Reset() (which does `delete _mes/_tru/_res/_fak`) vs. ~RooUnfoldResponseT()
  // side by side - this is a real leak in that RooUnfold build, not a false lead.
  // Without this, every toy iteration in a bootstrap loop (draw_toy_vs_analytic.C,
  // toy_iterations.C ("data"), toy_iterations.C ("resp") - anything calling this overload in a
  // loop) permanently leaks one response-matrix clone; at 10k toys that's what exhausted
  // memory and forced a machine restart, where 1k toys had stayed under the radar.
  response.Reset();
  return result;
}

TH1D * unfold_utility::unfoldOnce(RooUnfoldResponse * response, TH1D * flatMeasured, int niter, const char * name, bool includeSystematics) {
  RooUnfoldBayes unfold(response, flatMeasured, niter, 0, 1);
  // -1 fully silences RooUnfoldBayes's per-call console spam ("Now unfolding...",
  // "Iteration : N", "Chi^2 of change ...", "Calculating covariances...", the priors
  // vector dump, etc.) - the default verbosity of 1 makes toy/iteration-scan loops here
  // print thousands of lines per macro run.
  unfold.SetVerbose(-1);
  // Include the response-matrix-statistics term in the analytic covariance, not just the
  // data-measurement term RooUnfoldBayes computes by default (_dosys=kNoSystematics) -
  // see CLAUDE.md's "Local RooUnfold Patch" section: that term used to silently poison
  // the covariance into NaN for any response with fakes (this project's always does, via
  // handleFakes=true just above) because of a real bug in RooUnfoldBayes::getCovariance()
  // reading one column past Eresponse()'s real end - patched locally in
  // /home/samson72/RooUnfold/src/RooUnfoldBayes.cxx, so this is now safe to enable.
  // Expensive (~3s/call measured on this project's binning, vs a few ms without it) - see
  // includeSystematics's header comment. Only enabled when the caller actually wants the
  // resulting GetBinError(); toy bootstrap loops (draw_toy_vs_analytic.C,
  // toy_iterations.C ("data"), toy_iterations.C ("resp")) should pass includeSystematics=false -
  // they only ever read GetBinContent() per toy, so paying for this on every one of
  // thousands of toy calls was pure waste (and looked like the macro had frozen: the
  // per-toy cost went from milliseconds to ~3s, but the progress printout only fires
  // every 200 toys).
  if (includeSystematics) unfold.IncludeSystematics(RooUnfolding::kAll);
  // Hunfold()/Hreco() (RooUnfold/src/RooUnfold.cxx) allocates a brand-new histogram on
  // every call via createHist(...) that RooUnfoldBayes keeps no internal reference to -
  // it's the caller's to free. Clone it into `name` for the caller, then delete the
  // original here so it isn't silently leaked once per unfoldOnce() call.
  TH1D * raw = (TH1D*)unfold.Hreco();
  TH1D * result = (TH1D*)raw->Clone(name);
  delete raw;
  // RooUnfoldT::SetResponse() (RooUnfold/src/RooUnfold.cxx), called internally by the
  // RooUnfoldBayes constructor above since ownership isn't transferred, makes its OWN
  // deep copy of `response` ("_res = new RooUnfoldResponseT<Hist,Hist2D>(*res)") - a
  // second, independent full clone of the response matrix, private to this `unfold`
  // object, on top of the caller's own `response` (already cleaned up in the other
  // overload above). RooUnfoldT's destructor (RooUnfold.cxx) is completely empty - it
  // frees none of _res/_meas/_truth/_bkg - so this internal copy leaks every single call
  // too, independent of the outer response's leak. This is the dominant leak: it fires
  // on every unfoldOnce() call (both toy loops), growing resident memory by a fixed
  // amount each time, which is what produces a smooth iteration-by-iteration slowdown
  // rather than a one-time step change. Reach in and free it the same way as above.
  // RooUnfoldResponseT<TH1,TH2> has no public Reset() (that's only declared on the
  // derived, non-template RooUnfoldResponse - see RooUnfoldResponse.h - and `response`'s
  // internal copy above is constructed as the plain base type), so free the same four
  // histograms Reset() would via their public non-const accessors instead.
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
    // Unlike the |P_A-P_C|<minPurityDiff fallback below, this used to print
    // unconditionally regardless of `quiet` - harmless per-call (buildFullyCorrected
    // falls back to raw region A), but a per-toy bootstrap loop hitting this on a
    // known-empty pT slice (e.g. the last-bin migration buffer) then prints it once per
    // toy, thousands of times, same as the other warning quiet already covers.
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
    // Last bin (ipt==nPtBins-1, the 35-100 GeV migration-only buffer) has an empty
    // region B in this data sample - see purityCorrect's quiet parameter comment.
    bool quiet = quietAll || (ipt == ana::nPtBins - 1);
    TH1D * hcorr = purityCorrect(A, C, pA, pAErrLow, pAErrHigh, pC, pCErrLow, pCErrHigh, Form("htmpcorr_%s_pt%d", tag, ipt), nullptr, quiet);
    reflattenXj(hcorr ? hcorr : A, ipt, flatCorrected);
    delete A; delete C; if (hcorr) delete hcorr;
  }
  return flatCorrected;
}
