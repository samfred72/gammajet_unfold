#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/unfold_utility.h"
#include "../src/systematics_config.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Systematic uncertainties: Data's purity-corrected xJ is unfolded through each variation's own
// response and purity and compared with nominal (both shape-normalized, so a change in the
// accepted event count is not read as a shape change). herwig keeps Data at nominal and swaps
// the response generator; niterLow/High change only the iteration count; priorSensitivity is
// read from compute_prior_sensitivity.C's output (run that first, same radius).
//
// Total (per bin):
//   - JER/JES/emscale/EMR high/low: per-source, per-bin sign split (asymmetricSystematics).
//   - Purity: the five ABCD boundary sources in one symmetric quadrature sum (PPG12 Sec. 5.3).
//   - Unfolding: up = niterHigh (+) prior, down = niterLow (+) prior.
//   - threejet, timingwide, herwig: symmetric (timingwide is one-sided: a wider window, no direction).
// Writes the ratio and fractional difference per (source, pT bin), the combined Purity, Unfolding and
// total uncertainties, and the display-normalized nominal and variations to systematics.root.
// Calculation only: drawing/draw_systematics.C draws the pages. Sources: src/systematics_config.h.

// Jet radius index, set by compute_systematics(int); default R=0.4.
int ir = 2;
const int nPtBinsUsed = ana::nPtBinsUsed; // reported bins start at ana::firstUsedPtBin
// R=0.4 keeps the unsuffixed filenames; other radii get _<rname>.
string rootPathStr;
const char * rootPath;

// Data's purity-corrected xJ unfolded through systag's response (sim's generator) with niter
// iterations. Purity always comes from dataSystag.
TH1D * getUnfoldedData(string sim, string systag, int niter, const char * name) {
  drawer d(sim, systag);
  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i",ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i",ir), 1);
  TH1D * flatA = d.get(Form("hrecoxj%i_0",ir), 0);
  TH1D * flatC = d.get(Form("hrecoxj%i_2",ir), 0);
  TH1D * flatCorrected = unfold_utility::buildFullyCorrected(flatA, flatC, systag.c_str(), systag, ir);
  return unfold_utility::unfoldOnce(respRecoTemplate, respTruthTemplate, respMatrix2D, flatCorrected, niter, name);
}

void compute_systematics(int jetRadiusIndex = 2) {
  ir = jetRadiusIndex;
  rootPathStr = (ir == 2) ? ana::path("hists/systematics.root")
                          : Form("%s/hists/systematics_%s.root", ana::dir(), ana::rnames[ir]);
  rootPath = rootPathStr.c_str();

  TFile * fout = TFile::Open(rootPath, "RECREATE");

  TH1D * flatNominal = getUnfoldedData("pythia", "nominal", niterate, "hUnfoldedNominal");
  vector<TH1D*> nominalDisp(ana::nPtBins);
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hNom = unfold_utility::unflattenXj(flatNominal, ipt, Form("hNominal_pt%d", ipt));
    nominalDisp[ipt] = unfold_utility::densityForDisplay(hNom, Form("hNominalDisp_pt%d", ipt));
    nominalDisp[ipt]->Scale(1./nominalDisp[ipt]->Integral());
    nominalDisp[ipt]->GetYaxis()->SetTitle("Shape-normalized counts / bin width");
    fout->cd();
    nominalDisp[ipt]->Write();
  }

  // fracDiff[systag][ipt]: (variation-nominal)/nominal vs xJ.
  map<string, vector<TH1D*>> fracDiff;

  // Radius-specific (ir==2 keeps the unsuffixed name).
  string priorSensPathStr = (ir == 2) ? ana::path("hists/prior_sensitivity_nominal.root")
                                      : Form("%s/hists/prior_sensitivity_nominal_%s.root", ana::dir(), ana::rnames[ir]);
  TFile * fPriorSens = TFile::Open(priorSensPathStr.c_str());
  if (!fPriorSens || fPriorSens->IsZombie())
    cout << "WARNING: couldn't open " << priorSensPathStr << " - run "
         << "macros/compute_prior_sensitivity.C(\"nominal\", " << ir << ") first, or the "
         << "priorSensitivity source below will be empty." << endl;

  for (const string & systag : systematics) {
    bool isPriorSens = (systag == "priorSensitivity");
    const SystSource * src = isPriorSens ? nullptr : &systSources.at(systag);
    TH1D * flatVar = isPriorSens ? nullptr : getUnfoldedData(src->sim, src->dataSystag, src->niter, Form("hUnfolded_%s", systag.c_str()));
    fracDiff[systag] = vector<TH1D*>(ana::nPtBins);

    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      TH1D * hVar;
      if (isPriorSens) {
        // hw_pt<N>: the raw unfolded "w"-prior result, run through the same pipeline as every other
        // source. Not hPriorSensFracDiff_pt<N>, which is already a ratio-1.
        TH1D * hRaw = fPriorSens ? (TH1D*)fPriorSens->Get(Form("hw_pt%d", ipt)) : nullptr;
        if (!hRaw) {
          cout << "WARNING: hw_pt" << ipt << " missing from prior_sensitivity_nominal.root - skipping pt" << ipt << "." << endl;
          continue;
        }
        hVar = (TH1D*)hRaw->Clone(Form("h%s_pt%d", systag.c_str(), ipt));
      } else {
        hVar = unfold_utility::unflattenXj(flatVar, ipt, Form("h%s_pt%d", systag.c_str(), ipt));
      }
      TH1D * hVarDisp = unfold_utility::densityForDisplay(hVar, Form("h%s_pt%d_disp", systag.c_str(), ipt));
      hVarDisp->Scale(1./hVarDisp->Integral());
      hVarDisp->GetYaxis()->SetTitle("Shape-normalized counts / bin width");

      TH1D * hratio = (TH1D*)hVarDisp->Clone(Form("hratio_%s_pt%d", systag.c_str(), ipt));
      hratio->Divide(nominalDisp[ipt]);
      // Divide leaves 0 where the denominator is 0; NaN masking here blanked the page before.

      TH1D * hfrac = (TH1D*)hratio->Clone(Form("hfracdiff_%s_pt%d", systag.c_str(), ipt));
      for (int b = 1; b <= hfrac->GetNbinsX(); b++) {
        if (hratio->GetBinContent(b) != 0) hfrac->SetBinContent(b, hratio->GetBinContent(b) - 1);
      }
      hfrac->GetYaxis()->SetTitle("(Variation - Nominal) / Nominal");
      fracDiff[systag][ipt] = hfrac;

      fout->cd();
      hVarDisp->Write();
      hratio->Write();
      hfrac->Write();
    }
  }

  // Purity: symmetric quadrature sum of the five members.
  vector<TH1D*> purityUnc(ana::nPtBins);
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hPurity = (TH1D*)fracDiff[systematics[0]][ipt]->Clone(Form("hpurity_pt%d", ipt));
    hPurity->Reset("ICES");
    for (int b = 1; b <= hPurity->GetNbinsX(); b++) {
      double sumsq = 0;
      for (const string & systag : purityMembers) {
        double v = fracDiff[systag][ipt]->GetBinContent(b);
        sumsq += v*v;
      }
      hPurity->SetBinContent(b, sqrt(sumsq));
    }
    hPurity->GetYaxis()->SetTitle("Purity systematic uncertainty");
    purityUnc[ipt] = hPurity;
    fout->cd();
    hPurity->Write();
  }

  // Unfolding: up = niterHigh (+) prior, down = niterLow (+) prior. niter's direction is coherent
  // across xJ, so no per-bin sign test; the prior has no direction and enters both.
  vector<TH1D*> unfoldingUncUp(ana::nPtBins), unfoldingUncDown(ana::nPtBins);
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    // Null only if compute_prior_sensitivity.C has not been run for this radius.
    if (!fracDiff["priorSensitivity"][ipt]) {
      cout << "ERROR: fracDiff[\"priorSensitivity\"][" << ipt << "] is missing - run "
           << "macros/compute_prior_sensitivity.C(\"nominal\", " << ir << ") before compute_systematics.C(" << ir << ")." << endl;
      exit(1);
    }
    TH1D * hUp   = (TH1D*)fracDiff[systematics[0]][ipt]->Clone(Form("hunfolding_up_pt%d", ipt));
    TH1D * hDown = (TH1D*)fracDiff[systematics[0]][ipt]->Clone(Form("hunfolding_down_pt%d", ipt));
    hUp->Reset("ICES");
    hDown->Reset("ICES");
    for (int b = 1; b <= hUp->GetNbinsX(); b++) {
      double vHigh  = fracDiff["niterHigh"][ipt]->GetBinContent(b);
      double vLow   = fracDiff["niterLow"][ipt]->GetBinContent(b);
      double vPrior = fracDiff["priorSensitivity"][ipt]->GetBinContent(b);
      hUp->SetBinContent(b, sqrt(vHigh*vHigh + vPrior*vPrior));
      hDown->SetBinContent(b, sqrt(vLow*vLow + vPrior*vPrior));
    }
    hUp->GetYaxis()->SetTitle("Unfolding systematic uncertainty (up)");
    hDown->GetYaxis()->SetTitle("Unfolding systematic uncertainty (down)");
    unfoldingUncUp[ipt] = hUp;
    unfoldingUncDown[ipt] = hDown;
    fout->cd();
    hUp->Write();
    hDown->Write();
  }

  // Total: quadrature sum. Symmetric sources (including combined Purity) feed both totals;
  // asymmetricSystematics sources feed the total matching their sign in that bin.
  vector<TH1D*> totalUncUp(ana::nPtBins), totalUncDown(ana::nPtBins);
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hTotalUp   = (TH1D*)fracDiff[systematics[0]][ipt]->Clone(Form("hquadsum_up_pt%d", ipt));
    TH1D * hTotalDown = (TH1D*)fracDiff[systematics[0]][ipt]->Clone(Form("hquadsum_down_pt%d", ipt));
    hTotalUp->Reset("ICES");
    hTotalDown->Reset("ICES");
    for (int b = 1; b <= hTotalUp->GetNbinsX(); b++) {
      double sumsqUp = 0, sumsqDown = 0;
      for (const string & systag : systematics) {
        if (purityMembers.count(systag)) continue; // folded into purityUnc instead
        if (unfoldingMembers.count(systag)) continue; // folded into unfoldingUncUp/Down instead
        double v = fracDiff[systag][ipt]->GetBinContent(b);
        if (asymmetricSystematics.count(systag)) {
          if (v > 0) sumsqUp   += v*v;
          else       sumsqDown += v*v;
        } else {
          sumsqUp   += v*v;
          sumsqDown += v*v;
        }
      }
      double p = purityUnc[ipt]->GetBinContent(b);
      sumsqUp   += p*p;
      sumsqDown += p*p;
      double uUp   = unfoldingUncUp[ipt]->GetBinContent(b);
      double uDown = unfoldingUncDown[ipt]->GetBinContent(b);
      sumsqUp   += uUp*uUp;
      sumsqDown += uDown*uDown;
      hTotalUp->SetBinContent(b, sqrt(sumsqUp));
      hTotalDown->SetBinContent(b, sqrt(sumsqDown));
    }
    hTotalUp->GetYaxis()->SetTitle("Total systematic uncertainty (up)");
    hTotalDown->GetYaxis()->SetTitle("Total systematic uncertainty (down)");
    totalUncUp[ipt] = hTotalUp;
    totalUncDown[ipt] = hTotalDown;
    fout->cd();
    hTotalUp->Write();
    hTotalDown->Write();
  }

  fout->Close();
  cout << "Done. Wrote " << rootPath << endl;
}
