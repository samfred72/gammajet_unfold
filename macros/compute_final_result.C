#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Final result: nominal unfolded (1/N)dN/dxJ per pT bin vs #gamma+jet MC truth.
// Statistical errors are RooUnfoldBayes's analytic covariance (cross-checked by the toy macros),
// on the points. Systematics (hquadsum_up/down_pt<N> from compute_systematics.C) are asymmetric,
// kept as a separate box (TGraphAsymmErrors), also as a data/truth ratio box.
// Calculation only: drawing/draw_final_result.C draws the pages from final_result.root.

// Jet radius index, set by compute_final_result(int); default R=0.4.
int ir = 2;
const int nPtBinsUsed = ana::nPtBinsUsed; // reported bins start at ana::firstUsedPtBin
const int niterate = 2; // nominal iteration count
// R=0.4 keeps the unsuffixed filenames; other radii get _<rname> (must match compute_systematics).
string rootPathStr, systRootPathStr;
const char * rootPath;
const char * systRootPath;

void compute_final_result(int jetRadiusIndex = 2) {
  ir = jetRadiusIndex;
  rootPathStr     = (ir == 2) ? ana::path("hists/final_result.root")
                              : Form("%s/hists/final_result_%s.root", ana::dir(), ana::rnames[ir]);
  systRootPathStr = (ir == 2) ? ana::path("hists/systematics.root")
                              : Form("%s/hists/systematics_%s.root", ana::dir(), ana::rnames[ir]);
  rootPath     = rootPathStr.c_str();
  systRootPath = systRootPathStr.c_str();

  TFile * fsyst = TFile::Open(systRootPath);
  if (!fsyst || fsyst->IsZombie()) {
    cout << "ERROR: couldn't open " << systRootPath << " - run macros/compute_systematics.C first." << endl;
    return;
  }

  TFile * fout = TFile::Open(rootPath, "RECREATE");

  drawer d("pythia", "nominal");
  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i",ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i",ir), 1);
  TH1D * flatA = d.get(Form("hrecoxj%i_0",ir), 0);
  TH1D * flatC = d.get(Form("hrecoxj%i_2",ir), 0);
  TH1D * flatCorrected = unfold_utility::buildFullyCorrected(flatA, flatC, "data", "nominal", ir);
  TH1D * flatUnfolded  = unfold_utility::unfoldOnce(respRecoTemplate, respTruthTemplate, respMatrix2D, flatCorrected, niterate, "hUnfoldedFinal");

  cout << "pT bin, xJ bin: nominal, stat frac, syst frac, total frac" << endl;

  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hUnfold     = unfold_utility::unflattenXj(flatUnfolded, ipt, Form("hFinal_pt%d", ipt));
    TH1D * hUnfoldDisp = unfold_utility::densityForDisplay(hUnfold, Form("hFinalDisp_pt%d", ipt));
    hUnfoldDisp->Scale(1./hUnfoldDisp->Integral());
    hUnfoldDisp->GetYaxis()->SetTitle("(1/N) dN/dx_{J#gamma}");

    TH1D * hTruth     = unfold_utility::unflattenXj(respTruthTemplate, ipt, Form("hFinalTruth_pt%d", ipt));
    TH1D * hTruthDisp = unfold_utility::densityForDisplay(hTruth, Form("hFinalTruthDisp_pt%d", ipt));
    hTruthDisp->Scale(1./hTruthDisp->Integral());

    // Systematic box, separate from the statistical error bars; fractional, asymmetric.
    TH1D * hquadUp   = (TH1D*)fsyst->Get(Form("hquadsum_up_pt%d", ipt));
    TH1D * hquadDown = (TH1D*)fsyst->Get(Form("hquadsum_down_pt%d", ipt));
    int nb = hUnfoldDisp->GetNbinsX();
    TGraphAsymmErrors * gSystBox = new TGraphAsymmErrors(nb);
    gSystBox->SetName(Form("gFinalSystBox_pt%d", ipt));
    for (int b = 1; b <= nb; b++) {
      double content = hUnfoldDisp->GetBinContent(b);
      double x = hUnfoldDisp->GetBinCenter(b);
      double halfw = hUnfoldDisp->GetBinWidth(b) / 2.;
      double statFrac    = content > 0 ? hUnfoldDisp->GetBinError(b) / content : 0;
      double systFracUp   = (content > 0 && hquadUp)   ? hquadUp->GetBinContent(b)   : 0;
      double systFracDown = (content > 0 && hquadDown) ? hquadDown->GetBinContent(b) : 0;
      gSystBox->SetPoint(b-1, x, content);
      gSystBox->SetPointError(b-1, halfw, halfw, systFracDown*content, systFracUp*content);
      if (content > 0)
        cout << "  pt" << ipt << " bin" << b << ": " << content
             << ", stat=" << statFrac << ", syst up=" << systFracUp << ", syst down=" << systFracDown
             << ", total up=" << sqrt(statFrac*statFrac + systFracUp*systFracUp)
             << ", total down=" << sqrt(statFrac*statFrac + systFracDown*systFracDown) << endl;
    }

    // Ratio box: gSystBox scaled by the truth content. Only bins with both data and truth.
    TGraphAsymmErrors * gRatioSystBox = new TGraphAsymmErrors();
    gRatioSystBox->SetName(Form("gFinalRatioSystBox_pt%d", ipt));
    for (int b = 1; b <= nb; b++) {
      double truth = hTruthDisp->GetBinContent(b);
      double x, y;
      gSystBox->GetPoint(b-1, x, y);
      if (truth <= 0 || y <= 0) continue;
      int n = gRatioSystBox->GetN();
      gRatioSystBox->SetPoint(n, x, y/truth);
      gRatioSystBox->SetPointError(n, gSystBox->GetErrorXlow(b-1), gSystBox->GetErrorXhigh(b-1),
                                   gSystBox->GetErrorYlow(b-1)/truth, gSystBox->GetErrorYhigh(b-1)/truth);
    }
    TH1D * hRatio = (TH1D*)hUnfoldDisp->Clone(Form("hFinalRatio_pt%d", ipt));
    hRatio->Divide(hTruthDisp);
    for (int b = 1; b <= nb; b++) {
      if (hUnfoldDisp->GetBinContent(b) <= 0 || hTruthDisp->GetBinContent(b) <= 0) {
        hRatio->SetBinContent(b, 0);
        hRatio->SetBinError(b, 0);
      }
    }


    fout->cd();
    hUnfoldDisp->Write();
    gSystBox->Write();
    hTruthDisp->Write();
    hRatio->Write();
    gRatioSystBox->Write();
  }

  fout->Close();
  cout << "Done. Wrote " << rootPath << endl;
}
