#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Final result: nominal unfolded (1/N)dN/dxJ per pT bin vs #gamma+jet MC truth.
// Statistical errors are RooUnfoldBayes's analytic covariance (cross-checked by the toy macros),
// on the points. Systematics (hquadsum_up/down_pt<N> from draw_systematics.C) are asymmetric,
// drawn as a separate box (TGraphAsymmErrors). Ratio panel below.

// Jet radius index, set by draw_final_result(int); default R=0.4.
int ir = 2;
const int nPtBinsUsed = ana::nPtBinsUsed; // reported bins start at ana::firstUsedPtBin
const int niterate = 2; // nominal iteration count
// R=0.4 keeps the unsuffixed filenames; other radii get _<rname> (must match draw_systematics).
string pdfPathStr, rootPathStr, systRootPathStr;
const char * pdfPath;
const char * rootPath;
const char * systRootPath;

void draw_final_result(int jetRadiusIndex = 2) {
  gStyle->SetOptStat(0);

  ir = jetRadiusIndex;
  pdfPathStr      = (ir == 2) ? ana::path("pdfs/final_result.pdf")
                              : Form("%s/pdfs/final_result_%s.pdf", ana::dir(), ana::rnames[ir]);
  rootPathStr     = (ir == 2) ? ana::path("hists/final_result.root")
                              : Form("%s/hists/final_result_%s.root", ana::dir(), ana::rnames[ir]);
  systRootPathStr = (ir == 2) ? ana::path("hists/systematics.root")
                              : Form("%s/hists/systematics_%s.root", ana::dir(), ana::rnames[ir]);
  pdfPath      = pdfPathStr.c_str();
  rootPath     = rootPathStr.c_str();
  systRootPath = systRootPathStr.c_str();

  TFile * fsyst = TFile::Open(systRootPath);
  if (!fsyst || fsyst->IsZombie()) {
    cout << "ERROR: couldn't open " << systRootPath << " - run drawing/draw_systematics.C first." << endl;
    return;
  }

  TFile * fout = TFile::Open(rootPath, "RECREATE");
  TCanvas * c = new TCanvas("c","",700,900);
  c->SaveAs(Form("%s[", pdfPath));

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

    c->Clear();
    c->cd();
    TPad * p1 = new TPad(Form("p1_%d",ipt),"",0,.35,1,1);
    TPad * p2 = new TPad(Form("p2_%d",ipt),"",0,0,1,.35);
    p1->Draw();
    p2->Draw();

    p1->cd();
    p1->SetBottomMargin(0.02);
    p1->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    // Frame-only histogram for the axes.
    TH1D * hFrame1 = (TH1D*)hUnfoldDisp->Clone(Form("hFinalFrame1_pt%d", ipt));
    hFrame1->Reset("ICES");
    hFrame1->SetLineColor(kWhite);
    hFrame1->GetXaxis()->SetLabelSize(0);
    hFrame1->GetXaxis()->SetTitle("");
    hFrame1->GetYaxis()->SetRangeUser(0, std::max(hUnfoldDisp->GetMaximum(), hTruthDisp->GetMaximum())*1.5);
    hFrame1->Draw("p");
    gSystBox->SetFillColorAlpha(kAzure+1, 0.35);
    gSystBox->SetLineColor(kWhite);
    gSystBox->Draw("2 same");
    hUnfoldDisp->SetLineColor(kBlack);
    hUnfoldDisp->SetMarkerColor(kBlack);
    hUnfoldDisp->SetMarkerStyle(20);
    hUnfoldDisp->SetLineWidth(2);
    hUnfoldDisp->Draw("p e same");
    hTruthDisp->SetLineColor(kRed);
    hTruthDisp->SetMarkerColor(kRed);
    hTruthDisp->SetMarkerStyle(24);
    hTruthDisp->SetLineWidth(2);
    hTruthDisp->Draw("p e same");
    TLegend * l = new TLegend(.55,.62,.85,.85);
    l->SetLineWidth(0);
    l->SetTextSize(0.032);
    l->AddEntry(hUnfoldDisp, "p+p Run24 Data", "lp");
    l->AddEntry(gSystBox, "Systematic uncertainty", "f");
    l->AddEntry(hTruthDisp, "Pythia8 #gamma+jet (truth)", "lp");
    l->Draw();
    d.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir])}, .18, .85, 14, gPad->GetWh()*0.8);

    p2->cd();
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.3);
    p2->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
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

    TH1D * hFrame2 = (TH1D*)hRatio->Clone(Form("hFinalFrame2_pt%d", ipt));
    hFrame2->Reset("ICES");
    hFrame2->SetLineColor(kWhite);
    const double ratioMin = 0.5, ratioMax = 1.5;
    hFrame2->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
    hFrame2->GetYaxis()->SetTitle("Data / Pythia Truth");
    hFrame2->GetYaxis()->SetTitleSize(0.09);
    hFrame2->GetYaxis()->SetTitleOffset(0.7);
    hFrame2->GetYaxis()->SetLabelSize(0.08);
    hFrame2->GetXaxis()->SetTitle("x_{J#gamma}");
    hFrame2->GetXaxis()->SetTitleSize(0.09);
    hFrame2->GetXaxis()->SetLabelSize(0.08);
    hFrame2->Draw("p");
    gRatioSystBox->SetFillColorAlpha(kAzure+1, 0.35);
    gRatioSystBox->SetLineColor(kWhite);
    // "2 0" (graph) and "e0" (histogram) keep off-frame points' errors and boxes, clipped; an arrow
    // at the frame edge marks each off-scale central value.
    gRatioSystBox->Draw("2 0 same");
    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetMarkerStyle(20);
    hRatio->Draw("p e0 same");
    const double arrowLen = 0.12 * (ratioMax - ratioMin);
    for (int b = 1; b <= nb; b++) {
      double r = hRatio->GetBinContent(b);
      if (hUnfoldDisp->GetBinContent(b) <= 0 || hTruthDisp->GetBinContent(b) <= 0) continue;
      if (r >= ratioMin && r <= ratioMax) continue;
      double xc = hRatio->GetBinCenter(b);
      double yEdge = (r > ratioMax) ? ratioMax : ratioMin;
      double yTail = (r > ratioMax) ? ratioMax - arrowLen : ratioMin + arrowLen;
      TArrow * arr = new TArrow(xc, yTail, xc, yEdge, 0.015, "|>");
      arr->SetLineColor(kBlack);
      arr->SetFillColor(kBlack);
      arr->SetLineWidth(2);
      arr->Draw();
    }
    TLine * line = new TLine(ana::unfoldXjBins[0],1,ana::unfoldXjBins[ana::nUnfoldXjBins],1);
    line->SetLineStyle(9);
    line->Draw("same");
    c->SaveAs(pdfPath);

    fout->cd();
    hUnfoldDisp->Write();
    gSystBox->Write();
    hTruthDisp->Write();
  }

  c->SaveAs(Form("%s]", pdfPath));
  fout->Close();
  cout << "Done. Wrote " << pdfPath << " and " << rootPath << endl;
}
