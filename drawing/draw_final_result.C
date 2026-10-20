#include "../src/ana.h"
#include "../src/drawer.h"
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Final-result pages from macros/compute_final_result.C's final_result.root: per pT bin, the unfolded
// Data (statistical error bars) with its systematic box vs Pythia8 truth, and the data/truth ratio;
// then all pT bins side by side on one page. Drawing only.
// Output: pdfs/final_result[_<rname>].pdf, pdfs/final_result_combined[_<rname>].pdf

void draw_final_result(int ir = 2) {
  gStyle->SetOptStat(0);
  const int nPtBinsUsed = ana::nPtBinsUsed;
  // R=0.4 keeps the unsuffixed filenames; other radii get _<rname>.
  string pdfPathStr  = (ir == 2) ? ana::path("pdfs/final_result.pdf")
                                 : Form("%s/pdfs/final_result_%s.pdf", ana::dir(), ana::rnames[ir]);
  string rootPathStr = (ir == 2) ? ana::path("hists/final_result.root")
                                 : Form("%s/hists/final_result_%s.root", ana::dir(), ana::rnames[ir]);
  string combinedPdfPathStr = (ir == 2) ? ana::path("pdfs/final_result_combined.pdf")
                                         : Form("%s/pdfs/final_result_combined_%s.pdf", ana::dir(), ana::rnames[ir]);
  const char * pdfPath = pdfPathStr.c_str();
  TFile * fin = TFile::Open(rootPathStr.c_str());
  if (!fin || fin->IsZombie()) {
    cout << "ERROR: cannot open " << rootPathStr << " - run macros/compute_final_result.C(" << ir << ") first." << endl;
    return;
  }

  drawer d;
  TCanvas * c = new TCanvas("c","",700,900);
  c->SaveAs(Form("%s[", pdfPath));

  const double ratioMin = 0.5, ratioMax = 1.5;
  // Arrow at the ratio frame edge for each off-scale central value.
  auto drawOffScaleArrows = [&](TH1D * hRatio, TH1D * hData, TH1D * hTruth) {
    const double arrowLen = 0.12 * (ratioMax - ratioMin);
    for (int b = 1; b <= hRatio->GetNbinsX(); b++) {
      double r = hRatio->GetBinContent(b);
      if (hData->GetBinContent(b) <= 0 || hTruth->GetBinContent(b) <= 0) continue;
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
  };
  // Kept for the combined three-panel page.
  vector<TH1D*> vData, vTruth, vRatio;
  vector<TGraphAsymmErrors*> vSystBox, vRatioSystBox;

  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hUnfoldDisp = (TH1D*)fin->Get(Form("hFinalDisp_pt%d", ipt));
    TH1D * hTruthDisp  = (TH1D*)fin->Get(Form("hFinalTruthDisp_pt%d", ipt));
    TH1D * hRatio      = (TH1D*)fin->Get(Form("hFinalRatio_pt%d", ipt));
    TGraphAsymmErrors * gSystBox      = (TGraphAsymmErrors*)fin->Get(Form("gFinalSystBox_pt%d", ipt));
    TGraphAsymmErrors * gRatioSystBox = (TGraphAsymmErrors*)fin->Get(Form("gFinalRatioSystBox_pt%d", ipt));

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
    l->AddEntry(hUnfoldDisp, ana::dataLabel(), "lp");
    l->AddEntry(gSystBox, "Systematic uncertainty", "f");
    l->AddEntry(hTruthDisp, ana::paperLabels() ? ana::mcLabel() : "Pythia8 #gamma+jet (truth)", "lp");
    l->Draw();
    d.drawAll({ana::dataLabel()},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir])}, .18, .85, 14, gPad->GetWh()*0.8);

    p2->cd();
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.3);
    p2->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    TH1D * hFrame2 = (TH1D*)hRatio->Clone(Form("hFinalFrame2_pt%d", ipt));
    hFrame2->Reset("ICES");
    hFrame2->SetLineColor(kWhite);
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
    hRatio->SetLineWidth(2); // as the data points
    hRatio->Draw("p e0 same");
    drawOffScaleArrows(hRatio, hUnfoldDisp, hTruthDisp);
    TLine * line = new TLine(ana::unfoldXjBins[0],1,ana::unfoldXjBins[ana::nUnfoldXjBins],1);
    line->SetLineStyle(9);
    line->Draw("same");
    c->SaveAs(pdfPath);
    vData.push_back(hUnfoldDisp);
    vTruth.push_back(hTruthDisp);
    vRatio.push_back(hRatio);
    vSystBox.push_back(gSystBox);
    vRatioSystBox.push_back(gRatioSystBox);
  }

  c->SaveAs(Form("%s]", pdfPath));

  // Combined page: the pT bins side by side, sPHENIX label in the first, legend in the second.
  // Pixel-size (font 43) text throughout, so it does not scale with the pad heights.
  TCanvas * cc = new TCanvas("ccombined","",1500,900);
  const double ySplit = .32;
  double ymaxAll = 0;
  for (int i = 0; i < nPtBinsUsed; i++)
    for (int b = 1; b <= vData[i]->GetNbinsX(); b++)
      ymaxAll = std::max({ymaxAll, vData[i]->GetBinContent(b) + vData[i]->GetBinError(b),
                          vTruth[i]->GetBinContent(b) + vTruth[i]->GetBinError(b)});
  auto pixelAxes = [](TH1D * h) {
    for (TAxis * a : {h->GetXaxis(), h->GetYaxis()}) {
      a->SetTitleFont(43);
      a->SetTitleSize(28);
      a->SetLabelFont(43);
      a->SetLabelSize(25);
      a->SetLabelOffset(a == h->GetXaxis() ? 0.03 : 0.012);
    }
  };
  for (int i = 0; i < nPtBinsUsed; i++) {
    int ipt = ana::firstUsedPtBin + i;
    double x1 = (double)i/nPtBinsUsed, x2 = (double)(i+1)/nPtBinsUsed;
    cc->cd();
    TPad * q1 = new TPad(Form("q1_%d",i),"",x1,ySplit,x2,1);
    TPad * q2 = new TPad(Form("q2_%d",i),"",x1,0,x2,ySplit);
    q1->Draw();
    q2->Draw();

    q1->cd();
    q1->SetBottomMargin(0.02);
    q1->SetTopMargin(0.03);
    q1->SetLeftMargin(.2);
    q1->SetRightMargin(.03);
    gPad->SetTicks(1,1);
    TH1D * hFrame1 = (TH1D*)vData[i]->Clone(Form("hCombFrame1_%d", i));
    hFrame1->Reset("ICES");
    hFrame1->SetLineColor(kWhite);
    pixelAxes(hFrame1);
    hFrame1->GetXaxis()->SetLabelSize(0);
    hFrame1->GetXaxis()->SetTitle("");
    // Range only clears the tallest error bar; the taller pad leaves room for the text above the data.
    hFrame1->GetYaxis()->SetRangeUser(0, ymaxAll*1.08);
    hFrame1->GetYaxis()->SetTitleOffset(1.75);
    hFrame1->Draw("p");
    vSystBox[i]->Draw("2 same");
    vData[i]->Draw("p e same");
    vTruth[i]->Draw("p e same");
    if (i == 0)
      d.drawAll({ana::dataLabel()}, {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ana::ptBins[ipt], ana::ptBins[ipt+1]),
                                     Form("Jet R=%.1f", ana::JetRs[ir]),
                                     Form("p_{T}^{jet} > %.0f GeV", ana::jet_calib_pt_cut[ir])},
                .25, .9, 22, gPad->GetWh()*(1-ySplit));  // spacing scales with the pad, not the canvas
    if (i == 1) {
      // Two columns keep it short enough to clear the peak.
      TLegend * l = new TLegend(.23,.79,.99,.895);
      l->SetLineWidth(0);
      l->SetFillStyle(0);
      l->SetNColumns(2);
      l->SetMargin(0.14);
      l->SetTextFont(43);
      l->SetTextSize(20);
      l->AddEntry(vData[i], ana::dataLabel(), "lp");
      l->AddEntry(vSystBox[i], "Syst. uncertainty", "f");
      l->AddEntry(vTruth[i], ana::paperLabels() ? ana::mcLabel() : "Pythia8 #gamma+jet (truth)", "lp");
      l->Draw();
    }
    // pT range: in the label block for the first panel, top left for the others
    if (i > 0) {
      TLatex * t = new TLatex(.25, .94, Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ana::ptBins[ipt], ana::ptBins[ipt+1]));
      t->SetNDC();
      t->SetTextFont(43);
      t->SetTextSize(21);
      t->SetTextAlign(13);
      t->Draw();
    }

    q2->cd();
    q2->SetTopMargin(0.02);
    q2->SetBottomMargin(0.4);
    q2->SetLeftMargin(.2);
    q2->SetRightMargin(.03);
    gPad->SetTicks(1,1);
    TH1D * hFrame2 = (TH1D*)vRatio[i]->Clone(Form("hCombFrame2_%d", i));
    hFrame2->Reset("ICES");
    hFrame2->SetLineColor(kWhite);
    pixelAxes(hFrame2);
    hFrame2->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
    hFrame2->GetYaxis()->SetTitle("Data / Pythia Truth");
    hFrame2->GetYaxis()->SetTitleOffset(1.75);
    hFrame2->GetYaxis()->SetNdivisions(505);
    hFrame2->GetXaxis()->SetTitle("x_{J#gamma}");
    hFrame2->GetXaxis()->SetTitleOffset(1.9);
    hFrame2->Draw("p");
    vRatioSystBox[i]->Draw("2 0 same");
    vRatio[i]->Draw("p e0 same");
    drawOffScaleArrows(vRatio[i], vData[i], vTruth[i]);
    TLine * line = new TLine(ana::unfoldXjBins[0],1,ana::unfoldXjBins[ana::nUnfoldXjBins],1);
    line->SetLineStyle(9);
    line->Draw("same");
  }
  cc->SaveAs(combinedPdfPathStr.c_str());
  cout << "Wrote " << pdfPathStr << " and " << combinedPdfPathStr << endl;
}
