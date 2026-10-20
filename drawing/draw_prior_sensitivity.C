#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/unfold_utility.h"
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Prior-sensitivity pages from macros/compute_prior_sensitivity.C's output: per pT bin, the nominal
// unfolded result vs the sqrt(w), w, w^1.5 reweighted-prior unfoldings with their ratios; the w
// variant's fractional shift (the systematic, with its negation); its chi2/NDF vs niter. Drawing only.
// Output: pdfs/draw_prior_sensitivity_<systag>[_<rname>].pdf

const int niterate = 2; // nominal iteration count (dashed line on the last page)
const vector<string> variantNames = {"sqrtw", "w", "w1p5"};
const vector<string> variantLabels = {"#sqrt{w} prior", "w prior", "w^{1.5} prior"};
const vector<int> variantColors = {kAzure+2, kGreen+2, kRed};

void draw_prior_sensitivity(string systag = "nominal", int ir = 2) {
  gStyle->SetOptStat(0);
  const int nPtBinsUsed = ana::nPtBinsUsed;
  drawer d;
  // R=0.4 keeps the unsuffixed filenames; other radii get their own.
  string pdfPath  = (ir == 2) ? Form("%s/pdfs/draw_prior_sensitivity_%s.pdf", ana::dir(), systag.c_str())
                              : Form("%s/pdfs/draw_prior_sensitivity_%s_%s.pdf", ana::dir(), systag.c_str(), ana::rnames[ir]);
  string rootPath = (ir == 2) ? Form("%s/hists/prior_sensitivity_%s.root", ana::dir(), systag.c_str())
                              : Form("%s/hists/prior_sensitivity_%s_%s.root", ana::dir(), systag.c_str(), ana::rnames[ir]);
  TFile * fin = TFile::Open(rootPath.c_str());
  if (!fin || fin->IsZombie()) {
    cout << "ERROR: cannot open " << rootPath << " - run macros/compute_prior_sensitivity.C first." << endl;
    return;
  }

  TCanvas * c = new TCanvas("c","",700,900);
  c->SaveAs(Form("%s[", pdfPath.c_str()));

  // Pages 1..nPtBinsUsed: nominal vs the three variants (absolute densities) with ratios.
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hNominal = (TH1D*)fin->Get(Form("hNominal_pt%d", ipt));
    TH1D * hNominalDisp = unfold_utility::densityForDisplay(hNominal, Form("hNominalDisp_pt%d", ipt));
    vector<TH1D*> hVariantDisp(variantNames.size());
    double ymax = hNominalDisp->GetMaximum();
    for (unsigned iv = 0; iv < variantNames.size(); iv++) {
      TH1D * hVariant = (TH1D*)fin->Get(Form("h%s_pt%d", variantNames[iv].c_str(), ipt));
      hVariantDisp[iv] = unfold_utility::densityForDisplay(hVariant, Form("h%s_pt%d_disp", variantNames[iv].c_str(), ipt));
      ymax = std::max(ymax, hVariantDisp[iv]->GetMaximum());
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
    hNominalDisp->SetLineColor(kBlack);
    hNominalDisp->SetMarkerColor(kBlack);
    hNominalDisp->SetMarkerStyle(20);
    hNominalDisp->SetLineWidth(2);
    hNominalDisp->GetXaxis()->SetLabelSize(0);
    hNominalDisp->GetXaxis()->SetTitle("");
    hNominalDisp->GetYaxis()->SetRangeUser(0, ymax*1.4);
    hNominalDisp->Draw("p e");
    TLegend * l = new TLegend(.45,.58,.85,.80);
    l->SetLineWidth(0);
    l->SetTextSize(0.03);
    l->AddEntry(hNominalDisp, "Nominal (un-reweighted prior)", "lp");
    for (unsigned iv = 0; iv < variantNames.size(); iv++) {
      hVariantDisp[iv]->SetLineColor(variantColors[iv]);
      hVariantDisp[iv]->SetLineWidth(2);
      hVariantDisp[iv]->Draw("hist same");
      l->AddEntry(hVariantDisp[iv], variantLabels[iv].c_str(), "l");
    }
    l->Draw();
    d.drawAll({ana::dataLabel()},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir])}, .18, .85, 14, gPad->GetWh()*0.8);

    p2->cd();
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.3);
    p2->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    TH1D * hratioFrame = (TH1D*)hNominalDisp->Clone(Form("hratioFrame_pt%d", ipt));
    hratioFrame->Reset("ICES");
    hratioFrame->SetLineColor(kWhite);
    hratioFrame->SetMarkerColor(kWhite);
    hratioFrame->GetYaxis()->SetRangeUser(0.5,1.5);
    hratioFrame->GetYaxis()->SetTitle("Variant / Nominal");
    hratioFrame->GetYaxis()->SetTitleSize(0.09);
    hratioFrame->GetYaxis()->SetTitleOffset(0.7);
    hratioFrame->GetYaxis()->SetLabelSize(0.08);
    hratioFrame->GetXaxis()->SetTitle("x_{J#gamma}");
    hratioFrame->GetXaxis()->SetTitleSize(0.09);
    hratioFrame->GetXaxis()->SetLabelSize(0.08);
    hratioFrame->Draw("p");
    for (unsigned iv = 0; iv < variantNames.size(); iv++) {
      TH1D * hratio = (TH1D*)fin->Get(Form("hratio_%s_pt%d", variantNames[iv].c_str(), ipt));
      hratio->SetLineColor(variantColors[iv]);
      hratio->SetMarkerColor(variantColors[iv]);
      hratio->SetMarkerStyle(20);
      hratio->SetLineWidth(2); // as the variant curves
      hratio->Draw("p e same");
    }
    TLine * line = new TLine(ana::unfoldXjBins[0],1,ana::unfoldXjBins[ana::nUnfoldXjBins],1);
    line->SetLineStyle(9);
    line->Draw("same");
    c->SaveAs(pdfPath.c_str());
  }

  // Summary: the w variant's fractional shift (the systematic), drawn with its negation as a
  // symmetric source.
  c->Clear();
  c->cd();
  vector<TPad*> pads(nPtBinsUsed);
  for (int idisplay = 0; idisplay < nPtBinsUsed; idisplay++) {
    double y1 = 1.0 - (idisplay+1)*(1.0/nPtBinsUsed);
    double y2 = 1.0 - idisplay*(1.0/nPtBinsUsed);
    pads[idisplay] = new TPad(Form("psum_%d",idisplay),"",0,y1,1,y2);
    pads[idisplay]->Draw();
  }
  for (int idisplay = 0; idisplay < nPtBinsUsed; idisplay++) {
    int ipt = idisplay + ana::firstUsedPtBin;
    pads[idisplay]->cd();
    pads[idisplay]->SetLeftMargin(.15);
    pads[idisplay]->SetBottomMargin(idisplay == nPtBinsUsed-1 ? 0.2 : 0.02);
    pads[idisplay]->SetTopMargin(0.05);
    gPad->SetTicks(1,1);

    TH1D * hFrac    = (TH1D*)fin->Get(Form("hPriorSensFracDiff_pt%d", ipt));
    TH1D * hFracNeg = (TH1D*)hFrac->Clone(Form("hFracDiffNeg_pt%d", ipt));
    hFracNeg->Scale(-1);
    double ymax = 0.05;
    for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) ymax = std::max(ymax, std::abs(hFrac->GetBinContent(ixj+1)));
    hFrac->GetYaxis()->SetRangeUser(-ymax*1.3, ymax*1.3);
    hFrac->GetYaxis()->SetTitle("(w - Nom.)/Nom.");
    hFrac->GetYaxis()->SetTitleSize(0.08);
    hFrac->GetYaxis()->SetTitleOffset(0.8);
    hFrac->GetYaxis()->SetLabelSize(0.07);
    hFrac->GetXaxis()->SetTitle(idisplay == nPtBinsUsed-1 ? "x_{J#gamma}" : "");
    hFrac->GetXaxis()->SetLabelSize(idisplay == nPtBinsUsed-1 ? 0.07 : 0);
    hFrac->GetXaxis()->SetTitleSize(0.08);
    hFrac->SetLineColor(kBlack);
    hFrac->SetLineWidth(2);
    hFrac->Draw("hist");
    hFracNeg->SetLineColor(kGray+1);
    hFracNeg->SetLineStyle(2);
    hFracNeg->SetLineWidth(2);
    hFracNeg->Draw("hist same");
    TLine * zero = new TLine(ana::unfoldXjBins[0],0,ana::unfoldXjBins[ana::nUnfoldXjBins],0);
    zero->SetLineStyle(9);
    zero->Draw("same");
    TLatex * t = new TLatex(.18,.85,Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]));
    t->SetNDC();
    t->SetTextFont(43);
    t->SetTextSize(16);
    t->Draw();
  }
  c->SaveAs(pdfPath.c_str());

  // Mean chi2/NDF of the w variant vs nominal at each niter (should fall with niter).
  TGraph * gChi2 = (TGraph*)fin->Get("gPriorSensChi2");
  c->Clear();
  c->cd();
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(.15);
  gChi2->SetMarkerStyle(20);
  gChi2->SetMarkerColor(kBlack);
  gChi2->SetLineColor(kBlack);
  gChi2->SetLineWidth(2);
  gChi2->GetXaxis()->SetTitle("Bayesian unfolding iterations");
  gChi2->GetYaxis()->SetTitle("Prior-sensitivity #chi^{2}/NDF (w vs nominal)");
  double ymaxg = 0;
  for (int k = 0; k < gChi2->GetN(); k++) { double x,y; gChi2->GetPoint(k,x,y); ymaxg = std::max(ymaxg,y); }
  gChi2->SetMinimum(0);
  gChi2->SetMaximum(ymaxg*1.3);
  gChi2->Draw("APL");
  TLine * lnom = new TLine(niterate, 0, niterate, ymaxg*1.3);
  lnom->SetLineStyle(9);
  lnom->SetLineColor(kRed);
  lnom->Draw("same");
  d.drawAll({ana::dataLabel()},{Form("Jet R=%.1f",ana::JetRs[ir]),
      Form("Nominal: %d iterations (dashed line)",niterate)}, .5, .85, 16, 700);
  c->SaveAs(pdfPath.c_str());

  c->SaveAs(Form("%s]", pdfPath.c_str()));
  cout << "Wrote " << pdfPath << endl;
}
