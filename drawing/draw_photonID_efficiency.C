#include "../src/ana.h"
#include "../src/drawer.h"
#include "TLegend.h"
#include "TLine.h"
#include <algorithm>
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Does the fixed BDT/isolation cut give the same true-photon efficiency in every pT bin? Uses
// unfolder's hphoIDeff_bdt/iso (truth-matched photons, no ID cut). BDT passes above the cut,
// isolation below (as ana::findabcdBin).
// Returns the cumulative per-bin fraction (no interpolation).
TH1D * efficiencyCurve(TH1D * hscore, bool passAboveCut, const char * name) {
  int n = hscore->GetNbinsX();
  double total = hscore->Integral(0, n + 1);
  TH1D * heff = (TH1D*)hscore->Clone(name);
  heff->Reset("ICES");
  for (int i = 1; i <= n; i++) {
    double passing = passAboveCut ? hscore->Integral(i, n + 1) : hscore->Integral(0, i - 1);
    heff->SetBinContent(i, total > 0 ? passing / total : 0);
  }
  return heff;
}

void draw_photonID_efficiency(string systag = "nominal") {
  gStyle->SetOptStat(0);
  drawer d("pythia", systag);

  // Photon5/10/20 combined, as puritymaker.C.
  TH2D * hbdt = d.get2d("hphoIDeff_bdt", 1);
  TH2D * hiso = d.get2d("hphoIDeff_iso", 1);

  struct Source { TH2D * h; bool passAboveCut; double nominalCut; const char * label; const char * pdfName; };
  vector<Source> sources = {
    {hbdt, true,  ana::bdtGoodLow[0], "BDT score",           "photonID_efficiency_bdt"},
    {hiso, false, ana::isoBins[0],    "Isolation E_{T} [GeV]","photonID_efficiency_iso"},
  };

  int colors[ana::nPtBins] = {kGray+2, kBlue, kGreen+2, kOrange+1, kRed};

  for (auto& src : sources) {
    TCanvas * c = new TCanvas("c", "", 700, 700);
    gPad->SetLeftMargin(.15);
    gPad->SetBottomMargin(.13);
    gPad->SetTicks(1, 1);

    TH1F * frame = c->DrawFrame(src.h->GetXaxis()->GetBinLowEdge(1) < 0 ? -2 : 0,
        0, 1, 1.05);
    frame->GetXaxis()->SetTitle(src.label);
    frame->GetYaxis()->SetTitle("True-photon efficiency");
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);

    TLegend * leg = new TLegend(0.20, 0.55-0.2, 0.88-0.35, 0.88-0.2);
    leg->SetLineWidth(0);

    printf("\n=== %s: efficiency at nominal cut (%.3f) vs photon pT ===\n", src.label, src.nominalCut);
    double effAtNominal[ana::nPtBins];
    for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
      TH1D * hproj = src.h->ProjectionY(Form("%s_pt%i", src.h->GetName(), ipt), ipt + 1, ipt + 1);
      // Sanitize NaN/Inf before drawing.
      for (int b = 0; b <= hproj->GetNbinsX() + 1; b++) {
        double content = hproj->GetBinContent(b);
        if (std::isnan(content) || std::isinf(content)) hproj->SetBinContent(b, 0);
      }
      TH1D * heff = efficiencyCurve(hproj, src.passAboveCut, Form("heff_%s_pt%i", src.h->GetName(), ipt));
      heff->SetLineColor(colors[ipt]);
      heff->SetLineWidth(2);
      heff->Draw("hist same");
      leg->AddEntry(heff, Form("%.0f < p_{T}^{#gamma} < %.0f GeV%s",
          ana::ptBins[ipt], ana::ptBins[ipt + 1],
          (ipt == 0 || ipt == ana::nPtBins - 1) ? " (buffer bin)" : ""), "l");

      // Exact bin lookup.
      effAtNominal[ipt] = heff->GetBinContent(heff->GetXaxis()->FindBin(src.nominalCut));
      printf("  pT [%.0f,%.0f): efficiency = %.3f  (n_truth-matched = %.0f)\n",
          ana::ptBins[ipt], ana::ptBins[ipt + 1], effAtNominal[ipt], hproj->Integral(0, hproj->GetNbinsX() + 1));
    }
    double effMin = *std::min_element(effAtNominal, effAtNominal + ana::nPtBins);
    double effMax = *std::max_element(effAtNominal, effAtNominal + ana::nPtBins);
    printf("  --> spread across pT bins: %.3f - %.3f (Delta = %.3f)\n", effMin, effMax, effMax - effMin);

    TLine * lcut = new TLine(src.nominalCut, 0, src.nominalCut, 1.05);
    lcut->SetLineStyle(2);
    lcut->SetLineColor(kBlack);
    lcut->Draw();

    leg->Draw();
    d.drawAll({"Pythia8 #gamma+jet MC"}, {Form("systag: %s", systag.c_str()), "Truth-matched photons"}, .18, .3, 16, 700);

    const char * pdfPath = Form("%s/pdfs/%s_%s.pdf", ana::dir(), src.pdfName, systag.c_str());
    c->SaveAs(pdfPath);
    cout << "Wrote " << pdfPath << endl;
  }
}
