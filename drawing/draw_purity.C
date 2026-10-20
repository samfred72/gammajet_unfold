#include "../src/drawer.h"
#include "../src/ana.h"
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Purity and MC leakage-fraction pages for one systag, one page per jet radius, from
// macros/puritymaker.C's ana::purityFilename(systag). Drawing only.
// Output: pdfs/purity_<systag>.pdf, pdfs/purity_leakage_<systag>.pdf
void draw_purity(string systag = "nominal") {
  gStyle->SetOptStat(0);
  drawer d;

  string purityFile = ana::purityFilename(systag);
  TFile * fin = TFile::Open(purityFile.c_str());
  if (!fin || fin->IsZombie()) {
    cout << "ERROR: cannot open " << purityFile << " - run macros/puritymaker.C(\"" << systag << "\") first." << endl;
    return;
  }

  TCanvas * cf = new TCanvas("cf","",700,700);
  TCanvas * cu = new TCanvas("cu","",700,700);
  string purityPdfPath = Form("%s/pdfs/purity_%s.pdf", ana::dir(), systag.c_str());
  string leakagePdfPath = Form("%s/pdfs/purity_leakage_%s.pdf", ana::dir(), systag.c_str());
  cu->SaveAs(Form("%s[", purityPdfPath.c_str()));
  cf->SaveAs(Form("%s[", leakagePdfPath.c_str()));

  const int colors[4] = {kBlack, kRed, kBlue, kOrange};
  const char * letters[4] = {"A","B","C","D"};
  for (int ir = 0; ir < ana::nJetR; ir++) {
    // MC leakage fractions f^B, f^C, f^D.
    cf->cd();
    cf->Clear();
    gPad->SetTicks();
    gPad->SetLeftMargin(.15);
    TH1F * framef = cf->DrawFrame(ana::ptBinsUsed[0], 0, ana::ptBinsUsed[ana::nPtBinsUsed], 1.2);
    framef->GetXaxis()->SetTitle("Leading cluster p_{T} [GeV]");
    framef->GetYaxis()->SetTitle("f^{X} = N^{X}_{sig}/N^{A}_{sig}");
    TLegend * lf = new TLegend(0.5, 0.65, 0.8, 0.8);
    for (int i = 1; i < 4; i++) {
      TH1D * fp = (TH1D*)fin->Get(Form("%s/leakage_%s", ana::rnames[ir], letters[i]));
      fp->SetLineColor(colors[i]);
      fp->SetLineWidth(2);
      fp->Draw("hist e same");
      lf->AddEntry(fp,Form("N_{%s}/N_{A}",letters[i]));
    }
    lf->SetLineWidth(0);
    lf->Draw();
    d.drawAll({ana::mcLabel()},{Form("systag: %s",systag.c_str()),Form("Jet R=%.1f",ana::JetRs[ir]),"truth-matched leakage fractions"},.5,.55,16,700);
    cf->SaveAs(leakagePdfPath.c_str());

    // Purity vs. photon pT, reported bins, with region C alongside.
    TGraphAsymmErrors * odUsed  = (TGraphAsymmErrors*)fin->Get(Form("%s/combined_used", ana::rnames[ir]));
    TGraphAsymmErrors * odCUsed = (TGraphAsymmErrors*)fin->Get(Form("%s/combined_C_used", ana::rnames[ir]));
    cu->cd();
    cu->Clear();
    gPad->SetTicks();
    gPad->SetLeftMargin(.15);
    // Headroom above P_A for the label (top left) and legend (top right).
    TH1F * frameu = cu->DrawFrame(ana::ptBinsUsed[0], 0, ana::ptBinsUsed[ana::nPtBinsUsed], 1.35);
    frameu->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
    frameu->GetYaxis()->SetTitle("Purity");
    odUsed->SetLineColor(kBlack);
    odUsed->SetMarkerColor(kBlack);
    odUsed->SetMarkerSize(1);
    odUsed->SetMarkerStyle(20);
    odUsed->SetLineWidth(2);
    odUsed->Draw("p same");
    odCUsed->SetLineColor(kAzure+2);
    odCUsed->SetMarkerColor(kAzure+2);
    odCUsed->SetMarkerSize(1);
    odCUsed->SetMarkerStyle(21);
    odCUsed->SetLineWidth(2);
    odCUsed->Draw("p same");
    TLegend * lu = new TLegend(.6,.74,.88,.88);
    lu->SetLineWidth(0);
    lu->AddEntry(odUsed,  "P_{A} (region A)");
    lu->AddEntry(odCUsed, "P_{C} (region C)");
    lu->Draw();
    d.drawAll({ana::dataLabel()},{Form("Jet R=%.1f",ana::JetRs[ir]),"leakage correction applied"},.19,.86,16,700);
    cu->SaveAs(purityPdfPath.c_str());
  }

  cu->SaveAs(Form("%s]", purityPdfPath.c_str()));
  cout << "Wrote " << purityPdfPath << endl;
  cf->SaveAs(Form("%s]", leakagePdfPath.c_str()));
  cout << "Wrote " << leakagePdfPath << endl;
}
