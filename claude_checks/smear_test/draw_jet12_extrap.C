#include "../../src/ana.h"
#include "../../src/drawer.h"
// Explicit load - see drawing/draw_final_result.C. Run interpreted, never with ACLiC "+".
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Drawing half of fill_jet12_extrap.C (which runs on SDCC; see its header for the four jet-pT
// definitions, the selection and the missing-jet estimate).
//
//   pdfs/jet12_extrap_widths.pdf          left: sigma/pT for template, lin15, lin5
//                                         right: lin5 expected missing jets (storage floor) vs the
//                                         lin5 jet spectra
//   pdfs/jet12_extrap_jetpt_R<RR>_<f>.pdf leading / subleading / subsubleading jet pT, four definitions,
//                                         ratios to the template smear

void draw_jet12_extrap(int radius = 4, bool useTruthRef = true) {
  gStyle->SetOptStat(0);
  drawer d("pythia", "nominal");
  const float R = radius / 10.0f;
  const char *flav = useTruthRef ? "truth" : "reco";

  TFile fjer("jerband_smearing_templates.root");
  TH1D *hNom = (TH1D*)fjer.Get("h_jer_smear_r04_pileup_EMfracJES_nominal");
  TH1D *hUp  = (TH1D*)fjer.Get("h_jer_smear_r04_pileup_EMfracJES_sysup");
  hNom->SetDirectory(0); hUp->SetDirectory(0);
  fjer.Close();
  const double p15 = 15, w15 = hNom->Interpolate(p15);
  const double s15 = (hNom->Interpolate(p15 + 0.01) - hNom->Interpolate(p15 - 0.01)) / 0.02;
  const double c1 = hNom->GetBinCenter(1), w5 = hNom->GetBinContent(1);
  const double s5 = (hNom->GetBinContent(2) - hNom->GetBinContent(1)) / (hNom->GetBinCenter(2) - c1);
  auto width = [&](int sc, double pt) -> double {
    if (sc == 1) return hNom->Interpolate(pt);
    if (sc == 2) return pt >= p15 ? hNom->Interpolate(pt) : w15 + s15*(pt - p15);
    if (sc == 3) return pt >= c1 ? hNom->Interpolate(pt) : w5 + s5*(pt - c1);
    return hUp->Interpolate(pt);
  };

  const int nScen = 4, nRank = 3;
  const char *scenNames[nScen] = {"calib", "tmpl", "lin15", "lin5"};
  const char *rankNames[nRank] = {"leading", "subleading", "subsubleading"};
  const char *scenLabels[nScen] = {"Unsmeared (p_{T}^{calib})", "Template (treemaker, nominal)",
                                   "Tangent at 15 GeV, below 15", "Slope at 5 GeV, below 5.19"};
  const int col[nScen] = {kGray+1, kBlack, kAzure+1, kRed+1};

  TFile *fin = TFile::Open(Form("jet12_extrap_R%02d_%s.root", radius, flav));
  TH1D *h[nScen][nRank], *hMiss5;
  for (int sc = 0; sc < nScen; sc++)
    for (int r = 0; r < nRank; r++) {
      h[sc][r] = (TH1D*)fin->Get(Form("h_%s_%s", scenNames[sc], rankNames[r]));
      h[sc][r]->SetDirectory(0);
    }
  hMiss5 = (TH1D*)fin->Get("h_missing_lin5"); hMiss5->SetDirectory(0);
  const double truthLeadMin = ((TParameter<double>*)fin->Get("truthLeadMin"))->GetVal();
  fin->Close();
  for (int sc = 0; sc < nScen; sc++) for (int r = 0; r < nRank; r++)
    for (int b = 0; b <= h[sc][r]->GetNbinsX() + 1; b++)
      if (!std::isfinite(h[sc][r]->GetBinContent(b))) h[sc][r]->SetBinContent(b, 0);
  for (int b = 0; b <= hMiss5->GetNbinsX() + 1; b++)
    if (!std::isfinite(hMiss5->GetBinContent(b))) hMiss5->SetBinContent(b, 0);

  // ---- widths, and the lin5 missing-jet estimate ----
  {
    TCanvas *c = new TCanvas("cw", "", 1400, 650);
    c->Divide(2, 1);
    TPad *p = (TPad*)c->cd(1);
    p->SetLeftMargin(.14); p->SetRightMargin(.04); p->SetBottomMargin(.12); p->SetTopMargin(.05);
    TH1F *fr = p->DrawFrame(0, 0, 20, 1.4);
    fr->GetXaxis()->SetTitle("p_{T}^{ref} [GeV]"); fr->GetYaxis()->SetTitle("#sigma / p_{T}");
    fr->GetXaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleOffset(1.4);
    TLegend *leg = new TLegend(.4, .5, .93, .7);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);
    for (int sc : {1, 2, 3, 0}) {
      TGraph *g = new TGraph();
      for (int k = 0; k <= 400; k++) { double pt = 20.0*k/400; g->SetPoint(k, pt, width(sc, pt)); }
      g->SetLineColor(sc == 0 ? kBlack : col[sc]); g->SetLineWidth(2); g->SetLineStyle(sc == 0 ? 3 : 1);
      g->Draw("L same");
      leg->AddEntry(g, sc == 0 ? "Template sys. up (storage floor)" : scenLabels[sc], "l");
    }
    leg->Draw();
    d.drawAll({"Pythia8 Jet12 MC"}, {"JER template r04_pileup_EMfracJES", "Nominal widths"}, .42, .88, 22, 650);

    p = (TPad*)c->cd(2);
    p->SetLeftMargin(.14); p->SetRightMargin(.04); p->SetBottomMargin(.12); p->SetTopMargin(.05); p->SetLogy();
    TH1D *fr2 = (TH1D*)h[3][2]->Clone("fr2"); fr2->Reset();
    fr2->GetXaxis()->SetRangeUser(0, 15); fr2->SetMinimum(10); fr2->SetMaximum(h[3][1]->GetMaximum() * 2000);
    fr2->GetYaxis()->SetTitle("Jets"); fr2->GetXaxis()->SetTitleSize(0.045); fr2->GetYaxis()->SetTitleSize(0.045);
    fr2->Draw("axis");
    const int rcol[nRank] = {kRed+3, kRed+1, kRed-7};
    TLegend *leg2 = new TLegend(.36, .7, .93, .93);
    leg2->SetBorderSize(0); leg2->SetFillStyle(0); leg2->SetTextSize(0.03);
    for (int r = 0; r < nRank; r++) {
      h[3][r]->SetLineColor(rcol[r]); h[3][r]->SetLineWidth(2); h[3][r]->SetMarkerSize(0);
      h[3][r]->Draw("hist same");
      leg2->AddEntry(h[3][r], Form("Slope at 5 GeV, %s", rankNames[r]), "l");
    }
    hMiss5->SetLineColor(kBlack); hMiss5->SetFillColor(kGray); hMiss5->SetFillStyle(1001); hMiss5->SetLineWidth(1);
    hMiss5->Draw("hist same");
    for (int r = 0; r < nRank; r++) h[3][r]->Draw("hist same");
    p->RedrawAxis();
    leg2->AddEntry(hMiss5, "Est. missing (not stored), all ranks", "f");
    leg2->Draw();
    // a missing jet (pt_calib c >= 1 GeV, z below the sys-up threshold) lands at most at c + (4 - c) w5(c)/w_up(c)
    TLatex tm; tm.SetNDC(); tm.SetTextSize(0.032);
    tm.DrawLatex(.4, .64, "Missing jets lie at 4 - 6.4 GeV only;");
    tm.DrawLatex(.4, .60, "spectra above that are complete");
    c->SaveAs("pdfs/jet12_extrap_widths.pdf");
  }

  // ---- jet pT spectra and ratios ----
  TCanvas *c = new TCanvas("c", "", 1800, 800);
  for (int r = 0; r < nRank; r++) {
    const double x0 = r / 3.0, x1 = (r + 1) / 3.0;
    TPad *pt = new TPad(Form("pt%d", r), "", x0, 0.32, x1, 1);
    TPad *pb = new TPad(Form("pb%d", r), "", x0, 0, x1, 0.32);
    for (TPad *p : {pt, pb}) { p->SetLeftMargin(.15); p->SetRightMargin(.04); }
    pt->SetBottomMargin(0.02); pt->SetTopMargin(.05); pb->SetTopMargin(0.02); pb->SetBottomMargin(.3);
    c->cd(); pt->Draw(); pb->Draw();

    const double xmax = (r == 0) ? 50 : (r == 1) ? 40 : 25;
    pt->cd(); pt->SetLogy();
    double ymax = 0;
    for (int sc = 0; sc < nScen; sc++) ymax = std::max(ymax, h[sc][r]->GetMaximum());
    TH1D *fr = (TH1D*)h[1][r]->Clone(Form("fr%d", r)); fr->Reset();
    fr->GetXaxis()->SetRangeUser(0, xmax); fr->SetMinimum(0.5); fr->SetMaximum(ymax * 300);
    fr->GetXaxis()->SetLabelSize(0); fr->GetYaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetLabelSize(0.045);
    fr->GetYaxis()->SetTitleOffset(1.3);
    fr->Draw("axis");
    for (int sc = 0; sc < nScen; sc++) {
      h[sc][r]->SetLineColor(col[sc]); h[sc][r]->SetLineWidth(2); h[sc][r]->SetMarkerSize(0);
      h[sc][r]->Draw("hist same");
    }
    TLine *lf = new TLine(4, 0.5, 4, ymax * 300); lf->SetLineStyle(3); lf->SetLineColor(kGray+2); lf->Draw();
    // selection thresholds: leading jet 20 GeV, subleading 7 GeV (multijet/analysis.cc)
    const double cutLine = (r == 0) ? 20 : (r == 1) ? 7 : -1;
    if (cutLine > 0) { TLine *lc = new TLine(cutLine, 0.5, cutLine, ymax * 300); lc->SetLineStyle(3); lc->SetLineWidth(2); lc->Draw(); }
    if (r == 0)
      d.drawAll({"Pythia8 Jet12 MC"}, {Form("Jet R=%.1f, |#eta_{det}| < %.1f", R, 1.1 - R),
                 Form("p_{T}^{truth,lead} > %.0f GeV", truthLeadMin), Form("p_{T}^{ref}: %s", flav)}, .53, .9, 18, 544);
    if (r == 1) {
      TLegend *leg = new TLegend(.42, .66, .95, .93);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);
      for (int sc = 0; sc < nScen; sc++) leg->AddEntry(h[sc][r], scenLabels[sc], "l");
      leg->Draw();
    }
    if (r == 2) {
      TLatex tw; tw.SetNDC(); tw.SetTextSize(0.035); tw.SetTextColor(kRed+1);
      tw.DrawLatex(.4, .88, "Slope at 5 GeV: jets missing at 4-6.4 GeV");
      tw.DrawLatex(.4, .84, "(storage floor), see widths pdf");
    }
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.05);
    tl.DrawLatex(.2, .08, rankNames[r]);

    pb->cd();
    TH1D *rf = (TH1D*)fr->Clone(Form("rf%d", r));
    rf->SetMinimum(0.5); rf->SetMaximum(1.5);
    rf->GetYaxis()->SetTitle("Ratio to template"); rf->GetYaxis()->SetNdivisions(505);
    rf->GetXaxis()->SetLabelSize(0.1); rf->GetXaxis()->SetTitleSize(0.11); rf->GetXaxis()->SetTitleOffset(1.1);
    rf->GetYaxis()->SetLabelSize(0.09); rf->GetYaxis()->SetTitleSize(0.09); rf->GetYaxis()->SetTitleOffset(0.75);
    rf->Draw("axis");
    TLine *l1 = new TLine(0, 1, xmax, 1); l1->SetLineStyle(2); l1->Draw();
    if (cutLine > 0) { TLine *lc = new TLine(cutLine, 0.5, cutLine, 1.5); lc->SetLineStyle(3); lc->SetLineWidth(2); lc->Draw(); }
    for (int sc : {0, 2, 3}) {
      TH1D *ra = (TH1D*)h[sc][r]->Clone(Form("ra%d_%d", r, sc));
      ra->Divide(h[1][r]);
      for (int b = 0; b <= ra->GetNbinsX() + 1; b++)
        if (!std::isfinite(ra->GetBinContent(b)) || h[1][r]->GetBinContent(b) <= 0) { ra->SetBinContent(b, 0); ra->SetBinError(b, 0); }
      ra->SetLineColor(col[sc]); ra->SetLineWidth(2);
      ra->Draw("hist same"); // same events and z in all: bin errors would be meaningless
    }
  }
  c->SaveAs(Form("pdfs/jet12_extrap_jetpt_R%02d_%s.pdf", radius, flav));
}
