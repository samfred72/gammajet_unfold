#include "../../src/ana.h"
#include "../../src/drawer.h"
// Explicit load - see drawing/draw_final_result.C. Run interpreted, never with ACLiC "+".
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Drawing half (Alpine login node is fine); filling is fill_jet12_smear.C on SDCC.
//
// Jet12 (Pythia) leading / subleading / subsubleading jet pT with the treemaker's JER smearing
// vs. a version whose width is linearly extrapolated below 10 GeV.
//
// Treemaker smearing (multiJet DijetTreeMaker.cc:1739-1760, smear_pt() at :2373; same scheme as
// gammajet CaloAna.cc): smeared = pt_calib + z * pt_ref * w(pt_ref), one standard-normal z per jet
// shared by nominal/up/down, w = h_jer_smear_r04_pileup_EMfracJES_*->Interpolate(pt_ref). The
// template starts at 5 GeV and Interpolate holds w flat below 5.19 GeV (see draw_smear_function.C).
// The tree stores pt_calib, z (jet_smear_z_R) and the matched truth pT (jet_truth_pt_R, -1 if
// unmatched), so each jet is re-smeared here with the same z and only w changed:
//
//   template  w(pt) = Interpolate(pt)                                (as in the treemaker)
//   linear    w(pt) = Interpolate(pt) for pt >= 10 GeV,
//             w(10) + w'(10) * (pt - 10) below, w'(10) the template's slope at 10 GeV
//
// separately for nominal, sys up and sys down. useTruthRef picks the "_truth" flavour (pt_ref = matched
// truth pT, else pt_calib; what CaloAna uses for MC) or "_reco" (pt_ref = pt_calib). The recomputed
// template smear is checked against the stored jet_pt_smear_{truth,reco}_R branch.
//
// Selection: |z_vtx| < 60 cm and leading truth jet pT above the Jet12 full-efficiency threshold, as
// multijet/analysis.cc; jets with |eta_det| < 1.1 - R, ranked by pT separately in each scenario.
// Storage floor: the treemaker only stores a jet if raw pT, pt_calib or one of its six template smears
// reaches 4 GeV (DijetTreeMaker.cc:1767). With the shared z, a linear smear can only pass where a
// template smear did not if its width is larger than the template sys-up width: the linear nominal and
// down never are, and linear up only below pT_ref ~ 0.9 GeV (not reachable to 4 GeV). So no stored-jet
// loss affects these spectra.

void draw_jet12_smear(int radius = 4, bool useTruthRef = true) {
  gStyle->SetOptStat(0);
  drawer d("pythia", "nominal");
  const float R = radius / 10.0f;
  const double linPivot = 10;
  const char *flav = useTruthRef ? "truth" : "reco";

  TFile fjer("jerband_smearing_templates.root");
  const char *jerNames[3] = {"nominal", "sysup", "sysdown"};
  TH1D *hW[3];
  double w10[3], slope[3];
  for (int s = 0; s < 3; s++) {
    hW[s] = (TH1D*)fjer.Get(Form("h_jer_smear_r04_pileup_EMfracJES_%s", jerNames[s]));
    hW[s]->SetDirectory(0);
    w10[s] = hW[s]->Interpolate(linPivot);
    slope[s] = (hW[s]->Interpolate(linPivot + 0.01) - hW[s]->Interpolate(linPivot - 0.01)) / 0.02;
  }
  fjer.Close();
  auto wTmpl = [&](int s, double pt) { return hW[s]->Interpolate(pt); };
  auto wLin  = [&](int s, double pt) { return pt >= linPivot ? hW[s]->Interpolate(pt) : w10[s] + slope[s]*(pt - linPivot); };

  const int nScen = 7, nRank = 3;
  const char *scenNames[nScen] = {"calib", "tmpl_nom", "tmpl_up", "tmpl_down", "lin_nom", "lin_up", "lin_down"};
  const char *rankNames[nRank] = {"leading", "subleading", "subsubleading"};
  TFile *fin = TFile::Open(Form("jet12_smear_R%02d_%s.root", radius, flav));
  TH1D *h[nScen][nRank];
  for (int sc = 0; sc < nScen; sc++)
    for (int r = 0; r < nRank; r++) {
      h[sc][r] = (TH1D*)fin->Get(Form("h_%s_%s", scenNames[sc], rankNames[r]));
      h[sc][r]->SetDirectory(0);
      for (int b = 0; b <= h[sc][r]->GetNbinsX() + 1; b++)
        if (!std::isfinite(h[sc][r]->GetBinContent(b))) h[sc][r]->SetBinContent(b, 0);
    }
  const double truthLeadMin = ((TParameter<double>*)fin->Get("truthLeadMin"))->GetVal();
  fin->Close();

  // ---- width functions ----
  {
    TCanvas *c = new TCanvas("cw", "", 750, 650);
    c->SetLeftMargin(.14); c->SetRightMargin(.04); c->SetBottomMargin(.12); c->SetTopMargin(.05);
    TH1F *fr = c->DrawFrame(0, 0, 20, 0.9);
    fr->GetXaxis()->SetTitle("p_{T}^{ref} [GeV]"); fr->GetYaxis()->SetTitle("#sigma / p_{T}");
    fr->GetXaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleOffset(1.4);
    const int col[3] = {kBlack, kRed+1, kBlue+1};
    TLegend *leg = new TLegend(.5, .45, .93, .7);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);
    for (int s = 0; s < 3; s++) {
      TGraph *gt = new TGraph(), *gl = new TGraph();
      for (int k = 0; k <= 400; k++) { double pt = 20.0*k/400; gt->SetPoint(k, pt, wTmpl(s, pt)); gl->SetPoint(k, pt, wLin(s, pt)); }
      gt->SetLineColor(col[s]); gt->SetLineWidth(2);
      gl->SetLineColor(col[s]); gl->SetLineWidth(2); gl->SetLineStyle(2);
      gt->Draw("L same"); gl->Draw("L same");
      leg->AddEntry(gt, Form("Template, %s", jerNames[s]), "l");
      leg->AddEntry(gl, Form("Linear < %.0f GeV, %s", linPivot, jerNames[s]), "l");
    }
    leg->Draw();
    d.drawAll({"Pythia8 Jet12 MC"}, {"JER template r04_pileup_EMfracJES", "Dashed: tangent at 10 GeV"}, .45, .88, 22, 650);
    c->SaveAs("pdfs/smear_width_linear.pdf");
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
    auto style = [](TH1D *hh, int col, int ls, int lw) { hh->SetLineColor(col); hh->SetLineStyle(ls); hh->SetLineWidth(lw); hh->SetMarkerSize(0); };
    style(h[0][r], kGray+1, 1, 2);
    style(h[1][r], kBlack, 1, 2);  style(h[2][r], kBlack, 2, 1); style(h[3][r], kBlack, 3, 1);
    style(h[4][r], kRed+1, 1, 2);  style(h[5][r], kRed+1, 2, 1); style(h[6][r], kRed+1, 3, 1);
    for (int sc : {0, 2, 3, 5, 6, 1, 4}) h[sc][r]->Draw("hist same");
    TLine *lf = new TLine(4, 0.5, 4, ymax * 300); lf->SetLineStyle(3); lf->SetLineColor(kGray+2); lf->Draw();
    if (r == 0) {
      d.drawAll({"Pythia8 Jet12 MC"}, {Form("Jet R=%.1f, |#eta_{det}| < %.1f", R, 1.1 - R),
                 Form("p_{T}^{truth,lead} > %.0f GeV", truthLeadMin), Form("p_{T}^{ref}: %s", flav)}, .45, .9, 18, 544);
    }
    if (r == 1) {
      TLegend *leg = new TLegend(.5, .62, .95, .93);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);
      leg->AddEntry(h[0][r], "Unsmeared (p_{T}^{calib})", "l");
      leg->AddEntry(h[1][r], "Template (treemaker)", "l");
      leg->AddEntry(h[4][r], Form("Linear below %.0f GeV", linPivot), "l");
      leg->AddEntry(h[2][r], "Sys. up (dashed)", "l");
      leg->AddEntry(h[3][r], "Sys. down (dotted)", "l");
      leg->Draw();
    }
    TLatex tl; tl.SetNDC(); tl.SetTextSize(0.05);
    tl.DrawLatex(.2, .08, rankNames[r]);

    pb->cd();
    TH1D *rf = (TH1D*)fr->Clone(Form("rf%d", r));
    rf->SetMinimum(0.5); rf->SetMaximum(1.5);
    rf->GetYaxis()->SetTitle("Linear / Template"); rf->GetYaxis()->SetNdivisions(505);
    rf->GetXaxis()->SetLabelSize(0.1); rf->GetXaxis()->SetTitleSize(0.11); rf->GetXaxis()->SetTitleOffset(1.1);
    rf->GetYaxis()->SetLabelSize(0.09); rf->GetYaxis()->SetTitleSize(0.09); rf->GetYaxis()->SetTitleOffset(0.75);
    rf->Draw("axis");
    TLine *l1 = new TLine(0, 1, xmax, 1); l1->SetLineStyle(2); l1->Draw();
    for (int s = 0; s < 3; s++) {
      TH1D *ra = (TH1D*)h[4+s][r]->Clone(Form("ra%d_%d", r, s));
      ra->Divide(h[1+s][r]);
      for (int b = 0; b <= ra->GetNbinsX() + 1; b++)
        if (!std::isfinite(ra->GetBinContent(b)) || h[1+s][r]->GetBinContent(b) <= 0) { ra->SetBinContent(b, 0); ra->SetBinError(b, 0); }
      ra->SetLineColor(s == 0 ? kRed+1 : s == 1 ? kRed-7 : kRed+3); ra->SetLineStyle(s == 0 ? 1 : s == 1 ? 2 : 3);
            ra->SetLineWidth(s == 0 ? 2 : 1);
      ra->Draw("hist same"); // same events and z in both: bin errors would be meaningless
    }
  }
  c->SaveAs(Form("pdfs/jet12_jetpt_R%02d_%s.pdf", radius, flav));
}
