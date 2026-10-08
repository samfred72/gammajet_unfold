#include "../../src/ana.h"
#include "../../src/drawer.h"
// Explicit load - see drawing/draw_final_result.C. Run interpreted, never with ACLiC "+".
R__LOAD_LIBRARY(libgammajet_unfold.so);

// The JER smearing width the treemaker applies to MC jets (gammajet treemaking,
// CaloAna.cc:378-384 and smear_pt() at :908): smeared pT = Gaus(pt_calib, pt_ref * w(pt_ref)),
// w = h_jer_smear_r04_pileup_EMfracJES_{nominal,sysup,sysdown}->Interpolate(pt_ref).
// The same R=0.4 template is used for all radii.
//
// The templates (jerband_smearing_templates.root, copied from treemaking/macros/) are
// binned on [5, 80] GeV with empty underflow. TH1::Interpolate does not extrapolate:
// below the first bin centre (5.1875 GeV) it returns the first bin's content, so w is
// flat there. The curves below are w evaluated exactly as CaloAna does, through
// Interpolate, so the flat region is what the treemaker actually uses.
//
//   left   fractional width w = sigma/pT
//   right  absolute width sigma = pT * w (GeV)

void draw_smear_function() {
  gStyle->SetOptStat(0);
  drawer d("pythia", "nominal");

  TFile f("jerband_smearing_templates.root");
  const char *names[3]  = {"h_jer_smear_r04_pileup_EMfracJES_nominal",
                           "h_jer_smear_r04_pileup_EMfracJES_sysup",
                           "h_jer_smear_r04_pileup_EMfracJES_sysdown"};
  const char *labels[3] = {"Nominal", "Sys. up", "Sys. down"};
  const int colors[3]   = {kBlack, kRed+1, kBlue+1};
  const int styles[3]   = {1, 2, 2};

  const double ptMin = 0, ptMax = 20;
  const int nPts = 401;
  TH1D *h[3];
  TGraph *gFrac[3], *gAbs[3];
  for (int i = 0; i < 3; i++) {
    h[i] = (TH1D*)f.Get(names[i]);
    if (!h[i]) { printf("missing %s\n", names[i]); return; }
    gFrac[i] = new TGraph(nPts);
    gAbs[i]  = new TGraph(nPts);
    for (int k = 0; k < nPts; k++) {
      double pt = ptMin + (ptMax - ptMin) * k / (nPts - 1);
      double w = h[i]->Interpolate(pt);
      gFrac[i]->SetPoint(k, pt, w);
      gAbs[i] ->SetPoint(k, pt, pt * w);
    }
    for (TGraph *g : {gFrac[i], gAbs[i]}) {
      g->SetLineColor(colors[i]);
      g->SetLineStyle(styles[i]);
      g->SetLineWidth(2);
    }
  }
  const double firstCentre = h[0]->GetXaxis()->GetBinCenter(1);
  const double tmplLow = h[0]->GetXaxis()->GetXmin();

  TCanvas *c = new TCanvas("c", "", 1400, 650);
  c->Divide(2, 1);

  auto drawPad = [&](int ipad, TGraph **g, const char *ytitle, double ymax, bool label) {
    TPad *p = (TPad*)c->cd(ipad);
    p->SetLeftMargin(.14); p->SetRightMargin(.04);
    p->SetBottomMargin(.12); p->SetTopMargin(.05);
    TH1F *fr = p->DrawFrame(ptMin, 0, ptMax, ymax);
    fr->GetXaxis()->SetTitle("p_{T}^{ref} [GeV]");
    fr->GetYaxis()->SetTitle(ytitle);
    fr->GetXaxis()->SetTitleSize(0.045); fr->GetYaxis()->SetTitleSize(0.045);
    fr->GetYaxis()->SetTitleOffset(1.4);

    // template's lower edge and the point below which Interpolate returns a constant
    TBox *box = new TBox(ptMin, 0, firstCentre, ymax);
    box->SetFillColor(kGray); box->SetFillStyle(1001);
    box->Draw();
    TLine *l = new TLine(tmplLow, 0, tmplLow, ymax);
    l->SetLineStyle(3); l->SetLineColor(kGray+2);
    l->Draw();
    p->RedrawAxis();

    for (int i = 0; i < 3; i++) g[i]->Draw("L same");

    if (!label) {
      TLatex t; t.SetNDC(); t.SetTextSize(0.032); t.SetTextColor(kGray+3);
      t.DrawLatex(.40, .88, Form("Grey: p_{T}^{ref} < %.4g GeV (first bin centre),", firstCentre));
      t.DrawLatex(.40, .83, "Interpolate returns a constant #sigma/p_{T}");
      t.DrawLatex(.40, .78, Form("Dotted: template lower edge (%.0f GeV)", tmplLow));
    } else {
      d.drawAll({"Pythia8 #gamma+jet MC"},
                {"JER template r04_pileup_EMfracJES", "R=0.4 template, applied to all R"},
                .42, .88, 22, 650);
      TLegend *leg = new TLegend(.55, .5, .9, .67);
      leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.04);
      for (int i = 0; i < 3; i++) leg->AddEntry(g[i], labels[i], "l");
      leg->Draw();
    }
  };

  drawPad(1, gFrac, "#sigma / p_{T}", 0.9, true);
  drawPad(2, gAbs,  "#sigma [GeV]",   4.5, false);

  c->SaveAs("pdfs/smear_function.pdf");

  printf("pT   nominal  up      down   (sigma/pT)\n");
  for (double pt : {1., 2., 3., 4., 5., 5.1875, 6., 8., 10., 15., 20.})
    printf("%5.2f  %.4f  %.4f  %.4f\n", pt, h[0]->Interpolate(pt), h[1]->Interpolate(pt), h[2]->Interpolate(pt));
}
