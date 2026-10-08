#include "../../src/ana.h"
#include "../../src/drawer.h"
// Explicit load - see drawing/draw_final_result.C. Run interpreted, never with ACLiC "+".
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Drawing half of fill_jet12_xj.C (see its header for the width functions, the two smearing
// configurations and the selection). One PDF per configuration:
//   pdfs/jet12_xj_perjet_R<RR>.pdf   every jet smeared on its own
//   pdfs/jet12_xj_truth2_R<RR>.pdf   nominal: jets 2+3 belonging to one truth jet smeared as one recoil
// Columns: leading-jet pT bins 20-25, 25-30, 30-35 GeV. Top: unit-normalised xJ; bottom: ratio to the
// template smear (same events and deviates in all, so the ratio is drawn without bin errors).

void draw_jet12_xj(int radius = 4) {
  gStyle->SetOptStat(0);
  drawer d("pythia", "nominal");
  const float R = radius / 10.0f;
  const int nCfg = 2, nW = 4, nPt = 3;
  const char *cfgNames[nCfg] = {"perjet", "truth2"};
  const char *cfgLabels[nCfg] = {"Per-jet smearing only", "Nominal: truth2 recoil smearing"};
  const char *wNames[nW] = {"calib", "tmpl", "lin15", "lin5"};
  const char *wLabels[nW] = {"Unsmeared", "Template", "Tangent at 15 GeV", "Slope at 5 GeV"};
  const int col[nW] = {kGray+1, kBlack, kAzure+1, kRed+1};
  const float ptBins[nPt+1] = {20, 25, 30, 35};

  TFile *fin = TFile::Open(Form("jet12_xj_R%02d.root", radius));
  const double truthLeadMin = ((TParameter<double>*)fin->Get("truthLeadMin"))->GetVal();
  TH1D *h[nCfg][nW][nPt];
  for (int c = 0; c < nCfg; c++) for (int w = 0; w < nW; w++) for (int k = 0; k < nPt; k++) {
    h[c][w][k] = (TH1D*)fin->Get(Form("hxj_%s_%s_%.0f", cfgNames[c], wNames[w], ptBins[k]));
    h[c][w][k]->SetDirectory(0);
    for (int b = 0; b <= h[c][w][k]->GetNbinsX() + 1; b++)
      if (!std::isfinite(h[c][w][k]->GetBinContent(b))) { h[c][w][k]->SetBinContent(b, 0); h[c][w][k]->SetBinError(b, 0); }
    if (h[c][w][k]->Integral() > 0) h[c][w][k]->Scale(1.0 / h[c][w][k]->Integral(), "width");
  }
  fin->Close();

  printf("<xJ> (unit-normalised histograms)\n");
  for (int c = 0; c < nCfg; c++) for (int k = 0; k < nPt; k++) {
    printf("%-7s %.0f-%.0f:", cfgNames[c], ptBins[k], ptBins[k+1]);
    for (int w = 0; w < nW; w++) printf("  %s %.4f", wNames[w], h[c][w][k]->GetMean());
    printf("   lin15-tmpl %+.4f  lin5-tmpl %+.4f\n", h[c][2][k]->GetMean() - h[c][1][k]->GetMean(), h[c][3][k]->GetMean() - h[c][1][k]->GetMean());
  }

  for (int c = 0; c < nCfg; c++) {
    TCanvas *cv = new TCanvas(Form("c%d", c), "", 1800, 800);
    for (int k = 0; k < nPt; k++) {
      const double x0 = k / 3.0, x1 = (k + 1) / 3.0;
      TPad *pt = new TPad(Form("pt%d_%d", c, k), "", x0, 0.32, x1, 1);
      TPad *pb = new TPad(Form("pb%d_%d", c, k), "", x0, 0, x1, 0.32);
      for (TPad *p : {pt, pb}) { p->SetLeftMargin(.15); p->SetRightMargin(.04); }
      pt->SetBottomMargin(0.02); pt->SetTopMargin(.05); pb->SetTopMargin(0.02); pb->SetBottomMargin(.3);
      cv->cd(); pt->Draw(); pb->Draw();

      pt->cd();
      double ymax = 0;
      for (int w = 0; w < nW; w++) ymax = std::max(ymax, h[c][w][k]->GetMaximum());
      TH1D *fr = (TH1D*)h[c][1][k]->Clone(Form("fr%d_%d", c, k)); fr->Reset();
      fr->GetXaxis()->SetRangeUser(0.4, 2.65); fr->SetMinimum(0); fr->SetMaximum(ymax * 1.6);
      fr->GetYaxis()->SetTitle("#frac{1}{N} #frac{dN}{dx_{J}}");
      fr->GetXaxis()->SetLabelSize(0); fr->GetYaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetLabelSize(0.045);
      fr->GetYaxis()->SetTitleOffset(1.4);
      fr->Draw("axis");
      for (int w : {0, 1, 2, 3}) {
        h[c][w][k]->SetLineColor(col[w]); h[c][w][k]->SetMarkerColor(col[w]);
        h[c][w][k]->SetLineWidth(2); h[c][w][k]->SetMarkerSize(0);
        h[c][w][k]->Draw("hist same");
      }
      if (k == 0)
        d.drawAll({"Pythia8 Jet12 MC"}, {Form("Jet R=%.1f, %.0f < p_{T,1} < %.0f GeV", R, ptBins[k], ptBins[k+1]),
                   Form("p_{T}^{truth,lead} > %.0f GeV", truthLeadMin), cfgLabels[c]}, .45, .9, 18, 544);
      else
        d.drawAll({"Pythia8 Jet12 MC"}, {Form("Jet R=%.1f, %.0f < p_{T,1} < %.0f GeV", R, ptBins[k], ptBins[k+1])}, .45, .9, 18, 544);
      if (k == 1) {
        TLegend *leg = new TLegend(.42, .4, .97, .7);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.035);
        leg->SetTextSize(0.032);
        for (int w = 0; w < nW; w++) leg->AddEntry(h[c][w][k], Form("%s, #LTx_{J}#GT = %.3f", wLabels[w], h[c][w][k]->GetMean()), "l");
        leg->Draw();
      } else {
        TLatex tl; tl.SetNDC(); tl.SetTextSize(0.035);
        for (int w = 0; w < nW; w++) { tl.SetTextColor(col[w]); tl.DrawLatex(.55, .6 - 0.05*w, Form("#LTx_{J}#GT = %.3f", h[c][w][k]->GetMean())); }
      }

      pb->cd();
      TH1D *rf = (TH1D*)fr->Clone(Form("rf%d_%d", c, k));
      rf->SetMinimum(0.5); rf->SetMaximum(1.5);
      rf->GetXaxis()->SetTitle("x_{J} = p_{T,1} / |p_{T,2} + p_{T,3}|");
      rf->GetYaxis()->SetTitle("Ratio to template"); rf->GetYaxis()->SetNdivisions(505);
      rf->GetXaxis()->SetLabelSize(0.1); rf->GetXaxis()->SetTitleSize(0.1); rf->GetXaxis()->SetTitleOffset(1.2);
      rf->GetYaxis()->SetLabelSize(0.09); rf->GetYaxis()->SetTitleSize(0.09); rf->GetYaxis()->SetTitleOffset(0.75);
      rf->Draw("axis");
      TLine *l1 = new TLine(0.4, 1, 2.65, 1); l1->SetLineStyle(2); l1->Draw();
      for (int w : {0, 2, 3}) {
        TH1D *ra = (TH1D*)h[c][w][k]->Clone(Form("ra%d_%d_%d", c, w, k));
        ra->Divide(h[c][1][k]);
        for (int b = 0; b <= ra->GetNbinsX() + 1; b++)
          if (!std::isfinite(ra->GetBinContent(b)) || h[c][1][k]->GetBinContent(b) <= 0) { ra->SetBinContent(b, 0); ra->SetBinError(b, 0); }
        ra->SetLineColor(col[w]); ra->SetLineWidth(2);
        ra->Draw("hist same");
      }
    }
    cv->SaveAs(Form("pdfs/jet12_xj_%s_R%02d.pdf", cfgNames[c], radius));
  }
}
