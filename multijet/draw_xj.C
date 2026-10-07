#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLatex.h"
#include <algorithm>
#include <cmath>

// Mean multijet balance x_j = pT,lead/|pT,sub + pT,subsub| vs leading-jet pT, Data vs MC, from
// analysis.cc's multijet_analysis_<sim>.root (sys = RECO/HIGH/LOW). A second PDF overlays the
// Data and MC x_j distributions per leading-pT bin.
// tag: output-name suffix of the analysis run (e.g. "_tight" for analysis --tight).
// meansPdf/distsPdf: print the page into this multi-page PDF ("file.pdf(" opens it, "file.pdf)"
// closes it) instead of writing one file per radius and JER variation; "-" skips that page
// (see draw_xj_all below).
void draw_xj(int radius = 4, const char * sys = "RECO", const char * sim = "pythia", const char * tag = "",
             const char * meansPdf = "", const char * distsPdf = "") {

  gStyle->SetOptStat(0);

  // ---------------------------------------------------------
  // Open file
  // ---------------------------------------------------------

  TFile *f = TFile::Open(Form("multijet_analysis_%s%s.root", sim, tag));

  if (!f || f->IsZombie()) {
    std::cout << "Could not open file." << std::endl;
    return;
  }

  // ---------------------------------------------------------
  // Histograms for mean values
  // ---------------------------------------------------------

  const int nPtBins = 4;
  float ptBins[nPtBins + 1] = {20, 25, 30, 35, 50};
  TH1D *hDataMean = new TH1D(
      Form("hDataMean_r%d_%s", radius, sys),
      ";p_{T}^{lead} [GeV];Mean x_{j}",
      nPtBins,
      ptBins);

  TH1D *hSimMean = new TH1D(
      Form("hSimMean_r%d_%s", radius, sys),
      ";p_{T}^{lead} [GeV];Mean x_{j}",
      nPtBins,
      ptBins);

  TH1D *hRatio = new TH1D(
      Form("hRatio_r%d_%s", radius, sys),
      ";p_{T}^{lead} [GeV];Data / Sim",
      nPtBins,
      ptBins);

  // ---------------------------------------------------------
  // Loop over pT bins
  // ---------------------------------------------------------

  TH1D *hDist[nPtBins][2] = {{nullptr}}; // [pT bin][0 data, 1 sim]

  for (int i = 0; i < nPtBins; i++) {

    // pT bin labels: 20,25,30,35,40,50,60 ...
    int pt = (int)ptBins[i];

    TString hnameData = Form("hxj_r%d_%i_data", radius, pt); // data has no JER variations
    TString hnameSim  = Form("hxj_r%d_%i_%s_sim", radius, pt, sys);

    TH1D *hData = (TH1D*)f->Get(hnameData);
    TH1D *hSim  = (TH1D*)f->Get(hnameSim);

    if (!hData || !hSim) {
      std::cout << "Missing histogram: "
                << hnameData << " or "
                << hnameSim << std::endl;
      continue;
    }

    hDist[i][0] = (TH1D*)hData->Clone(Form("hdist_data_%d_r%d_%s", i, radius, sys));
    hDist[i][1] = (TH1D*)hSim->Clone(Form("hdist_sim_%d_r%d_%s", i, radius, sys));

    double meanData = hData->GetMean();
    double errData  = hData->GetMeanError();

    double meanSim = hSim->GetMean();
    double errSim  = hSim->GetMeanError();

    cout << meanData << " " << meanSim << " ";

    hDataMean->SetBinContent(i+1, meanData);
    hDataMean->SetBinError(i+1, errData);

    hSimMean->SetBinContent(i+1, meanSim);
    hSimMean->SetBinError(i+1, errSim);

    // Ratio + propagated uncertainty
    double ratio = 0.0;
    double ratioErr = 0.0;

    if (meanSim != 0.0) {

      ratio = meanData / meanSim;
      cout << ratio << endl;

      double relErrData =
          (meanData != 0.0) ? errData / meanData : 0.0;

      double relErrSim = errSim / meanSim;

      ratioErr = ratio * std::sqrt(
          relErrData*relErrData +
          relErrSim*relErrSim);
    }

    hRatio->SetBinContent(i+1, ratio);
    hRatio->SetBinError(i+1, ratioErr);
  }

  // ---------------------------------------------------------
  // Style
  // ---------------------------------------------------------

  hDataMean->SetLineColor(kBlack);
  hDataMean->SetMarkerColor(kBlack);
  hDataMean->SetMarkerStyle(20);
  hDataMean->SetMarkerSize(1.2);

  hSimMean->SetLineColor(kRed+1);
  hSimMean->SetMarkerColor(kRed+1);
  hSimMean->SetMarkerStyle(21);
  hSimMean->SetMarkerSize(1.2);

  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerColor(kBlack);
  hRatio->SetMarkerStyle(20);
  hRatio->SetMarkerSize(1.2);

  // ---------------------------------------------------------
  // Canvas + pads
  // ---------------------------------------------------------

  TCanvas *c = new TCanvas(Form("c_r%d_%s", radius, sys),"c",700,700);

  TPad *pTop = new TPad(
      "pTop","pTop",
      0.0,0.5,1.0,1.0);

  TPad *pBot = new TPad(
      "pBot","pBot",
      0.0,0.0,1.0,0.5);

  pTop->SetBottomMargin(0.02);

  pBot->SetTopMargin(0.02);
  pBot->SetBottomMargin(0.15);

  pTop->Draw();
  pBot->Draw();

  // ---------------------------------------------------------
  // Top panel
  // ---------------------------------------------------------

  pTop->cd();

  hDataMean->SetTitle("");

  hDataMean->GetYaxis()->SetTitle("<x_{j}>");
  hDataMean->GetYaxis()->SetTitleSize(0.05);
  hDataMean->GetYaxis()->SetLabelSize(0.04);

  hDataMean->GetXaxis()->SetLabelSize(0);

  hDataMean->Draw("E1");
  hSimMean->Draw("E1 SAME");

  TLegend *leg = new TLegend(
      0.62,0.72,0.88,0.88);

  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  leg->AddEntry(hDataMean,"Data","pl");
  leg->AddEntry(hSimMean,"Simulation","pl");

  leg->Draw();

  TLatex label;
  label.SetNDC();
  label.SetTextSize(0.05);
  label.DrawLatex(0.15, 0.82, Form("Jet R = 0.%d, JER %s, %s", radius, sys, sim));

  // ---------------------------------------------------------
  // Bottom panel
  // ---------------------------------------------------------

  pBot->cd();

  hRatio->SetTitle("");

  hRatio->GetYaxis()->SetTitle("Data/Sim");
  hRatio->GetYaxis()->SetTitleSize(0.06);
  hRatio->GetYaxis()->SetTitleOffset(0.8);
  hRatio->GetYaxis()->SetLabelSize(0.05);

  hRatio->GetXaxis()->SetTitle("p_{T}^{lead} [GeV]");
  hRatio->GetXaxis()->SetTitleSize(0.06);
  hRatio->GetXaxis()->SetLabelSize(0.05);

  hRatio->SetMinimum(0.95);
  hRatio->SetMaximum(1.1);

  hRatio->Draw("E1");

  TLine *line = new TLine(
      ptBins[0],
      1.0,
      ptBins[nPtBins],
      1.0);

  line->SetLineStyle(2);
  line->Draw("SAME");

  if (strcmp(meansPdf, "-")) c->Print(meansPdf[0] ? meansPdf : Form("pdfs/hxj_means_r%d_%s_%s%s.pdf", radius, sys, sim, tag));

  // ---------------------------------------------------------
  // Data and MC x_j distributions per leading-pT bin, each normalized to unit area
  // ---------------------------------------------------------

  TCanvas *cDist = new TCanvas(Form("cDist_r%d_%s", radius, sys),"cDist",1600,800);
  cDist->Divide(3,2);
  TLatex tex;
  tex.SetNDC();
  for (int i = 0; i < nPtBins; i++) {
    if (!hDist[i][0]) continue;
    TVirtualPad *p = cDist->cd(i+1);
    p->SetLeftMargin(0.15);
    p->SetRightMargin(0.04);
    p->SetTopMargin(0.06);
    p->SetBottomMargin(0.13);

    double ymax = 0;
    for (int k = 0; k < 2; k++) {
      TH1D *h = hDist[i][k];
      for (int b = 0; b <= h->GetNbinsX()+1; b++) {
        if (!std::isfinite(h->GetBinContent(b)) || !std::isfinite(h->GetBinError(b))) {
          h->SetBinContent(b, 0);
          h->SetBinError(b, 0);
        }
      }
      if (h->Integral() > 0) h->Scale(1.0/h->Integral());
      ymax = std::max(ymax, h->GetMaximum());
    }

    TH1D *hD = hDist[i][0];
    TH1D *hS = hDist[i][1];
    hS->SetTitle(";x_{j};Normalized counts");
    hS->SetLineColor(kRed+1);
    hS->SetLineWidth(2);
    hS->SetMinimum(0);
    hS->SetMaximum(1.45*ymax);
    hS->GetXaxis()->SetTitleSize(0.055);
    hS->GetXaxis()->SetLabelSize(0.045);
    hS->GetYaxis()->SetTitleSize(0.055);
    hS->GetYaxis()->SetLabelSize(0.045);
    hS->GetYaxis()->SetTitleOffset(1.3);
    hD->SetLineColor(kBlack);
    hD->SetMarkerColor(kBlack);
    hD->SetMarkerStyle(20);
    hD->SetMarkerSize(0.7);

    hS->Draw("HIST");
    hD->Draw("E1 SAME");
    tex.SetTextSize(0.055);
    tex.DrawLatex(0.45, 0.86, Form("%.0f < p_{T}^{lead} < %.0f GeV", ptBins[i], ptBins[i+1]));
  }

  cDist->cd(6);
  TLegend *legDist = new TLegend(0.1, 0.35, 0.9, 0.75);
  legDist->SetBorderSize(0);
  legDist->SetFillStyle(0);
  legDist->SetTextSize(0.07);
  legDist->SetHeader(Form("Jet R = 0.%d, JER %s", radius, sys));
  for (int i = 0; i < nPtBins; i++) {
    if (!hDist[i][0]) continue;
    legDist->AddEntry(hDist[i][0], "Data", "pl");
    legDist->AddEntry(hDist[i][1], sim, "l");
    break;
  }
  legDist->Draw();

  if (strcmp(distsPdf, "-")) cDist->Print(distsPdf[0] ? distsPdf : Form("pdfs/hxj_dists_r%d_%s_%s%s.pdf", radius, sys, sim, tag));
  delete c;
  delete cDist;
}

// Every radius and JER variation in two multi-page PDFs, pdfs/hxj_means_<sim><tag>.pdf and
// pdfs/hxj_dists_<sim><tag>.pdf. Run: root -b -l -q -e '.L draw_xj.C' -e 'draw_xj_all("pythia")'
// The two files are written one after the other: ROOT keeps only one multi-page PDF open at a time.
void draw_xj_all(const char * sim = "pythia", const char * tag = "") {
  const char * systs[3] = {"RECO", "HIGH", "LOW"};
  for (int pass = 0; pass < 2; pass++) {
    TString out = Form("pdfs/hxj_%s_%s%s.pdf", pass ? "dists" : "means", sim, tag);
    for (int r = 2; r <= 8; r++) {
      for (int j = 0; j < 3; j++) {
        TString page = out + (r == 2 && j == 0 ? "(" : (r == 8 && j == 2 ? ")" : ""));
        draw_xj(r, systs[j], sim, tag, pass ? "-" : page.Data(), pass ? page.Data() : "-");
      }
    }
  }
}
