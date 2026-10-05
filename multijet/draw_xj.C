#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMath.h"

// Mean multijet balance x_j = pT,lead/|pT,sub + pT,subsub| vs leading-jet pT, Data vs MC, from
// analysis.cc's multijet_analysis_<sim>.root (sys = RECO/HIGH/LOW).
void draw_xj(int radius = 4, const char * sys = "RECO", const char * sim = "pythia") {

  gStyle->SetOptStat(0);

  // ---------------------------------------------------------
  // Open file
  // ---------------------------------------------------------

  TFile *f = TFile::Open(Form("multijet_analysis_%s.root", sim));

  if (!f || f->IsZombie()) {
    std::cout << "Could not open file." << std::endl;
    return;
  }

  // ---------------------------------------------------------
  // Histograms for mean values
  // ---------------------------------------------------------

  const int nPtBins = 7;
  float ptBins[nPtBins + 1] = {20, 25, 30, 35, 40, 50, 60, 70};
  TH1D *hDataMean = new TH1D(
      "hDataMean",
      ";p_{T}^{lead} [GeV];Mean x_{j}",
      nPtBins,
      ptBins);

  TH1D *hSimMean = new TH1D(
      "hSimMean",
      ";p_{T}^{lead} [GeV];Mean x_{j}",
      nPtBins,
      ptBins);

  TH1D *hRatio = new TH1D(
      "hRatio",
      ";p_{T}^{lead} [GeV];Data / Sim",
      nPtBins,
      ptBins);

  // ---------------------------------------------------------
  // Loop over pT bins
  // ---------------------------------------------------------

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

  TCanvas *c = new TCanvas("c","c",700,700);

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

  c->SaveAs(Form("pdfs/hxj_means_r%d_%s_%s.pdf", radius, sys, sim));
}
