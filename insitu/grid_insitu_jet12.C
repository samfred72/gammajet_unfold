#include "../src/ana.h"
#include "../src/insitu_utility.h"
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <cfloat>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
using namespace std;

// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Non-purity-corrected cross-check of grid_insitu.C: Data region A vs the region A of the
// Jet12_long QCD Pythia8 sample (mostly fake-photon dijets), same scan and conventions. No
// region C or purity here.

const char * insitu_input_dir  = ana::path("insitu/inputs");
const char * insitu_output_dir = ana::path("insitu/output");
const char * insitu_pdf_dir    = ana::path("insitu/pdfs");

// Reported pT bins only (ana::ptBinsUsed).
const int nPtBinsUsed = ana::nPtBinsUsed;
// Shape method: last 3 sparse xJ bins dropped; pT bin 2 merged in pairs.
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;
const int coarseRebinPtBin = 2;
const int coarseGroupSize = 2;
bool shapeMethod = false; // set by grid_insitu_jet12(..., method)

// MC reference sample.
const string mcTrigger = "Jet12_long";

// Jet12 cross-section weight (drawer.h's scalemap); cancels with one sample.
map<int,double> jet_scale = {{12,1.4903e+06}};

// Mean xJ vs pT for MC and raw Data (top); raw and corrected ratios (bottom).
void drawJESPage(TCanvas * c, const char * pdfPath, const char * label, int ir,
    TGraphErrors * gMC, TGraphErrors * gDataRaw, TGraphErrors * gRatioRaw, TGraphErrors * gRatioCorr,
    float pa, float paErrLow, float paErrHigh) {
  c->Clear();
  TPad * p1 = new TPad("p1","",0,.5,1,1);
  TPad * p2 = new TPad("p2","",0,0,1,.5);
  p1->Draw();
  p2->Draw();

  p1->cd();
  p1->SetBottomMargin(0.02);
  p1->SetLeftMargin(.15);
  gPad->SetTicks(1,1);
  TH1F * frame1 = p1->DrawFrame(ana::ptBinsUsed[0], 0.6, ana::ptBinsUsed[nPtBinsUsed], 1.0);
  frame1->GetYaxis()->SetTitle("<x_{J#gamma}>");
  frame1->GetXaxis()->SetLabelSize(0);
  gMC->SetLineColor(kMagenta+1);
  gMC->SetMarkerColor(kMagenta+1);
  gMC->SetMarkerStyle(21);
  gMC->SetLineWidth(2);
  gMC->Draw("p same");
  gDataRaw->SetLineColor(kBlue);
  gDataRaw->SetMarkerColor(kBlue);
  gDataRaw->SetMarkerStyle(20);
  gDataRaw->SetLineWidth(2);
  gDataRaw->Draw("p same");
  TLegend * l1 = new TLegend(.55,.1,.85,.3);
  l1->SetLineWidth(0);
  l1->AddEntry(gMC, "Pythia8 Jet12 (reco, Region A)");
  l1->AddEntry(gDataRaw, "Data (reco, Region A)");
  l1->Draw();
  insitu_utility::drawSPhenixLabel({label}, {
      Form("Jet R=%.1f", ana::JetRs[ir]),
      Form("|#eta^{#gamma}|<%.1f, |#eta^{jet}|<%.1f", ana::photonEtaCut, ana::etacut-ana::JetRs[ir]),
      Form("#Delta#phi>%.0f#pi/%.0f", ana::oppnum, ana::oppden)
    }, .18, .85, 16, p1->GetWh()/1.5);

  p2->cd();
  p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.2);
  p2->SetLeftMargin(.15);
  gPad->SetTicks(1,1);
  TH1F * frame2 = p2->DrawFrame(ana::ptBinsUsed[0], 0.80, ana::ptBinsUsed[nPtBinsUsed], 1.10);
  frame2->GetYaxis()->SetTitle("Data/MC");
  frame2->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  frame2->GetYaxis()->SetTitleSize(0.06);
  frame2->GetYaxis()->SetTitleOffset(1.1);
  frame2->GetYaxis()->SetLabelSize(0.05);
  frame2->GetXaxis()->SetTitleSize(0.06);
  frame2->GetXaxis()->SetLabelSize(0.05);
  gRatioRaw->SetLineColor(kBlack);
  gRatioRaw->SetMarkerColor(kBlack);
  gRatioRaw->SetMarkerStyle(20);
  gRatioRaw->Draw("p same");
  gRatioCorr->SetLineColor(kBlack);
  gRatioCorr->SetMarkerColor(kBlack);
  gRatioCorr->SetMarkerStyle(24); // open circle
  gRatioCorr->Draw("p same");
  TLine * line = new TLine(ana::ptBinsUsed[0],1,ana::ptBinsUsed[nPtBinsUsed],1);
  line->SetLineStyle(9);
  line->Draw("same");
  TLegend * l2 = new TLegend(.55,.7,.85,.9);
  l2->SetLineWidth(0);
  l2->AddEntry(gRatioRaw,  "Raw ratio");
  l2->AddEntry(gRatioCorr, "Corrected ratio");
  l2->Draw();

  float paErr = (paErrLow+paErrHigh)/2.0;
  TLatex jestext;
  jestext.SetNDC();
  jestext.SetTextColor(kRed);
  jestext.DrawLatex(.18,.28, Form("Data to MC (Jet12) JES%s = %.4f #pm %.4f", shapeMethod ? " (shape #chi^{2})" : "", pa, paErr));

  c->SaveAs(pdfPath);
}

// xJ for one pT bin: Jet12 reference, raw Data, Data at the best-fit pa (shape-normalized densities).
void drawXjPage(TCanvas * c, const char * pdfPath, const char * label, int ir, float ptlow, float pthigh,
    TH1D * hMC, TH1D * hDataRaw, TH1D * hDataCorr) {
  c->Clear();
  c->cd();
  gPad->SetLeftMargin(.15);
  gPad->SetTicks();

  TH1D * hMCdisp = (TH1D*)hMC->Clone(Form("%s_disp_mc", hDataCorr->GetName()));
  hMCdisp->Scale(1., "width");
  if (hMCdisp->Integral() > 0) hMCdisp->Scale(1./hMCdisp->Integral());
  TH1D * hRawdisp = (TH1D*)hDataRaw->Clone(Form("%s_disp_raw", hDataCorr->GetName()));
  hRawdisp->Scale(1., "width");
  if (hRawdisp->Integral() > 0) hRawdisp->Scale(1./hRawdisp->Integral());
  TH1D * hCorrdisp = (TH1D*)hDataCorr->Clone(Form("%s_disp", hDataCorr->GetName()));
  hCorrdisp->Scale(1., "width");
  if (hCorrdisp->Integral() > 0) hCorrdisp->Scale(1./hCorrdisp->Integral());

  hMCdisp->SetLineColor(kMagenta+1);
  hMCdisp->SetMarkerColor(kMagenta+1);
  hMCdisp->SetLineWidth(2);
  hMCdisp->GetXaxis()->SetTitle("x_{J#gamma}");
  hMCdisp->GetYaxis()->SetTitle("Shape-normalized counts / bin width");
  hMCdisp->GetYaxis()->SetRangeUser(0, std::max({hMCdisp->GetMaximum(), hRawdisp->GetMaximum(), hCorrdisp->GetMaximum()})*1.4);
  hMCdisp->Draw("hist");

  hRawdisp->SetLineColor(kGray+2);
  hRawdisp->SetMarkerColor(kGray+2);
  hRawdisp->SetMarkerStyle(24); // open circle
  hRawdisp->SetLineWidth(2);
  hRawdisp->Draw("p e same");

  hCorrdisp->SetLineColor(kBlue);
  hCorrdisp->SetMarkerColor(kBlue);
  hCorrdisp->SetMarkerStyle(20);
  hCorrdisp->SetLineWidth(2);
  hCorrdisp->Draw("p e same");

  TLegend * l = new TLegend(.55,.5,.85,.7);
  l->SetLineWidth(0);
  l->AddEntry(hMCdisp,   "Pythia8 Jet12 (reco, Region A)");
  l->AddEntry(hRawdisp,  "Data (raw, Region A)");
  l->AddEntry(hCorrdisp, "Data (JES-corrected, Region A)");
  l->Draw();

  insitu_utility::drawSPhenixLabel({label, Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ptlow, pthigh)}, {
      Form("Jet R=%.1f", ana::JetRs[ir]),
      Form("p_{T}^{jet} > %.0f GeV", ana::jet_calib_pt_cut[ir]),
      Form("|#eta^{#gamma}|<%.1f, |#eta^{jet}|<%.1f", ana::photonEtaCut, ana::etacut-ana::JetRs[ir]),
      Form("#Delta#phi>%.0f#pi/%.0f", ana::oppnum, ana::oppden)
    }, .18, .85, 16, gPad->GetWh());

  c->SaveAs(pdfPath);
}

// method = "mean" or "shape".
void grid_insitu_jet12(string systag = "nominal", string method = "mean") {
  if (method != "mean" && method != "shape") { cout << "method must be \"mean\" or \"shape\"" << endl; return; }
  shapeMethod = (method == "shape");
  const char * tag = shapeMethod ? "jet12_shapechi2" : "jet12";
  // Don't auto-register into the current radius directory.
  TH1::AddDirectory(kFALSE);

  // One file and one PDF per systag; each radius in its ana::rnames[ir] subdirectory.
  string pdfPathStr = Form("%s/grid_insitu_%s_%s.pdf", insitu_pdf_dir, tag, systag.c_str());
  TCanvas * c = new TCanvas("c","",700,700);
  c->SaveAs(Form("%s[", pdfPathStr.c_str()));

  string outfilename = Form("%s/grid_insitu_%s_%s.root", insitu_output_dir, tag, systag.c_str());
  TFile * fout = TFile::Open(outfilename.c_str(), "RECREATE");

  for (int ir = 0; ir < ana::nJetR; ir++) {

  const char * dataFile = Form("%s/Data_%s_insitu.root", insitu_input_dir, systag.c_str());
  vector<DataEvent> dataA = insitu_utility::cacheDataEvents(dataFile, 0, ir);
  cout << "Cached Data events: region A=" << dataA.size() << endl;

  // Low-xJ floor per used pT bin.
  float lowXj[nPtBinsUsed];
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) lowXj[ipt] = insitu_utility::lowXjFloor(ir, ana::ptBinsUsed[ipt]);

  vector<pair<string,double>> mcSamples = {
    {Form("%s/%s_pythia_%s_insitu.root", insitu_input_dir, mcTrigger.c_str(), systag.c_str()), jet_scale[12]},
  };

  float refMean[nPtBinsUsed], refMeanErr[nPtBinsUsed];
  insitu_utility::referenceMeans(mcSamples, 0, ir, refMean, refMeanErr, lowXj);
  vector<vector<double>> refFrac, refFracErr;
  if (shapeMethod) insitu_utility::referenceShape(mcSamples, 0, ir, refFrac, refFracErr, lowXj);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    cout << "Jet12 MC reference <x_J> pt bin " << ipt << " [" << ana::ptBinsUsed[ipt] << "," << ana::ptBinsUsed[ipt+1]
         << "): " << refMean[ipt] << " +/- " << refMeanErr[ipt] << endl;
  }

  // -----------------------------
  // Grid scan over pa (window in insitu_utility.h)
  // -----------------------------
  const int na = insitu_utility::scanN;
  const float lowa = insitu_utility::scanLow, higha = insitu_utility::scanHigh;

  TGraph * gchisq = new TGraph(na);
  gchisq->SetName("gchisq_regionA_jet12ref");
  gchisq->SetTitle(shapeMethod ? ";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});Shape #chi^{2}" : ";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});#chi^{2}");

  float minchisq = FLT_MAX, minpa = 1;
  int ibest = 0;

  for (int ia = 0; ia < na; ia++) {
    float pa = lowa + ia*(higha-lowa)/na;

    float chisq = 0;
    if (shapeMethod) {
      vector<vector<double>> countXj(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
      for (auto & ev : dataA) {
        float x = (ev.jet_pt/pa)/ev.pho_pt;
        if (x < lowXj[ev.ptbin]) continue;
        int ixj = ana::findUnfoldXjBin(x);
        if (ixj < 0 || ixj >= ana::nUnfoldXjBins) continue;
        countXj[ev.ptbin][ixj] += 1;
      }
      for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
        double NA = 0;
        for (double n : countXj[ipt]) NA += n;
        if (NA <= 0) continue;
        vector<double> useA, useRefFrac, useRefFracErr;
        int nBins;
        if (ipt == coarseRebinPtBin) {
          int startBin = ana::findUnfoldXjBin(lowXj[ipt]);
          useA          = insitu_utility::coarsenSum(countXj[ipt], startBin, nXjBinsForChi2, coarseGroupSize);
          useRefFrac    = insitu_utility::coarsenSum(refFrac[ipt], startBin, nXjBinsForChi2, coarseGroupSize);
          useRefFracErr = insitu_utility::coarsenQuadrature(refFracErr[ipt], startBin, nXjBinsForChi2, coarseGroupSize);
          nBins = (int)useA.size();
        } else {
          useA = countXj[ipt];
          useRefFrac = refFrac[ipt]; useRefFracErr = refFracErr[ipt];
          nBins = nXjBinsForChi2;
        }
        double errFloor = 1.0/NA;
        for (int ib = 0; ib < nBins; ib++) {
          double fData = useA[ib]/NA;
          double errData = sqrt(useA[ib])/NA;
          double errt = std::max(sqrt(errData*errData + useRefFracErr[ib]*useRefFracErr[ib]), errFloor);
          double diff = fData - useRefFrac[ib];
          chisq += diff*diff/(errt*errt);
        }
      }
    } else {
    vector<double> sumA(nPtBinsUsed,0), sumA2(nPtBinsUsed,0);
    vector<int> countA(nPtBinsUsed,0);
    for (auto & ev : dataA) {
      float x = (ev.jet_pt/pa)/ev.pho_pt;
      if (x < lowXj[ev.ptbin]) continue;
      sumA[ev.ptbin]  += x;
      sumA2[ev.ptbin] += x*x;
      countA[ev.ptbin]++;
    }
    for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
      if (refMean[ipt] <= 0 || countA[ipt] == 0) continue;
      double mean = sumA[ipt]/countA[ipt];
      double var  = sumA2[ipt]/countA[ipt] - mean*mean;
      double err  = sqrt(std::max(var,0.)/countA[ipt]);
      double diff = 1 - mean/refMean[ipt];
      double errt = sqrt((err*err)/(refMean[ipt]*refMean[ipt])
          + mean*mean*refMeanErr[ipt]*refMeanErr[ipt]/pow(refMean[ipt],4));
      if (errt > 0) chisq += diff*diff/(errt*errt);
    }
    }

    gchisq->SetPoint(ia, pa, chisq);
    if (chisq < minchisq) { minchisq = chisq; minpa = pa; ibest = ia; }
  }

  float errLow, errHigh;
  insitu_utility::findError(gchisq, ibest, minchisq, errLow, errHigh);

  cout << "\nFINAL RESULT (jet R=" << ana::JetRs[ir] << ", systag=" << systag << (shapeMethod ? ", shape chi2" : "") << ", non-purity-corrected, "
       << mcTrigger << " reference)\n";
  cout << "Region A (Data) vs. Region A (" << mcTrigger << " MC):  p_a = " << minpa
       << " +" << errHigh << "/-" << errLow << " (chi2=" << minchisq << ")" << endl;

  // -----------------------------
  // xJ histograms per pT bin: Jet12 reference, raw Data, Data at the best-fit pa
  // -----------------------------
  vector<TH1D*> hxjMC_pt    = insitu_utility::buildMCXjByPtBin(mcSamples, 0, ir, "hxjA_jet12", lowXj);
  vector<TH1D*> hxjA_raw_pt  = insitu_utility::buildXjByPtBin(dataA, 1.0,   nPtBinsUsed, "hxjA_data_raw", lowXj);
  vector<TH1D*> hxjA_corr_pt = insitu_utility::buildXjByPtBin(dataA, minpa, nPtBinsUsed, "hxjA_data_corr", lowXj);

  auto sumPtBins = [&](const vector<TH1D*> & h, const char * name) {
    TH1D * hsum = (TH1D*)h[0]->Clone(name);
    for (int ipt = 1; ipt < nPtBinsUsed; ipt++) hsum->Add(h[ipt]);
    return hsum;
  };
  TH1D * hxjA_jet12          = sumPtBins(hxjMC_pt,     "hxjA_jet12");
  TH1D * hxjA_data_raw       = sumPtBins(hxjA_raw_pt,  "hxjA_data_raw");
  TH1D * hxjA_data_bestscale = sumPtBins(hxjA_corr_pt, "hxjA_data_bestscale");

  // -----------------------------
  // Mean xJ vs pT page
  // -----------------------------
  gStyle->SetOptStat(0);

  float rawMeanA[nPtBinsUsed], rawErrA[nPtBinsUsed];
  float corrMeanA[nPtBinsUsed], corrErrA[nPtBinsUsed];
  insitu_utility::computeRegionAMeans(dataA, 1.0,   rawMeanA,  rawErrA,  lowXj);
  insitu_utility::computeRegionAMeans(dataA, minpa, corrMeanA, corrErrA, lowXj);

  TGraphErrors * gMC          = insitu_utility::meanGraph(refMean, refMeanErr, "gMeanMC_jet12");
  TGraphErrors * gDataRaw     = insitu_utility::meanGraph(rawMeanA,  rawErrA,  "gMeanData_regionA_raw");
  TGraphErrors * gRatioRaw    = insitu_utility::ratioGraph(rawMeanA,  rawErrA,  refMean, refMeanErr, "gRatio_regionA_raw");
  TGraphErrors * gRatioCorr   = insitu_utility::ratioGraph(corrMeanA, corrErrA, refMean, refMeanErr, "gRatio_regionA_corrected");

  const char * pdfPath = pdfPathStr.c_str();
  drawJESPage(c, pdfPath, Form("Region A, non-purity-corrected (%s ref.)", mcTrigger.c_str()), ir, gMC, gDataRaw, gRatioRaw, gRatioCorr, minpa, errLow, errHigh);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    drawXjPage(c, pdfPath, Form("Region A, non-purity-corrected (%s ref.)", mcTrigger.c_str()), ir, ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1],
        hxjMC_pt[ipt], hxjA_raw_pt[ipt], hxjA_corr_pt[ipt]);
  }

  // -----------------------------
  // Save
  // -----------------------------
  fout->cd();
  fout->mkdir(ana::rnames[ir])->cd();
  gchisq->Write();
  hxjA_jet12->Write();
  hxjA_data_raw->Write();
  hxjA_data_bestscale->Write();

  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    hxjMC_pt[ipt]->Write();
    hxjA_raw_pt[ipt]->Write();
    hxjA_corr_pt[ipt]->Write();
  }

  gMC->Write();
  gDataRaw->Write();
  gRatioRaw->Write();
  gRatioCorr->Write();

  TTree * wt = new TTree("results", shapeMethod ? "best-fit jet energy scale result (shape chi2, non-purity-corrected, Jet12 reference)" : "best-fit jet energy scale result (non-purity-corrected, Jet12 reference)");
  float wpa = minpa, wchisq = minchisq, werrLow = errLow, werrHigh = errHigh;
  wt->Branch("pa_regionA_jet12ref", &wpa);
  wt->Branch("chisq_regionA_jet12ref", &wchisq);
  wt->Branch("errLow_regionA_jet12ref", &werrLow);
  wt->Branch("errHigh_regionA_jet12ref", &werrHigh);
  wt->Fill();
  wt->Write();

  cout << "Finished ir=" << ir << " (" << ana::rnames[ir] << ")" << endl;
  } // end of ir loop

  c->SaveAs(Form("%s]", pdfPathStr.c_str()));
  cout << "Wrote " << pdfPathStr << endl;
  fout->Close();
  cout << "Wrote " << outfilename << endl;
}
