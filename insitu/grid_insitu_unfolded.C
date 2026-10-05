#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/insitu_utility.h"
#include "../src/unfold_utility.h"
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <cfloat>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
using namespace std;

// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);
R__LOAD_LIBRARY(libRooUnfold.so);

// Unfolded version of grid_insitu.C's purity-corrected scan: at each trial pa, unfold the
// purity-corrected Data through the nominal response (fixed) and compare its mean xJ (or xJ
// shape) with the Pythia8 truth, per photon-pT bin. Purity is fixed across the scan.

const char * insitu_input_dir  = ana::path("insitu/inputs");
const char * insitu_output_dir = ana::path("insitu/output");
const char * insitu_pdf_dir    = ana::path("insitu/pdfs");

const int nPtBinsUsed = ana::nPtBinsUsed;
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3; // shape method: last 3 sparse x_J bins dropped
bool shapeMethod = false; // set by grid_insitu_unfolded(..., method)

// Same iteration count as draw_purity_corrected.C.
const int niterate = 2;

// Keep every ana::ptBins bin (restrictToUsed=false): the unfolding needs the full flattened
// input for cross-bin migration.

// Purity-correct all pT bins, unfold, and return the mean xJ and error per used pT bin.
// unfoldedOut (if non-null) gets clones of the unfolded histograms (caller owns them).
void computeUnfoldedMeans(const vector<DataEvent> & dataA, const vector<DataEvent> & dataC, float pa,
    const float purity[], const float purityC[], RooUnfoldResponse * response, TH1D * respRecoTemplate,
    float mean[], float err[], const float lowXj[], vector<TH1D*> * unfoldedOut = nullptr) {
  vector<TH1D*> hA    = insitu_utility::buildXjByPtBin(dataA, pa, ana::nPtBins, "hUnfA_tmp", lowXj);
  vector<TH1D*> hC    = insitu_utility::buildXjByPtBin(dataC, pa, ana::nPtBins, "hUnfC_tmp", lowXj);
  // Purity fixed: zero error arrays.
  float zeroErr[ana::nPtBins] = {0};
  vector<TH1D*> hCorr = insitu_utility::purityCorrectByPtBin(hA, hC, ana::nPtBins,
      purity, zeroErr, zeroErr, purityC, zeroErr, zeroErr, "hUnfCorr_tmp");

  TH1D * flatCorrected = (TH1D*)respRecoTemplate->Clone("flatCorrected_tmp");
  flatCorrected->Reset("ICES");
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) unfold_utility::reflattenXj(hCorr[ipt], ipt, flatCorrected);

  // includeSystematics=false: called at every scan point (~3 s each otherwise). The errors keep
  // the data-statistics term only; this is a cross-check, not the headline JES.
  TH1D * flatUnfolded = unfold_utility::unfoldOnce(response, flatCorrected, niterate, "flatUnfolded_tmp", false);

  for (int k = 0; k < nPtBinsUsed; k++) {
    int ipt = ana::firstUsedPtBin + k;
    TH1D * hU = unfold_utility::unflattenXj(flatUnfolded, ipt, "hUnfoldedPt_tmp");
    mean[k] = hU->GetMean();
    err[k]  = hU->GetMeanError();
    if (unfoldedOut) (*unfoldedOut)[k] = (TH1D*)hU->Clone(Form("hUnfoldSave_pt%d", k));
    delete hU;
  }

  for (auto h : hA)    delete h;
  for (auto h : hC)    delete h;
  for (auto h : hCorr) delete h;
  delete flatCorrected;
  delete flatUnfolded;
}

// Shape method: chi2 of unfolded vs truth xJ bin fractions per used pT bin.
float computeUnfoldedShapeChi2(const vector<DataEvent> & dataA, const vector<DataEvent> & dataC, float pa,
    const float purity[], const float purityC[], RooUnfoldResponse * response, TH1D * respRecoTemplate,
    const vector<vector<double>> & truthFrac, const vector<vector<double>> & truthFracErr, const float lowXj[]) {
  vector<TH1D*> hA    = insitu_utility::buildXjByPtBin(dataA, pa, ana::nPtBins, "hUnfShapeA_tmp", lowXj);
  vector<TH1D*> hC    = insitu_utility::buildXjByPtBin(dataC, pa, ana::nPtBins, "hUnfShapeC_tmp", lowXj);
  float zeroErrShape[ana::nPtBins] = {0};
  vector<TH1D*> hCorr = insitu_utility::purityCorrectByPtBin(hA, hC, ana::nPtBins,
      purity, zeroErrShape, zeroErrShape, purityC, zeroErrShape, zeroErrShape, "hUnfShapeCorr_tmp");

  TH1D * flatCorrected = (TH1D*)respRecoTemplate->Clone("flatCorrectedShape_tmp");
  flatCorrected->Reset("ICES");
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) unfold_utility::reflattenXj(hCorr[ipt], ipt, flatCorrected);

  // includeSystematics=false, as above.
  TH1D * flatUnfolded = unfold_utility::unfoldOnce(response, flatCorrected, niterate, "flatUnfoldedShape_tmp", false);

  float chisq = 0;
  for (int k = 0; k < nPtBinsUsed; k++) {
    int ipt = ana::firstUsedPtBin + k;
    TH1D * hU = unfold_utility::unflattenXj(flatUnfolded, ipt, "hUnfoldedShapePt_tmp");
    double N = hU->Integral();
    if (N > 0) {
      // errt floored at 1/N, as in grid_insitu.C.
      double errFloor = 1.0/N;
      for (int ixj = 1; ixj <= nXjBinsForChi2; ixj++) {
        double fUnf   = hU->GetBinContent(ixj)/N;
        double errUnf = hU->GetBinError(ixj)/N;
        double errt = sqrt(errUnf*errUnf + truthFracErr[k][ixj-1]*truthFracErr[k][ixj-1]);
        errt = std::max(errt, errFloor);
        double diff = fUnf - truthFrac[k][ixj-1];
        chisq += diff*diff/(errt*errt);
      }
    }
    delete hU;
  }

  for (auto h : hA)    delete h;
  for (auto h : hC)    delete h;
  for (auto h : hCorr) delete h;
  delete flatCorrected;
  delete flatUnfolded;
  return chisq;
}

// Unfolded mean xJ vs pT for truth and raw Data (top); raw and corrected ratios (bottom).
void drawJESPage(TCanvas * c, const char * pdfPath, const char * label, int ir,
    TGraphErrors * gTruth, TGraphErrors * gUnfoldRaw, TGraphErrors * gRatioRaw, TGraphErrors * gRatioCorr,
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
  gTruth->SetLineColor(kBlack);
  gTruth->SetMarkerColor(kBlack);
  gTruth->SetMarkerStyle(20);
  gTruth->SetLineWidth(2);
  gTruth->Draw("p same");
  gUnfoldRaw->SetLineColor(kRed);
  gUnfoldRaw->SetMarkerColor(kRed);
  gUnfoldRaw->SetMarkerStyle(21);
  gUnfoldRaw->SetLineWidth(2);
  gUnfoldRaw->Draw("p same");
  TLegend * l1 = new TLegend(.55,.1,.85,.3);
  l1->SetLineWidth(0);
  l1->AddEntry(gTruth, "Pythia8 #gamma+jet (truth)");
  l1->AddEntry(gUnfoldRaw, "Data (unfolded, raw JES)");
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
  TH1F * frame2 = p2->DrawFrame(ana::ptBinsUsed[0], 0.80, ana::ptBinsUsed[nPtBinsUsed], 1.20);
  frame2->GetYaxis()->SetTitle("Unfolded/Truth");
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
  l2->AddEntry(gRatioCorr, "JES-corrected ratio");
  l2->Draw();

  float paErr = (paErrLow+paErrHigh)/2.0;
  TLatex jestext;
  jestext.SetNDC();
  jestext.SetTextColor(kRed);
  jestext.DrawLatex(.18,.28, Form("Data to MC JES (unfolded vs. truth%s) = %.4f #pm %.4f", shapeMethod ? ", shape #chi^{2}" : "", pa, paErr));

  c->SaveAs(pdfPath);
}

// xJ shape for one pT bin: truth, unfolded Data at pa=1 and at the best-fit pa (unit area).
void drawXjPage(TCanvas * c, const char * pdfPath, const char * label, int ir, float ptlow, float pthigh,
    TH1D * hTruth, TH1D * hUnfoldRaw, TH1D * hUnfoldCorr) {
  c->Clear();
  c->cd();
  gPad->SetLeftMargin(.15);
  gPad->SetTicks();

  auto shapeNorm = [](TH1D * h, const char * name) {
    TH1D * hd = (TH1D*)h->Clone(name);
    hd->Scale(1., "width");
    if (hd->Integral() > 0) hd->Scale(1./hd->Integral());
    return hd;
  };
  TH1D * hTruthdisp = shapeNorm(hTruth, Form("%s_disp_truth", hUnfoldCorr->GetName()));
  TH1D * hRawdisp   = shapeNorm(hUnfoldRaw, Form("%s_disp_raw", hUnfoldCorr->GetName()));
  TH1D * hCorrdisp  = shapeNorm(hUnfoldCorr, Form("%s_disp", hUnfoldCorr->GetName()));

  hTruthdisp->SetLineColor(kBlack);
  hTruthdisp->SetMarkerColor(kBlack);
  hTruthdisp->SetLineWidth(2);
  hTruthdisp->GetXaxis()->SetTitle("x_{J#gamma}");
  hTruthdisp->GetYaxis()->SetTitle("Shape-normalized counts / bin width");
  hTruthdisp->GetYaxis()->SetRangeUser(0, std::max({hTruthdisp->GetMaximum(), hRawdisp->GetMaximum(), hCorrdisp->GetMaximum()})*1.4);
  hTruthdisp->Draw("hist");

  hRawdisp->SetLineColor(kGray+2);
  hRawdisp->SetMarkerColor(kGray+2);
  hRawdisp->SetMarkerStyle(24); // open circle
  hRawdisp->SetLineWidth(2);
  hRawdisp->Draw("p e same");

  hCorrdisp->SetLineColor(kRed);
  hCorrdisp->SetMarkerColor(kRed);
  hCorrdisp->SetMarkerStyle(21);
  hCorrdisp->SetLineWidth(2);
  hCorrdisp->Draw("p e same");

  TLegend * l = new TLegend(.55,.5,.85,.7);
  l->SetLineWidth(0);
  l->AddEntry(hTruthdisp, "Pythia8 #gamma+jet (truth)");
  l->AddEntry(hRawdisp,   "Data (unfolded, raw JES)");
  l->AddEntry(hCorrdisp,  "Data (unfolded, JES-corrected)");
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
void grid_insitu_unfolded(string systag = "nominal", int na = insitu_utility::scanN, string method = "mean") {
  if (method != "mean" && method != "shape") { cout << "method must be \"mean\" or \"shape\"" << endl; return; }
  shapeMethod = (method == "shape");
  const char * tag = shapeMethod ? "unfolded_shapechi2" : "unfolded";
  // No auto-registration: the per-scan-point temporaries are deleted explicitly.
  TH1::AddDirectory(kFALSE);

  // Response source and output shared across radii.
  drawer d("pythia", systag);

  string pdfPathStr = Form("%s/grid_insitu_%s_%s.pdf", insitu_pdf_dir, tag, systag.c_str());
  TCanvas * c = new TCanvas("c","",700,700);
  c->SaveAs(Form("%s[", pdfPathStr.c_str()));

  string outfilename = Form("%s/grid_insitu_%s_%s.root", insitu_output_dir, tag, systag.c_str());
  TFile * fout = TFile::Open(outfilename.c_str(), "RECREATE");

  for (int ir = 0; ir < ana::nJetR; ir++) {

  string dataFile = insitu_utility::insituFilename(insitu_input_dir, "Data", "", systag);
  vector<DataEvent> dataA = insitu_utility::cacheDataEvents(dataFile.c_str(), 0, ir, false);
  vector<DataEvent> dataC = insitu_utility::cacheDataEvents(dataFile.c_str(), 2, ir, false);
  cout << "Cached Data events: region A=" << dataA.size() << " region C=" << dataC.size() << endl;

  // Low-xJ floor for every ana::ptBins bin.
  float lowXj[ana::nPtBins];
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) lowXj[ipt] = insitu_utility::lowXjFloor(ir, ana::ptBins[ipt]);

  // Purity for every ana::ptBins bin, fixed across the scan.
  float purity[ana::nPtBins], purityC[ana::nPtBins];
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    purity[ipt]  = ana::getPurity(ana::ptBins[ipt], ana::ptBins[ipt+1], systag, ir);
    purityC[ipt] = ana::getPurityC(ana::ptBins[ipt], ana::ptBins[ipt+1], systag, ir);
    cout << "Purity pt bin " << ipt << " [" << ana::ptBins[ipt] << "," << ana::ptBins[ipt+1]
         << "): P_A=" << purity[ipt] << " P_C=" << purityC[ipt] << endl;
  }

  // Response: Photon5/10/20 combined, nominal MC JES, fixed across the scan.
  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i", ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i", ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i", ir), 1);
  RooUnfoldResponse * response = new RooUnfoldResponse(respRecoTemplate, respTruthTemplate, respMatrix2D);

  // Truth mean xJ per used pT bin (not rescaled).
  float truthMean[nPtBinsUsed], truthMeanErr[nPtBinsUsed];
  vector<vector<double>> truthFrac(nPtBinsUsed), truthFracErr(nPtBinsUsed);
  for (int k = 0; k < nPtBinsUsed; k++) {
    int ipt = ana::firstUsedPtBin + k;
    TH1D * hTruth = unfold_utility::unflattenXj(respTruthTemplate, ipt, "hTruthMean_tmp");
    truthMean[k] = hTruth->GetMean();
    truthMeanErr[k] = hTruth->GetMeanError();
    double N = hTruth->Integral();
    truthFrac[k].assign(ana::nUnfoldXjBins, 0.);
    truthFracErr[k].assign(ana::nUnfoldXjBins, 0.);
    for (int ixj = 1; N > 0 && ixj <= ana::nUnfoldXjBins; ixj++) {
      truthFrac[k][ixj-1]    = hTruth->GetBinContent(ixj)/N;
      truthFracErr[k][ixj-1] = hTruth->GetBinError(ixj)/N;
    }
    cout << "Truth <x_J> pt bin " << k << " [" << ana::ptBinsUsed[k] << "," << ana::ptBinsUsed[k+1]
         << "): " << truthMean[k] << " +/- " << truthMeanErr[k] << endl;
    delete hTruth;
  }

  // -----------------------------
  // Grid scan over pa: purity-correct, unfold, chi2 against truth
  // -----------------------------
  const float lowa = insitu_utility::scanLow, higha = insitu_utility::scanHigh;

  TGraph * gchisqUnfold = new TGraph(na);
  gchisqUnfold->SetName("gchisq_unfolded");
  gchisqUnfold->SetTitle(shapeMethod ? ";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});Shape #chi^{2}" : ";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});#chi^{2}");

  float minchisqUnfold = FLT_MAX, minpaUnfold = 1;
  int ibestUnfold = 0;

  for (int ia = 0; ia < na; ia++) {
    float pa = lowa + ia*(higha-lowa)/na;

    float chisq = 0;
    if (shapeMethod) {
      chisq = computeUnfoldedShapeChi2(dataA, dataC, pa, purity, purityC, response, respRecoTemplate, truthFrac, truthFracErr, lowXj);
    } else {
    float mean[nPtBinsUsed], err[nPtBinsUsed];
    computeUnfoldedMeans(dataA, dataC, pa, purity, purityC, response, respRecoTemplate, mean, err, lowXj);
    for (int k = 0; k < nPtBinsUsed; k++) {
      if (truthMean[k] <= 0) continue;
      double diff = 1 - mean[k]/truthMean[k];
      double errt = sqrt((err[k]*err[k])/(truthMean[k]*truthMean[k])
          + mean[k]*mean[k]*truthMeanErr[k]*truthMeanErr[k]/pow(truthMean[k],4));
      if (errt > 0) chisq += diff*diff/(errt*errt);
    }
    }

    gchisqUnfold->SetPoint(ia, pa, chisq);
    if (chisq < minchisqUnfold) { minchisqUnfold = chisq; minpaUnfold = pa; ibestUnfold = ia; }

    if (ia % 100 == 0) cout << "  scan " << ia << "/" << na << ": pa=" << pa << " chi2=" << chisq << endl;
  }

  float errLowUnfold, errHighUnfold;
  insitu_utility::findError(gchisqUnfold, ibestUnfold, minchisqUnfold, errLowUnfold, errHighUnfold);

  cout << "\nFINAL RESULT (jet R=" << ana::JetRs[ir] << ", systag=" << systag << ", unfolded vs. truth" << (shapeMethod ? ", shape chi2" : "") << ")\n";
  cout << "Purity-corrected + unfolded:  p_a = " << minpaUnfold
       << " +" << errHighUnfold << "/-" << errLowUnfold << " (chi2=" << minchisqUnfold << ")" << endl;

  // -----------------------------
  // Final comparisons at pa=1, the best-fit pa, and truth
  // -----------------------------
  gStyle->SetOptStat(0);

  float rawMean[nPtBinsUsed], rawErr[nPtBinsUsed];
  float bestMean[nPtBinsUsed], bestErr[nPtBinsUsed];
  vector<TH1D*> hUnfoldRaw(nPtBinsUsed), hUnfoldBest(nPtBinsUsed);
  computeUnfoldedMeans(dataA, dataC, 1.0,         purity, purityC, response, respRecoTemplate, rawMean,  rawErr,  lowXj, &hUnfoldRaw);
  computeUnfoldedMeans(dataA, dataC, minpaUnfold, purity, purityC, response, respRecoTemplate, bestMean, bestErr, lowXj, &hUnfoldBest);

  vector<TH1D*> hTruthPt(nPtBinsUsed);
  for (int k = 0; k < nPtBinsUsed; k++) {
    int ipt = ana::firstUsedPtBin + k;
    hTruthPt[k] = unfold_utility::unflattenXj(respTruthTemplate, ipt, Form("hxjtruth_pt%d", k));
  }

  TGraphErrors * gTruth      = insitu_utility::meanGraph(truthMean, truthMeanErr, "gMeanTruth");
  TGraphErrors * gUnfoldRaw  = insitu_utility::meanGraph(rawMean,  rawErr,  "gMeanUnfold_raw");
  TGraphErrors * gRatioRaw   = insitu_utility::ratioGraph(rawMean,  rawErr,  truthMean, truthMeanErr, "gRatio_unfold_raw");
  TGraphErrors * gRatioCorr  = insitu_utility::ratioGraph(bestMean, bestErr, truthMean, truthMeanErr, "gRatio_unfold_corrected");

  const char * pdfPath = pdfPathStr.c_str();
  drawJESPage(c, pdfPath, "Purity-corrected, unfolded", ir, gTruth, gUnfoldRaw, gRatioRaw, gRatioCorr,
      minpaUnfold, errLowUnfold, errHighUnfold);
  for (int k = 0; k < nPtBinsUsed; k++) {
    drawXjPage(c, pdfPath, "Purity-corrected, unfolded", ir, ana::ptBinsUsed[k], ana::ptBinsUsed[k+1],
        hTruthPt[k], hUnfoldRaw[k], hUnfoldBest[k]);
  }

  // -----------------------------
  // Save
  // -----------------------------
  fout->cd();
  fout->mkdir(ana::rnames[ir])->cd();
  gchisqUnfold->Write();
  gTruth->Write();
  gUnfoldRaw->Write();
  gRatioRaw->Write();
  gRatioCorr->Write();
  for (int k = 0; k < nPtBinsUsed; k++) {
    hTruthPt[k]->Write();
    hUnfoldRaw[k]->Write();
    hUnfoldBest[k]->Write();
  }

  TTree * wt = new TTree("results", shapeMethod ? "best-fit jet energy scale result (unfolded vs. truth, shape chi2)" : "best-fit jet energy scale result (unfolded vs. truth)");
  float wpa = minpaUnfold, wchisq = minchisqUnfold, werrLow = errLowUnfold, werrHigh = errHighUnfold;
  wt->Branch("pa_unfolded", &wpa);
  wt->Branch("chisq_unfolded", &wchisq);
  wt->Branch("errLow_unfolded", &werrLow);
  wt->Branch("errHigh_unfolded", &werrHigh);
  wt->Fill();
  wt->Write();

  cout << "Finished ir=" << ir << " (" << ana::rnames[ir] << ")" << endl;
  } // end of ir loop

  c->SaveAs(Form("%s]", pdfPathStr.c_str()));
  cout << "Wrote " << pdfPathStr << endl;
  fout->Close();
  cout << "Wrote " << outfilename << endl;
}
