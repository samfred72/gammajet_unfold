#include "../src/ana.h"
#include "../src/insitu_utility.h"
#include "../src/unfold_utility.h"
#include "../src/purity_utility.h"
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <cfloat>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TBox.h"
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

// In-situ JES scan on the insitu trees.
//   mode "gammajet" (default): one constant scale pa (jet_pt/pa) chosen so Data's mean xJ (or xJ
//     shape) matches the Pythia8 reference per photon-pT bin, for region A alone and for the
//     purity-corrected A/C combination. draw_jes_summary.C turns this into ana.h's constants.
//   mode "combined": gamma+jet plus multijet/analysis.cc's balance, fitting f(pT) = pa + pb*pT
//     on a 2D grid (multijet fixes the slope, gamma+jet the normalization). Cross-check only.
// Regions A and C are scaled by the same pa. The Data purity depends on pa through the pairing (jet pT
// cut, x_J floor), so in gammajet mode it is recomputed from the Data insitu tree at the scan result and
// the scan repeated until the two agree (self-consistent; no re-unfolding needed). fixedPurity, and
// combined mode, use purity_<systag>.root, made at the pa of the last unfolding.

// insitu/: inputs/ (insitu trees), output/ (.root), pdfs/.
const char * insitu_input_dir  = ana::path("insitu/inputs");
const char * insitu_output_dir = ana::path("insitu/output");
const char * insitu_pdf_dir    = ana::path("insitu/pdfs");

// Reported pT bins only (ana::ptBinsUsed).
const int nPtBinsUsed = ana::nPtBinsUsed;
// Shape method: xJ bin fractions; last 3 sparse xJ bins dropped, pT bin 2 merged in pairs.
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;
const int coarseRebinPtBin = 2;
const int coarseGroupSize = 2;
bool shapeMethod = false; // set by grid_insitu(..., method)

// Photon5/10/20 cross-section weights (drawer.h's scalemap).
map<int,double> photon_scale = {{5,146359.3},{10,6944.675},{20,130.4461}};
map<int,double> photon_scale_herwig = {{5,6.48487e+07},{10,3.62808e+02},{20,5.34010e+01}};

// Combined mode: directory of multijet/analysis.cc's output (selection and weights are its).
string multijet_analysis_dir = ana::path("multijet"); // non-const so a test can point it elsewhere
// Combined-mode grid: pa coarser than the 1D scan, pb over +-0.005/GeV. 1 sigma: delta-chi2 < 2.30.
const int combinedNa = 200;
const int combinedNb = 100;
const float combinedLowb = -0.005, combinedHighb = 0.005; // a first test fit sat at +0.002
const float combinedDchi2 = 2.30;

// Mean xJ vs pT for MC and raw Data (top); raw and corrected Data/MC ratios (bottom). The
// corrected ratio is flat at 1 by construction.
void drawJESPage(TCanvas * c, const char * pdfPath, const char * label, int ir,
    TGraphErrors * gMC, TGraphErrors * gDataRaw, TGraphErrors * gRatioRaw, TGraphErrors * gRatioCorr,
    float pa, float paErrLow, float paErrHigh, const char * extraText = "") {
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
  l1->AddEntry(gMC, "Pythia8 #gamma+jet (reco)");
  l1->AddEntry(gDataRaw, "Data (reco)");
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
  jestext.DrawLatex(.18,.28, Form("Data to MC JES%s = %.4f #pm %.4f", shapeMethod ? " (shape #chi^{2})" : "", pa, paErr));
  if (extraText[0]) jestext.DrawLatex(.18,.22, extraText);

  c->SaveAs(pdfPath);
}

// xJ for one pT bin: MC reference, raw Data, and Data at the best-fit pa; shape-normalized
// densities. Purity-corrected histograms can go negative; normalize only when the integral
// is positive.
void drawXjPage(TCanvas * c, const char * pdfPath, const char * label, int ir, float ptlow, float pthigh,
    TH1D * hMC, TH1D * hDataRaw, TH1D * hDataCorr, const char * dataRawLabel, const char * dataCorrLabel) {
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
  l->AddEntry(hMCdisp,   "Pythia8 #gamma+jet (reco)");
  l->AddEntry(hRawdisp,  dataRawLabel);
  l->AddEntry(hCorrdisp, dataCorrLabel);
  l->Draw();

  insitu_utility::drawSPhenixLabel({label, Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ptlow, pthigh)}, {
      Form("Jet R=%.1f", ana::JetRs[ir]),
      Form("p_{T}^{jet} > %.0f GeV", ana::jet_calib_pt_cut[ir]),
      Form("|#eta^{#gamma}|<%.1f, |#eta^{jet}|<%.1f", ana::photonEtaCut, ana::etacut-ana::JetRs[ir]),
      Form("#Delta#phi>%.0f#pi/%.0f", ana::oppnum, ana::oppden)
    }, .18, .85, 16, gPad->GetWh());

  c->SaveAs(pdfPath);
}

void runCombined(TCanvas * c, const char * pdfPath, TFile * fout, int ir, const string & systag,
    const vector<DataEvent> & dataA, const vector<DataEvent> & dataC,
    const float refMean[], const float refMeanErr[], const float purity[], const float purityC[],
    const float lowXj[]);

// onlyIr >= 0: one radius only. method = "mean" or "shape" (gammajet mode only). fixedPurity: hold the
// purity at purity_<systag>.root instead of the self-consistent iteration. paStart > 0: start that
// iteration there instead of at the ana.h table value.
void grid_insitu(string systag = "nominal", string mode = "gammajet", int onlyIr = -1, string method = "mean",
    bool fixedPurity = false, float paStart = -1) {
  if (mode != "gammajet" && mode != "combined") {
    cout << "ERROR: mode must be \"gammajet\" or \"combined\", got \"" << mode << "\"" << endl;
    return;
  }
  const bool combined = (mode == "combined");
  if (method != "mean" && method != "shape") { cout << "ERROR: method must be \"mean\" or \"shape\"" << endl; return; }
  shapeMethod = (method == "shape");
  if (shapeMethod && combined) { cout << "ERROR: the shape method is gammajet-mode only" << endl; return; }
  const string tag = combined ? "combined_" + systag : (shapeMethod ? "shapechi2_" + systag : systag);
  // systag "herwig" (generator systematic, not an ana::systags entry): the MC reference is Herwig;
  // Data, purity and everything else are nominal.
  const bool herwig = (systag == "herwig");
  if (herwig && combined) { cout << "ERROR: herwig is gammajet-mode only" << endl; return; }
  const string inSystag = herwig ? "nominal" : systag; // input trees, purity, starting pa
  const char * mcSim = herwig ? "herwig" : "pythia";
  map<int,double> & mcScale = herwig ? photon_scale_herwig : photon_scale;
  // Don't auto-register new histograms into whichever radius directory is current.
  TH1::AddDirectory(kFALSE);

  // One file and one PDF per systag; each radius in its ana::rnames[ir] subdirectory.
  string pdfPathStr = Form("%s/grid_insitu_%s.pdf", insitu_pdf_dir, tag.c_str());
  TCanvas * c = new TCanvas("c","",700,700);
  c->SaveAs(Form("%s[", pdfPathStr.c_str()));

  string outfilename = Form("%s/grid_insitu_%s.root", insitu_output_dir, tag.c_str());
  TFile * fout = TFile::Open(outfilename.c_str(), "RECREATE");

  for (int ir = 0; ir < ana::nJetR; ir++) {
  if (onlyIr >= 0 && ir != onlyIr) continue;

  string dataFile = insitu_utility::insituFilename(insitu_input_dir, "Data", "", inSystag);
  vector<DataEvent> dataA = insitu_utility::cacheDataEvents(dataFile.c_str(), 0, ir);
  vector<DataEvent> dataC = insitu_utility::cacheDataEvents(dataFile.c_str(), 2, ir);
  cout << "Cached Data events: region A=" << dataA.size() << " region C=" << dataC.size() << endl;

  // Low-xJ floor per pT bin, as unfolder::check_pair applies (insitu_utility::lowXjFloor).
  float lowXj[nPtBinsUsed];
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) lowXj[ipt] = insitu_utility::lowXjFloor(ir, ana::ptBinsUsed[ipt]);

  float refMean[nPtBinsUsed], refMeanErr[nPtBinsUsed];
  vector<pair<string,double>> mcSamples = {
    {insitu_utility::insituFilename(insitu_input_dir, "Photon5",  mcSim, inSystag), mcScale[5]},
    {insitu_utility::insituFilename(insitu_input_dir, "Photon10", mcSim, inSystag), mcScale[10]},
    {insitu_utility::insituFilename(insitu_input_dir, "Photon20", mcSim, inSystag), mcScale[20]},
  };
  insitu_utility::referenceMeans(mcSamples, 0, ir, refMean, refMeanErr, lowXj);
  vector<vector<double>> refFrac, refFracErr;
  if (shapeMethod) insitu_utility::referenceShape(mcSamples, 0, ir, refFrac, refFracErr, lowXj);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    cout << "MC reference <x_J> pt bin " << ipt << " [" << ana::ptBinsUsed[ipt] << "," << ana::ptBinsUsed[ipt+1]
         << "): " << refMean[ipt] << " +/- " << refMeanErr[ipt] << endl;
  }

  // Purity per pT bin, fixed within one scan. Self-consistent (gammajet mode unless fixedPurity): Data
  // counts from the insitu tree at paPurity (purity_utility::dataCountsFromTree, the unfolder's pairing)
  // and the MC leakage fractions, solved as puritymaker.C. threejet at R=0.2 uses the nominal tree and
  // leakage, as ana::getPurity does (ana.cc).
  float purity[nPtBinsUsed], purityErrLow[nPtBinsUsed], purityErrHigh[nPtBinsUsed];
  float purityC[nPtBinsUsed], purityCErrLow[nPtBinsUsed], purityCErrHigh[nPtBinsUsed];
  const bool selfConsistent = !combined && !fixedPurity;
  const string puritySystag = (systag == "threejet" && ir == 0) ? "nominal" : inSystag;
  const string purityDataFile = insitu_utility::insituFilename(insitu_input_dir, "Data", "", puritySystag);
  TH1D * leak[4];
  TH1D * counts[4];
  if (selfConsistent) {
    purity_utility::leakageFractions(puritySystag, ir, leak);
    for (int j = 0; j < 4; j++) {
      counts[j] = new TH1D(Form("purity_counts_%d", j), "", ana::nPtBins, ana::ptBins);
      counts[j]->Sumw2();
    }
  }
  auto setPurity = [&](float paPur) {
    if (!selfConsistent) {
      for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
        purity[ipt]        = ana::getPurity(ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1], inSystag, ir);
        purityErrLow[ipt]  = ana::getPurityErrorLow(ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1], inSystag, ir);
        purityErrHigh[ipt] = ana::getPurityErrorHigh(ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1], inSystag, ir);
        purityC[ipt]        = ana::getPurityC(ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1], inSystag, ir);
        purityCErrLow[ipt]  = ana::getPurityCErrorLow(ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1], inSystag, ir);
        purityCErrHigh[ipt] = ana::getPurityCErrorHigh(ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1], inSystag, ir);
        cout << "Purity pt bin " << ipt << ": P_A=" << purity[ipt] << " P_C=" << purityC[ipt] << endl;
      }
      return;
    }
    purity_utility::dataCountsFromTree(purityDataFile, ir, paPur, counts);
    TGraphAsymmErrors * gC = nullptr;
    TGraphAsymmErrors * gA = purity_utility::combine(counts, leak, &gC, false);
    for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
      int i = ana::findPtBin((ana::ptBinsUsed[ipt] + ana::ptBinsUsed[ipt+1])/2.0); // as ana::getPurity
      purity[ipt]  = gA->GetPointY(i); purityErrLow[ipt]  = gA->GetErrorYlow(i); purityErrHigh[ipt]  = gA->GetErrorYhigh(i);
      purityC[ipt] = gC->GetPointY(i); purityCErrLow[ipt] = gC->GetErrorYlow(i); purityCErrHigh[ipt] = gC->GetErrorYhigh(i);
      cout << "Purity pt bin " << ipt << " at p_a=" << paPur << ": P_A=" << purity[ipt] << " P_C=" << purityC[ipt] << endl;
    }
    delete gA; delete gC;
  };
  float paPurity = (paStart > 0) ? paStart : ana::jesForSystag(inSystag, ir);
  setPurity(paPurity);

  if (combined) {
    runCombined(c, pdfPathStr.c_str(), fout, ir, systag, dataA, dataC, refMean, refMeanErr, purity, purityC, lowXj);
    continue;
  }

  // -----------------------------
  // Grid scan over pa (window in insitu_utility.h)
  // -----------------------------
  const int na = insitu_utility::scanN;
  const float lowa = insitu_utility::scanLow, higha = insitu_utility::scanHigh;

  // Self-consistency loop: scan at fixed purity; if the purity-corrected result is not the grid point the
  // purity was made at, remake the purity there and rescan. The counts are discrete, so purity(p_a) can
  // make the result flip between neighbouring grid points: on a cycle (a result already scanned), the
  // scan in the cycle with the lowest chi2 minimum is redone and its result reported.
  const int maxSelfConsistentIter = 10;
  vector<float> scannedPa, scannedChi2; // purity point of each scan and its purity-corrected chi2 minimum
  bool cycleRescan = false;
  TGraph * gchisqA = nullptr;
  TGraph * gchisqCorr = nullptr;
  float minchisqA, minpaA, minchisqCorr, minpaCorr;
  int ibestA, ibestCorr;
  for (int iterSC = 1; ; iterSC++) {
  delete gchisqA;
  delete gchisqCorr;
  gchisqA    = new TGraph(na);
  gchisqCorr = new TGraph(na);
  gchisqA->SetName("gchisq_regionA");
  gchisqA->SetTitle(shapeMethod ? ";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});Shape #chi^{2}" : ";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});#chi^{2}");
  gchisqCorr->SetName("gchisq_puritycorrected");
  gchisqCorr->SetTitle(shapeMethod ? ";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});Shape #chi^{2}" : ";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});#chi^{2}");

  minchisqA = FLT_MAX; minpaA = 1;
  minchisqCorr = FLT_MAX; minpaCorr = 1;
  ibestA = 0; ibestCorr = 0;

  for (int ia = 0; ia < na; ia++) {
    float pa = lowa + ia*(higha-lowa)/na;

    float chisqA = 0, chisqCorr = 0;
    if (shapeMethod) {
      vector<vector<double>> countA(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
      vector<vector<double>> countC(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
      for (auto & ev : dataA) {
        if (insitu_utility::vetoed(ev, pa)) continue;
        float x = (ev.jet_pt/pa)/ev.pho_pt;
        if (x < lowXj[ev.ptbin]) continue;
        int ixj = ana::findUnfoldXjBin(x);
        if (ixj < 0 || ixj >= ana::nUnfoldXjBins) continue;
        countA[ev.ptbin][ixj] += 1;
      }
      for (auto & ev : dataC) {
        if (insitu_utility::vetoed(ev, pa)) continue;
        float x = (ev.jet_pt/pa)/ev.pho_pt;
        if (x < lowXj[ev.ptbin]) continue;
        int ixj = ana::findUnfoldXjBin(x);
        if (ixj < 0 || ixj >= ana::nUnfoldXjBins) continue;
        countC[ev.ptbin][ixj] += 1;
      }

      for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
        double NA = 0;
        for (double n : countA[ipt]) NA += n;

        // pT bin 2 uses the coarsened bins.
        vector<double> useA, useC, useRefFrac, useRefFracErr;
        int nBins;
        if (ipt == coarseRebinPtBin) {
          int startBin = ana::findUnfoldXjBin(lowXj[ipt]);
          useA         = insitu_utility::coarsenSum(countA[ipt], startBin, nXjBinsForChi2, coarseGroupSize);
          useC         = insitu_utility::coarsenSum(countC[ipt], startBin, nXjBinsForChi2, coarseGroupSize);
          useRefFrac   = insitu_utility::coarsenSum(refFrac[ipt], startBin, nXjBinsForChi2, coarseGroupSize);
          useRefFracErr = insitu_utility::coarsenQuadrature(refFracErr[ipt], startBin, nXjBinsForChi2, coarseGroupSize);
          nBins = (int)useA.size();
        } else {
          useA = countA[ipt]; useC = countC[ipt];
          useRefFrac = refFrac[ipt]; useRefFracErr = refFracErr[ipt];
          nBins = nXjBinsForChi2;
        }

        // Region A: pull^2 per xJ bin between Data's and MC's bin fractions. errt is floored at 1/NA
        // (a zero-count bin otherwise has zero error and gives huge pulls).
        if (NA > 0) {
          double errFloor = 1.0/NA;
          for (int ib = 0; ib < nBins; ib++) {
            double fData = useA[ib]/NA;
            double errData = sqrt(useA[ib])/NA;
            double errt = sqrt(errData*errData + useRefFracErr[ib]*useRefFracErr[ib]);
            errt = std::max(errt, errFloor);
            double diff = fData - useRefFrac[ib];
            chisqA += diff*diff/(errt*errt);
          }
        }

        // Purity-corrected: the same linear coeffA/coeffC on per-bin counts; error from A/C counting
        // only (purity fixed).
        double NC = 0;
        for (double n : countC[ipt]) NC += n;
        if (NA > 0 && NC > 0) {
          float coeffA, coeffC;
          unfold_utility::purityCorrectCoeffs(purity[ipt], purityC[ipt], NA, NC, coeffA, coeffC);
          double Ncorr = coeffA*NA - coeffC*NC;
          if (Ncorr > 0) {
            // Same 1/Ncorr floor.
            double errFloor = 1.0/Ncorr;
            for (int ib = 0; ib < nBins; ib++) {
              double corrCount = coeffA*useA[ib] - coeffC*useC[ib];
              double fCorr = corrCount/Ncorr;
              double errCorrCount = sqrt(coeffA*coeffA*useA[ib] + coeffC*coeffC*useC[ib]);
              double errCorr = errCorrCount/Ncorr;
              double errt = sqrt(errCorr*errCorr + useRefFracErr[ib]*useRefFracErr[ib]);
              errt = std::max(errt, errFloor);
              double diff = fCorr - useRefFrac[ib];
              chisqCorr += diff*diff/(errt*errt);
            }
          }
        }
      }
    } else {
      vector<double> sumA(nPtBinsUsed,0), sumA2(nPtBinsUsed,0);
      vector<int> countA(nPtBinsUsed,0);
      for (auto & ev : dataA) {
        if (insitu_utility::vetoed(ev, pa)) continue;
        float x = (ev.jet_pt/pa)/ev.pho_pt;
        if (x < lowXj[ev.ptbin]) continue;
        sumA[ev.ptbin]  += x;
        sumA2[ev.ptbin] += x*x;
        countA[ev.ptbin]++;
      }
      vector<double> sumC(nPtBinsUsed,0), sumC2(nPtBinsUsed,0);
      vector<int> countC(nPtBinsUsed,0);
      for (auto & ev : dataC) {
        if (insitu_utility::vetoed(ev, pa)) continue;
        float x = (ev.jet_pt/pa)/ev.pho_pt;
        if (x < lowXj[ev.ptbin]) continue;
        sumC[ev.ptbin]  += x;
        sumC2[ev.ptbin] += x*x;
        countC[ev.ptbin]++;
      }

      for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
        if (refMean[ipt] <= 0) continue;

        // Region A only.
        if (countA[ipt] > 0) {
          double mean = sumA[ipt]/countA[ipt];
          double var  = sumA2[ipt]/countA[ipt] - mean*mean;
          double err  = sqrt(std::max(var,0.)/countA[ipt]);
          double diff = 1 - mean/refMean[ipt];
          double errt = sqrt((err*err)/(refMean[ipt]*refMean[ipt])
              + mean*mean*refMeanErr[ipt]*refMeanErr[ipt]/pow(refMean[ipt],4));
          if (errt > 0) chisqA += diff*diff/(errt*errt);
        }

        // Purity-corrected.
        if (countA[ipt] > 0 && countC[ipt] > 0) {
          double NA = countA[ipt], NC = countC[ipt];
          float coeffA, coeffC;
          unfold_utility::purityCorrectCoeffs(purity[ipt], purityC[ipt], NA, NC, coeffA, coeffC);
          double sumXcorr  = coeffA*sumA[ipt]  - coeffC*sumC[ipt];
          double sumX2corr = coeffA*sumA2[ipt] - coeffC*sumC2[ipt];
          double Ncorr = coeffA*NA - coeffC*NC;
          if (Ncorr > 0) {
            double mean = sumXcorr/Ncorr;
            double var  = sumX2corr/Ncorr - mean*mean;
            double err  = sqrt(std::max(var,0.)/Ncorr);
            double diff = 1 - mean/refMean[ipt];
            double errt = sqrt((err*err)/(refMean[ipt]*refMean[ipt])
                + mean*mean*refMeanErr[ipt]*refMeanErr[ipt]/pow(refMean[ipt],4));
            if (errt > 0) chisqCorr += diff*diff/(errt*errt);
          }
        }
      }
    }

    gchisqA->SetPoint(ia, pa, chisqA);
    gchisqCorr->SetPoint(ia, pa, chisqCorr);

    if (chisqA < minchisqA)       { minchisqA = chisqA;       minpaA = pa;       ibestA = ia; }
    if (chisqCorr < minchisqCorr) { minchisqCorr = chisqCorr; minpaCorr = pa; ibestCorr = ia; }
  }

  if (!selfConsistent) break;
  cout << "Self-consistent purity, iteration " << iterSC << ": purity at p_a=" << paPurity
       << " -> purity-corrected p_a=" << minpaCorr << endl;
  if (cycleRescan) break;
  const float halfStep = 0.5f*(higha-lowa)/na;
  if (fabs(minpaCorr - paPurity) < halfStep) break;  // same grid point
  scannedPa.push_back(paPurity);
  scannedChi2.push_back(minchisqCorr);
  int icycle = -1;
  for (int i = 0; i < (int)scannedPa.size(); i++) if (fabs(scannedPa[i] - minpaCorr) < halfStep) icycle = i;
  if (icycle >= 0) {
    int ibest = icycle;
    for (int i = icycle; i < (int)scannedPa.size(); i++) if (scannedChi2[i] < scannedChi2[ibest]) ibest = i;
    cout << "Self-consistent purity: cycle over " << scannedPa.size() - icycle << " scans; using the scan with purity at p_a="
         << scannedPa[ibest] << " (lowest chi2, " << scannedChi2[ibest] << ")" << endl;
    paPurity = scannedPa[ibest];
    setPurity(paPurity);
    cycleRescan = true;
    continue;
  }
  if (iterSC == maxSelfConsistentIter) {
    cout << "WARNING: self-consistent purity not converged after " << maxSelfConsistentIter
         << " iterations (last p_a " << paPurity << " -> " << minpaCorr << "); using the last scan." << endl;
    break;
  }
  paPurity = minpaCorr;
  setPurity(paPurity);
  }

  float errLowA, errHighA, errLowCorr, errHighCorr;
  insitu_utility::findError(gchisqA,    ibestA,    minchisqA,    errLowA,    errHighA);
  insitu_utility::findError(gchisqCorr, ibestCorr, minchisqCorr, errLowCorr, errHighCorr);

  cout << "\nFINAL RESULT (jet R=" << ana::JetRs[ir] << ", systag=" << systag << (shapeMethod ? ", shape chi2" : "") << ")\n";
  cout << "Region A only:        p_a = " << minpaA
       << " +" << errHighA << "/-" << errLowA << " (chi2=" << minchisqA << ")" << endl;
  cout << "Purity-corrected:     p_a = " << minpaCorr
       << " +" << errHighCorr << "/-" << errLowCorr << " (chi2=" << minchisqCorr << ")" << endl;

  // -----------------------------
  // xJ histograms per pT bin: MC reference, raw Data, Data at each best-fit pa
  // -----------------------------
  vector<TH1D*> hxjMC_pt       = insitu_utility::buildMCXjByPtBin(mcSamples, 0, ir, "hxjA_pythia", lowXj);
  vector<TH1D*> hxjA_raw_pt    = insitu_utility::buildXjByPtBin(dataA, 1.0,       nPtBinsUsed, "hxjA_data_raw", lowXj);
  vector<TH1D*> hxjA_corr_pt   = insitu_utility::buildXjByPtBin(dataA, minpaA,    nPtBinsUsed, "hxjA_data_corr", lowXj);
  vector<TH1D*> hxjC_raw_pt    = insitu_utility::buildXjByPtBin(dataC, 1.0,       nPtBinsUsed, "hxjC_data_raw", lowXj);
  vector<TH1D*> hxjA_atCorr_pt = insitu_utility::buildXjByPtBin(dataA, minpaCorr, nPtBinsUsed, "hxjA_data_atCorrScale", lowXj);
  vector<TH1D*> hxjC_atCorr_pt = insitu_utility::buildXjByPtBin(dataC, minpaCorr, nPtBinsUsed, "hxjC_data_atCorrScale", lowXj);
  // Purity-corrected, at pa=1 and at the best-fit pa.
  vector<TH1D*> hxjcorr_raw_pt  = insitu_utility::purityCorrectByPtBin(hxjA_raw_pt,    hxjC_raw_pt, nPtBinsUsed,
      purity, purityErrLow, purityErrHigh, purityC, purityCErrLow, purityCErrHigh, "hxjcorrected_data_raw");
  vector<TH1D*> hxjcorr_best_pt = insitu_utility::purityCorrectByPtBin(hxjA_atCorr_pt, hxjC_atCorr_pt, nPtBinsUsed,
      purity, purityErrLow, purityErrHigh, purityC, purityCErrLow, purityCErrHigh, "hxjcorrected_data_bestscale");

  auto sumPtBins = [&](const vector<TH1D*> & h, const char * name) {
    TH1D * hsum = (TH1D*)h[0]->Clone(name);
    for (int ipt = 1; ipt < nPtBinsUsed; ipt++) hsum->Add(h[ipt]);
    return hsum;
  };
  TH1D * hxjA_pythia                 = sumPtBins(hxjMC_pt,        "hxjA_pythia");
  TH1D * hxjA_data_raw               = sumPtBins(hxjA_raw_pt,     "hxjA_data_raw");
  TH1D * hxjA_data_bestscale         = sumPtBins(hxjA_corr_pt,    "hxjA_data_bestscale");
  TH1D * hxjcorrected_data_raw       = sumPtBins(hxjcorr_raw_pt,  "hxjcorrected_data_raw");
  TH1D * hxjcorrected_data_bestscale = sumPtBins(hxjcorr_best_pt, "hxjcorrected_data_bestscale");

  // -----------------------------
  // Mean xJ vs pT pages (region A and purity-corrected)
  // -----------------------------
  gStyle->SetOptStat(0);

  float rawMeanA[nPtBinsUsed], rawErrA[nPtBinsUsed];
  float corrMeanA[nPtBinsUsed], corrErrA[nPtBinsUsed]; // "corr" here = best-fit-scaled, not purity-corrected
  insitu_utility::computeRegionAMeans(dataA, 1.0,    rawMeanA,  rawErrA,  lowXj);
  insitu_utility::computeRegionAMeans(dataA, minpaA, corrMeanA, corrErrA, lowXj);

  float rawMeanCorr[nPtBinsUsed], rawErrCorr[nPtBinsUsed];
  float bestMeanCorr[nPtBinsUsed], bestErrCorr[nPtBinsUsed];
  insitu_utility::computeCorrectedMeans(dataA, dataC, 1.0,        purity, purityC, rawMeanCorr,  rawErrCorr,  lowXj);
  insitu_utility::computeCorrectedMeans(dataA, dataC, minpaCorr,  purity, purityC, bestMeanCorr, bestErrCorr, lowXj);

  TGraphErrors * gMC              = insitu_utility::meanGraph(refMean, refMeanErr, "gMeanMC");
  TGraphErrors * gDataRaw_A       = insitu_utility::meanGraph(rawMeanA,  rawErrA,  "gMeanData_regionA_raw");
  TGraphErrors * gRatioRaw_A      = insitu_utility::ratioGraph(rawMeanA,  rawErrA,  refMean, refMeanErr, "gRatio_regionA_raw");
  TGraphErrors * gRatioCorr_A     = insitu_utility::ratioGraph(corrMeanA, corrErrA, refMean, refMeanErr, "gRatio_regionA_corrected");

  TGraphErrors * gDataRaw_Corr    = insitu_utility::meanGraph(rawMeanCorr,  rawErrCorr,  "gMeanData_puritycorrected_raw");
  TGraphErrors * gRatioRaw_Corr   = insitu_utility::ratioGraph(rawMeanCorr,  rawErrCorr,  refMean, refMeanErr, "gRatio_puritycorrected_raw");
  TGraphErrors * gRatioCorr_Corr  = insitu_utility::ratioGraph(bestMeanCorr, bestErrCorr, refMean, refMeanErr, "gRatio_puritycorrected_corrected");

  const char * pdfPath = pdfPathStr.c_str();
  drawJESPage(c, pdfPath, "Region A", ir, gMC, gDataRaw_A, gRatioRaw_A, gRatioCorr_A, minpaA, errLowA, errHighA);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    drawXjPage(c, pdfPath, "Region A", ir, ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1],
        hxjMC_pt[ipt], hxjA_raw_pt[ipt], hxjA_corr_pt[ipt], "Data (raw)", "Data (JES-corrected)");
  }
  drawJESPage(c, pdfPath, "Purity-corrected", ir, gMC, gDataRaw_Corr, gRatioRaw_Corr, gRatioCorr_Corr, minpaCorr, errLowCorr, errHighCorr);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    drawXjPage(c, pdfPath, "Purity-corrected", ir, ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1],
        hxjMC_pt[ipt], hxjcorr_raw_pt[ipt], hxjcorr_best_pt[ipt], "Data (raw, purity-corr.)", "Data (JES+purity-corr.)");
  }

  // -----------------------------
  // Save (mkdir/cd before the TTree is constructed)
  // -----------------------------
  fout->cd();
  fout->mkdir(ana::rnames[ir])->cd();
  gchisqA->Write();
  gchisqCorr->Write();
  hxjA_pythia->Write();
  hxjA_data_raw->Write();
  hxjA_data_bestscale->Write();
  hxjcorrected_data_raw->Write();
  hxjcorrected_data_bestscale->Write();

  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    hxjMC_pt[ipt]->Write();
    hxjA_raw_pt[ipt]->Write();
    hxjA_corr_pt[ipt]->Write();
    hxjcorr_raw_pt[ipt]->Write();
    hxjcorr_best_pt[ipt]->Write();
  }

  gMC->Write();
  gDataRaw_A->Write();
  gRatioRaw_A->Write();
  gRatioCorr_A->Write();
  gDataRaw_Corr->Write();
  gRatioRaw_Corr->Write();
  gRatioCorr_Corr->Write();

  TTree * wt = new TTree("results", shapeMethod ? "best-fit jet energy scale results (shape chi2)" : "best-fit jet energy scale results");
  float wpaA = minpaA, wchisqA = minchisqA, werrLowA = errLowA, werrHighA = errHighA;
  float wpaCorr = minpaCorr, wchisqCorr = minchisqCorr, werrLowCorr = errLowCorr, werrHighCorr = errHighCorr;
  wt->Branch("pa_regionA", &wpaA);
  wt->Branch("chisq_regionA", &wchisqA);
  wt->Branch("errLow_regionA", &werrLowA);
  wt->Branch("errHigh_regionA", &werrHighA);
  wt->Branch("pa_puritycorrected", &wpaCorr);
  wt->Branch("chisq_puritycorrected", &wchisqCorr);
  wt->Branch("errLow_puritycorrected", &werrLowCorr);
  wt->Branch("errHigh_puritycorrected", &werrHighCorr);
  wt->Fill();
  wt->Write();

  cout << "Finished ir=" << ir << " (" << ana::rnames[ir] << ")" << endl;
  } // end of ir loop

  c->SaveAs(Form("%s]", pdfPathStr.c_str()));
  cout << "Wrote " << pdfPathStr << endl;
  fout->Close();
  cout << "Wrote " << outfilename << endl;
}

// =============================================================================
// Combined mode: gamma+jet + multijet balance, linear JES f(pT) = pa + pb*pT.
// =============================================================================

// Gamma+jet sums of x = (jet_pt/f(jet_pt))/pho_pt per pT bin, with the low-xJ floor.
static void gammaSumsLinear(const vector<DataEvent> & ev, double pa, double pb, const float lowXj[],
    vector<double> & sum, vector<double> & sum2, vector<int> & count) {
  sum.assign(nPtBinsUsed, 0); sum2.assign(nPtBinsUsed, 0); count.assign(nPtBinsUsed, 0);
  for (auto & e : ev) {
    if (insitu_utility::vetoed(e, pa + pb*e.third_pt)) continue;
    double x = (e.jet_pt/(pa + pb*e.jet_pt))/e.pho_pt;
    if (x < lowXj[e.ptbin]) continue;
    sum[e.ptbin] += x; sum2[e.ptbin] += x*x; count[e.ptbin]++;
  }
}

// Region-A and purity-corrected mean xJ per pT bin from those sums.
static void gammaMeansLinear(const vector<DataEvent> & dataA, const vector<DataEvent> & dataC,
    double pa, double pb, const float purity[], const float purityC[], const float lowXj[],
    float meanA[], float errA[], float meanCorr[], float errCorr[]) {
  vector<double> sA, sA2, sC, sC2; vector<int> nA, nC;
  gammaSumsLinear(dataA, pa, pb, lowXj, sA, sA2, nA);
  gammaSumsLinear(dataC, pa, pb, lowXj, sC, sC2, nC);
  for (int i = 0; i < nPtBinsUsed; i++) {
    meanA[i] = errA[i] = meanCorr[i] = errCorr[i] = 0;
    if (nA[i] > 0) {
      double m = sA[i]/nA[i];
      meanA[i] = m;
      errA[i] = sqrt(std::max(sA2[i]/nA[i] - m*m, 0.)/nA[i]);
    }
    if (nA[i] > 0 && nC[i] > 0) {
      float cA, cC;
      unfold_utility::purityCorrectCoeffs(purity[i], purityC[i], nA[i], nC[i], cA, cC);
      double N = cA*nA[i] - cC*nC[i];
      if (N > 0) {
        double m = (cA*sA[i] - cC*sC[i])/N;
        meanCorr[i] = m;
        errCorr[i] = sqrt(std::max((cA*sA2[i] - cC*sC2[i])/N - m*m, 0.)/N);
      }
    }
  }
}

// chi2 of measured means against a reference.
static double meanChi2(const float mean[], const float err[], const float ref[], const float refErr[], int n) {
  double chi2 = 0;
  for (int i = 0; i < n; i++) {
    if (ref[i] <= 0 || mean[i] <= 0) continue;
    double diff = 1 - mean[i]/ref[i];
    double errt = sqrt(err[i]*err[i]/(ref[i]*ref[i]) + mean[i]*mean[i]*refErr[i]*refErr[i]/pow(ref[i],4));
    if (errt > 0) chi2 += diff*diff/(errt*errt);
  }
  return chi2;
}

// Mean multijet balance vs leading-jet pT (top); raw and corrected Data/MC (bottom).
static void drawMultijetPage(TCanvas * c, const char * pdfPath, int ir,
    const float mcMean[], const float mcErr[], const float rawMean[], const float rawErr[],
    const float corrMean[], const float corrErr[], float pa, float pb) {
  const int n = insitu_utility::nMultijetPtBins;
  const double * b = insitu_utility::multijetPtBins;
  auto graph = [&](const float y[], const float ey[], const char * name) {
    TGraphErrors * g = new TGraphErrors();
    g->SetName(name);
    for (int i = 0; i < n; i++) {
      if (y[i] <= 0) continue;
      int k = g->GetN();
      g->SetPoint(k, 0.5*(b[i]+b[i+1]), y[i]);
      g->SetPointError(k, 0.5*(b[i+1]-b[i]), ey[i]);
    }
    return g;
  };
  auto ratio = [&](const float y[], const float ey[], const char * name) {
    float r[insitu_utility::nMultijetPtBins], er[insitu_utility::nMultijetPtBins];
    for (int i = 0; i < n; i++) {
      r[i] = (mcMean[i] > 0 && y[i] > 0) ? y[i]/mcMean[i] : 0;
      er[i] = r[i] > 0 ? r[i]*sqrt(pow(ey[i]/y[i],2) + pow(mcErr[i]/mcMean[i],2)) : 0;
    }
    return graph(r, er, name);
  };
  TGraphErrors * gMC = graph(mcMean, mcErr, "gMultijetMC");
  TGraphErrors * gRaw = graph(rawMean, rawErr, "gMultijetDataRaw");
  TGraphErrors * rRaw = ratio(rawMean, rawErr, "gMultijetRatioRaw");
  TGraphErrors * rCorr = ratio(corrMean, corrErr, "gMultijetRatioCorr");

  c->Clear();
  TPad * p1 = new TPad("p1","",0,.5,1,1);
  TPad * p2 = new TPad("p2","",0,0,1,.5);
  p1->Draw(); p2->Draw();
  p1->cd(); p1->SetBottomMargin(0.02); p1->SetLeftMargin(.15); gPad->SetTicks(1,1);
  TH1F * f1 = p1->DrawFrame(b[0], 0.8, b[n], 1.6);
  f1->GetYaxis()->SetTitle("<p_{T}^{lead}/|#vec{p}_{T}^{sub}+#vec{p}_{T}^{subsub}|>");
  f1->GetXaxis()->SetLabelSize(0);
  gMC->SetLineColor(kMagenta+1); gMC->SetMarkerColor(kMagenta+1); gMC->SetMarkerStyle(21); gMC->Draw("p same");
  gRaw->SetLineColor(kBlue); gRaw->SetMarkerColor(kBlue); gRaw->SetMarkerStyle(20); gRaw->Draw("p same");
  TLegend * l1 = new TLegend(.55,.1,.85,.3); l1->SetLineWidth(0);
  l1->AddEntry(gMC, "Pythia8 dijet (reco)"); l1->AddEntry(gRaw, "Data (reco)"); l1->Draw();
  insitu_utility::drawSPhenixLabel({"Multijet balance"}, {
      Form("Jet R=%.1f", ana::JetRs[ir]),
      "selection and MC weights: multijet/analysis.cc"
    }, .18, .85, 16, p1->GetWh()/1.5);
  p2->cd(); p2->SetTopMargin(0.02); p2->SetBottomMargin(0.2); p2->SetLeftMargin(.15); gPad->SetTicks(1,1);
  TH1F * f2 = p2->DrawFrame(b[0], 0.85, b[n], 1.15);
  f2->GetYaxis()->SetTitle("Data/MC");
  f2->GetXaxis()->SetTitle("p_{T}^{lead} [GeV]");
  f2->GetYaxis()->SetTitleSize(0.06); f2->GetYaxis()->SetTitleOffset(1.1); f2->GetYaxis()->SetLabelSize(0.05);
  f2->GetXaxis()->SetTitleSize(0.06); f2->GetXaxis()->SetLabelSize(0.05);
  rRaw->SetMarkerStyle(20); rRaw->Draw("p same");
  rCorr->SetMarkerStyle(24); rCorr->Draw("p same");
  TLine * line = new TLine(b[0],1,b[n],1); line->SetLineStyle(9); line->Draw("same");
  TLegend * l2 = new TLegend(.55,.7,.85,.9); l2->SetLineWidth(0);
  l2->AddEntry(rRaw, "Raw ratio"); l2->AddEntry(rCorr, "Corrected ratio"); l2->Draw();
  TLatex t; t.SetNDC(); t.SetTextColor(kRed);
  t.DrawLatex(.18,.28, Form("f(p_{T}) = %.4f %+.2e p_{T}", pa, pb));
  c->SaveAs(pdfPath);
  gMC->Write(); gRaw->Write(); rRaw->Write(); rCorr->Write();
}

void runCombined(TCanvas * c, const char * pdfPath, TFile * fout, int ir, const string & systag,
    const vector<DataEvent> & dataA, const vector<DataEvent> & dataC,
    const float refMean[], const float refMeanErr[], const float purity[], const float purityC[],
    const float lowXj[]) {
  const int nMJ = insitu_utility::nMultijetPtBins;

  // Multijet inputs: Data weight 1, MC weight from analysis.cc. JERhigh/low use its HIGH/LOW smear.
  const int radius = (int)std::lround(ana::JetRs[ir]*10);
  const string sys = insitu_utility::multijetSysName(systag);
  string mjFile = insitu_utility::multijetAnalysisFilename(multijet_analysis_dir.c_str(), "pythia");
  vector<MultijetEvent> mjData = insitu_utility::cacheMultijetEvents(mjFile.c_str(), Form("ttree_data_r%d", radius));
  vector<MultijetEvent> mjMC;
  for (int is = 0; is < insitu_utility::nMultijetMCSamples; is++) {
    const char * sample = insitu_utility::multijetMCSamples[is];
    vector<MultijetEvent> ev = insitu_utility::cacheMultijetEvents(mjFile.c_str(), Form("ttree_%s_r%d_%s", sample, radius, sys.c_str()));
    cout << "  multijet MC " << sample << " (" << sys << "): " << ev.size() << " events" << endl;
    mjMC.insert(mjMC.end(), ev.begin(), ev.end());
  }
  cout << "Cached multijet events: Data=" << mjData.size() << " MC=" << mjMC.size() << endl;
  if (mjData.empty() || mjMC.empty()) {
    cout << "ERROR: no multijet Data or MC events for R=" << ana::JetRs[ir] << " in " << mjFile
         << " - run multijet/analysis first (multijet/README.md); skipping this radius" << endl;
    return;
  }
  float mjRef[nMJ], mjRefErr[nMJ];
  insitu_utility::multijetMeans(mjMC, 1.0, 0.0, mjRef, mjRefErr, false); // MC trees: already selected

  // 2D grid scan.
  const float lowa = insitu_utility::scanLow, higha = insitu_utility::scanHigh;
  const float stepa = (higha-lowa)/combinedNa, stepb = (combinedHighb-combinedLowb)/combinedNb;
  TH2D * hA   = new TH2D("hchisq2d_regionA_multijet", ";p_{a};p_{b} [GeV^{-1}];#chi^{2}",
      combinedNa, lowa-0.5*stepa, higha-0.5*stepa, combinedNb, combinedLowb-0.5*stepb, combinedHighb-0.5*stepb);
  TH2D * hCorr = (TH2D*)hA->Clone("hchisq2d_puritycorrected_multijet");
  TH2D * hMJ   = (TH2D*)hA->Clone("hchisq2d_multijet_only");
  struct Pt { float pa, pb, cA, cCorr, cMJ; };
  vector<Pt> pts; pts.reserve(combinedNa*combinedNb);
  float gA[ana::nPtBinsUsed], eA[ana::nPtBinsUsed], gC[ana::nPtBinsUsed], eC[ana::nPtBinsUsed];
  float mj[nMJ], emj[nMJ];
  for (int ia = 0; ia < combinedNa; ia++) {
    float pa = lowa + ia*stepa;
    for (int ib = 0; ib < combinedNb; ib++) {
      float pb = combinedLowb + ib*stepb;
      gammaMeansLinear(dataA, dataC, pa, pb, purity, purityC, lowXj, gA, eA, gC, eC);
      insitu_utility::multijetMeans(mjData, pa, pb, mj, emj);
      float cMJ = meanChi2(mj, emj, mjRef, mjRefErr, nMJ);
      float cA = meanChi2(gA, eA, refMean, refMeanErr, nPtBinsUsed) + cMJ;
      float cCorr = meanChi2(gC, eC, refMean, refMeanErr, nPtBinsUsed) + cMJ;
      pts.push_back({pa, pb, cA, cCorr, cMJ});
      hA->SetBinContent(ia+1, ib+1, cA);
      hCorr->SetBinContent(ia+1, ib+1, cCorr);
      hMJ->SetBinContent(ia+1, ib+1, cMJ);
    }
  }

  // Best point and delta-chi2 < 2.30 region: pa/pb ranges and the f(pT) envelope.
  struct Fit { float pa, pb, chi2, paLo, paHi, pbLo, pbHi; TGraph * band[3]; bool edge; };
  auto fit = [&](float Pt::*chi, const char * name) {
    Fit r; r.chi2 = FLT_MAX;
    for (auto & p : pts) if (p.*chi < r.chi2) { r.chi2 = p.*chi; r.pa = p.pa; r.pb = p.pb; }
    r.paLo = r.paHi = r.pa; r.pbLo = r.pbHi = r.pb;
    const int nx = 56; // f(pT) band from 5 to 60 GeV
    vector<double> x(nx), lo(nx, 1e9), hi(nx, -1e9), best(nx);
    for (int i = 0; i < nx; i++) { x[i] = 5 + i; best[i] = r.pa + r.pb*x[i]; }
    for (auto & p : pts) {
      if (p.*chi >= r.chi2 + combinedDchi2) continue;
      r.paLo = std::min(r.paLo, p.pa); r.paHi = std::max(r.paHi, p.pa);
      r.pbLo = std::min(r.pbLo, p.pb); r.pbHi = std::max(r.pbHi, p.pb);
      for (int i = 0; i < nx; i++) { double f = p.pa + p.pb*x[i]; lo[i] = std::min(lo[i], f); hi[i] = std::max(hi[i], f); }
    }
    r.edge = (r.paLo <= lowa || r.paHi >= higha - stepa || r.pbLo <= combinedLowb || r.pbHi >= combinedHighb - stepb);
    r.band[0] = new TGraph(nx, x.data(), best.data()); r.band[0]->SetName(Form("gJES_best_%s", name));
    r.band[1] = new TGraph(nx, x.data(), lo.data());   r.band[1]->SetName(Form("gJES_low_%s", name));
    r.band[2] = new TGraph(nx, x.data(), hi.data());   r.band[2]->SetName(Form("gJES_high_%s", name));
    return r;
  };
  Fit fitA = fit(&Pt::cA, "regionA_multijet");
  Fit fitC = fit(&Pt::cCorr, "puritycorrected_multijet");

  cout << "\nCOMBINED RESULT (jet R=" << ana::JetRs[ir] << ", systag=" << systag << "), f(pT) = pa + pb*pT\n";
  for (auto p : {make_pair("Region A + multijet:        ", &fitA), make_pair("Purity-corrected + multijet:", &fitC)}) {
    const Fit & r = *p.second;
    cout << p.first << " pa = " << r.pa << " [" << r.paLo << ", " << r.paHi << "]"
         << "  pb = " << r.pb << " [" << r.pbLo << ", " << r.pbHi << "] /GeV  chi2 = " << r.chi2 << endl;
    if (r.edge) cout << "  WARNING: the delta-chi2 < 2.30 region touches the grid edge - widen the scan" << endl;
  }

  // Pages: gamma+jet at the best fits, multijet balance, f(pT) band, chi2 map.
  float rawA[ana::nPtBinsUsed], rawEA[ana::nPtBinsUsed], rawC[ana::nPtBinsUsed], rawEC[ana::nPtBinsUsed];
  float bA[ana::nPtBinsUsed], bEA[ana::nPtBinsUsed], bC[ana::nPtBinsUsed], bEC[ana::nPtBinsUsed], dum[ana::nPtBinsUsed], dumE[ana::nPtBinsUsed];
  gammaMeansLinear(dataA, dataC, 1.0, 0.0, purity, purityC, lowXj, rawA, rawEA, rawC, rawEC);
  gammaMeansLinear(dataA, dataC, fitA.pa, fitA.pb, purity, purityC, lowXj, bA, bEA, dum, dumE);
  gammaMeansLinear(dataA, dataC, fitC.pa, fitC.pb, purity, purityC, lowXj, dum, dumE, bC, bEC);

  fout->cd();
  fout->mkdir(ana::rnames[ir])->cd();

  TGraphErrors * gMC = insitu_utility::meanGraph(refMean, refMeanErr, "gMeanMC");
  drawJESPage(c, pdfPath, "Region A + multijet", ir, gMC,
      insitu_utility::meanGraph(rawA, rawEA, "gMeanData_regionA_raw"),
      insitu_utility::ratioGraph(rawA, rawEA, refMean, refMeanErr, "gRatio_regionA_raw"),
      insitu_utility::ratioGraph(bA, bEA, refMean, refMeanErr, "gRatio_regionA_multijet_corrected"),
      fitA.pa, fitA.pa - fitA.paLo, fitA.paHi - fitA.pa, Form("p_{b} = %.2e [%.2e, %.2e] GeV^{-1}", fitA.pb, fitA.pbLo, fitA.pbHi));
  drawJESPage(c, pdfPath, "Purity-corrected + multijet", ir, gMC,
      insitu_utility::meanGraph(rawC, rawEC, "gMeanData_puritycorrected_raw"),
      insitu_utility::ratioGraph(rawC, rawEC, refMean, refMeanErr, "gRatio_puritycorrected_raw"),
      insitu_utility::ratioGraph(bC, bEC, refMean, refMeanErr, "gRatio_puritycorrected_multijet_corrected"),
      fitC.pa, fitC.pa - fitC.paLo, fitC.paHi - fitC.pa, Form("p_{b} = %.2e [%.2e, %.2e] GeV^{-1}", fitC.pb, fitC.pbLo, fitC.pbHi));

  float mjRaw[nMJ], mjRawE[nMJ], mjCorr[nMJ], mjCorrE[nMJ];
  insitu_utility::multijetMeans(mjData, 1.0, 0.0, mjRaw, mjRawE);
  insitu_utility::multijetMeans(mjData, fitC.pa, fitC.pb, mjCorr, mjCorrE);
  drawMultijetPage(c, pdfPath, ir, mjRef, mjRefErr, mjRaw, mjRawE, mjCorr, mjCorrE, fitC.pa, fitC.pb);

  // f(pT) bands, with the gammajet-only scale for comparison if present.
  c->Clear(); c->cd(); gPad->SetLeftMargin(.15); gPad->SetBottomMargin(.12); gPad->SetTicks(1,1);
  TH1F * fr = gPad->DrawFrame(5, 0.75, 60, 1.05);
  fr->GetXaxis()->SetTitle("p_{T}^{jet} [GeV]");
  fr->GetYaxis()->SetTitle("f(p_{T}) = p_{a} + p_{b} p_{T}  (jet_{pt,corrected} = jet_{pt}/f)");
  int cols[2] = {kBlue, kRed};
  Fit * fits[2] = {&fitC, &fitA};
  const char * names[2] = {"Purity-corrected + multijet", "Region A + multijet"};
  TLegend * lb = new TLegend(.45,.15,.9,.35); lb->SetLineWidth(0);
  for (int k = 0; k < 2; k++) {
    for (int j = 0; j < 3; j++) { fits[k]->band[j]->SetLineColor(cols[k]); fits[k]->band[j]->SetLineWidth(j ? 1 : 2); fits[k]->band[j]->SetLineStyle(j ? 2 : 1); fits[k]->band[j]->Draw("l same"); fits[k]->band[j]->Write(); }
    lb->AddEntry(fits[k]->band[0], names[k], "l");
  }
  TFile * fg = TFile::Open(Form("%s/grid_insitu_%s.root", insitu_output_dir, systag.c_str()), "READ");
  if (fg && !fg->IsZombie()) {
    TTree * rt = (TTree*)fg->Get(Form("%s/results", ana::rnames[ir]));
    float paG = 0, eLo = 0, eHi = 0;
    if (rt) {
      rt->SetBranchAddress("pa_puritycorrected", &paG); rt->SetBranchAddress("errLow_puritycorrected", &eLo);
      rt->SetBranchAddress("errHigh_puritycorrected", &eHi); rt->GetEntry(0);
      TBox * bx = new TBox(5, paG - eLo, 60, paG + eHi); bx->SetFillColor(kGray); bx->SetFillStyle(1001); bx->Draw();
      TLine * lg = new TLine(5, paG, 60, paG); lg->SetLineColor(kBlack); lg->SetLineWidth(2); lg->Draw();
      lb->AddEntry(lg, "#gamma+jet only (constant, purity-corr.)", "l");
      for (int k = 0; k < 2; k++) for (int j = 0; j < 3; j++) fits[k]->band[j]->Draw("l same"); // redraw over the box
    }
    fg->Close();
  }
  fout->cd(ana::rnames[ir]);
  lb->Draw();
  insitu_utility::drawSPhenixLabel({"#gamma+jet + multijet"}, {Form("Jet R=%.1f", ana::JetRs[ir]), "solid: best fit, dashed: #Delta#chi^{2} < 2.30"}, .18, .85, 16, gPad->GetWh());
  c->SaveAs(pdfPath);

  // chi2 map (purity-corrected + multijet)
  c->Clear(); c->cd(); gPad->SetRightMargin(.15); gPad->SetLeftMargin(.15); gPad->SetBottomMargin(.12);
  TH2D * hDisp = (TH2D*)hCorr->Clone("hchisq2d_display");
  hDisp->SetMaximum(fitC.chi2 + 25); hDisp->SetMinimum(fitC.chi2);
  hDisp->Draw("colz");
  TGraph * gb = new TGraph(1); gb->SetPoint(0, fitC.pa, fitC.pb); gb->SetMarkerStyle(29); gb->SetMarkerSize(2); gb->SetMarkerColor(kRed); gb->Draw("p same");
  insitu_utility::drawSPhenixLabel({"#chi^{2}: purity-corrected #gamma+jet + multijet"}, {Form("Jet R=%.1f", ana::JetRs[ir])}, .2, .85, 14, gPad->GetWh());
  c->SaveAs(pdfPath);
  gPad->SetRightMargin(.05);

  fout->cd(ana::rnames[ir]);
  hA->Write(); hCorr->Write(); hMJ->Write();
  TTree * wt = new TTree("results", "best-fit linear JES f(pT) = pa + pb*pT");
  float w[2][7];
  for (int k = 0; k < 2; k++) { Fit & r = (k ? fitA : fitC); float v[7] = {r.pa, r.pb, r.chi2, r.paLo, r.paHi, r.pbLo, r.pbHi}; for (int j = 0; j < 7; j++) w[k][j] = v[j]; }
  const char * lab[7] = {"pa", "pb", "chisq", "paLow", "paHigh", "pbLow", "pbHigh"};
  for (int j = 0; j < 7; j++) {
    wt->Branch(Form("%s_puritycorrected_multijet", lab[j]), &w[0][j]);
    wt->Branch(Form("%s_regionA_multijet", lab[j]), &w[1][j]);
  }
  int nData = mjData.size(), nMC = mjMC.size();
  wt->Branch("nMultijetData", &nData); wt->Branch("nMultijetMC", &nMC);
  wt->Fill(); wt->Write();
  cout << "Finished combined ir=" << ir << " (" << ana::rnames[ir] << ")" << endl;
}
