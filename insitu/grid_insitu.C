#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/insitu_utility.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/unfold_utility.h"
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

// ana::findPtBin/getPurity/etc. live in ana.cc, compiled into libgammajet_unfold.so -
// load it explicitly (see draw_purity_corrected.C / draw_insitu_xj.C) so cling resolves
// the real compiled definitions instead of misbinding against the sibling gammajet
// project's own ana/drawer classes on the same library path.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// In-situ jet-energy-scale study, in the style of gammajet/macros/run_grid.sh's
// grid_insitu.C, but reading this project's insitutree files (see draw_insitu_xj.C)
// instead of the sibling gammajet project's tree_Data.root, and gammajet-only (no
// multijet/trijet reference - we don't have those trees here).
//
// Scans a single overall jet-energy-scale factor pa (jet_pt_corrected = jet_pt/pa,
// no pT-dependence) and finds the value that makes Data's mean(x_{J#gamma}) match the
// fixed Pythia8 gamma+jet MC reference mean, per photon-pT bin, two ways:
//   1. Region A alone (the naive fit - background-contaminated).
//   2. The purity-corrected combination of Region A and Region C, since the physical
//      photon signal is what the JES should actually be tuned against.
// Comparing the two best-fit scales is the actual "insitu correction study": it shows
// how much the region-C background pulls a naive region-A-only JES fit away from the
// purity-corrected answer.
//
// Region A and Region C jets are scaled by the *same* trial pa at every grid point
// ("scale all the jets in regions A and C, then do the correction"). The purity P is
// computed once, outside the pa loop, and held fixed across the whole scan - jet energy
// scale doesn't move photon isolation/BDT, so it can't move which events fall in A vs.
// C, or the measured purity fraction itself ("the purity is independent of jet energy
// scale, so we can use the same purity the whole way through").

// insitu/ is split into inputs/ (the raw Data/Photon insitu ntuples, written by
// unfolder.h's production pipeline), output/ (this and the other grid_insitu*.C
// macros' own .root output), and pdfs/ (their .pdf output).
const char * insitu_input_dir  = "/home/samson72/sphnx/gammajet_unfold/insitu/inputs";
const char * insitu_output_dir = "/home/samson72/sphnx/gammajet_unfold/insitu/output";
const char * insitu_pdf_dir    = "/home/samson72/sphnx/gammajet_unfold/insitu/pdfs";

// Only ana::ptBinsUsed (15-20, 20-25, 25-35 GeV) is used for every calculation and
// plot below - both the low-pT migration-only buffer bin (13-15 GeV, ana::ptBins[0])
// and the high-pT overflow bin (35-100 GeV) are excluded, since neither is a reported
// physics bin (see ana.h's ptBins/ptBinsUsed/firstUsedPtBin comment) and the top one
// also has too few Data events for a meaningful in-situ point.
const int nPtBinsUsed = ana::nPtBinsUsed;

// Cross-section weights for combining the Photon5/10/20 MC samples - same numbers as
// drawer.h's scalemap[isphoton=1][sample] for sim="pythia".
map<int,double> photon_scale = {{5,146359.3},{10,6944.675},{20,130.4461}};

// referenceMeans, computeRegionAMeans, computeCorrectedMeans, buildXjByPtBin,
// buildMCXjByPtBin, and purityCorrectByPtBin now live in src/insitu_utility.h/.cc
// (insitu_utility:: namespace) - moved there after being found copy-pasted
// byte-for-byte across all six grid_insitu*.C macros (see that header's comment).

// One comparison page: top panel is mean(x_J) vs pT for MC and raw Data; bottom panel
// is the raw ratio (raw Data/MC) and the corrected ratio (best-fit-scaled Data/MC) -
// the latter should sit flat at 1 by construction, since pa was fit to make it so.
// pa/paErrLow/paErrHigh are the best-fit in-situ jet-energy-scale factor for this page
// (Region A or purity-corrected), stamped on the bottom panel in red.
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
  l1->AddEntry(gMC, "Pythia8 #gamma+jet (reco)");
  l1->AddEntry(gDataRaw, "Data (reco)");
  l1->Draw();
  insitu_utility::drawSPhenixLabel({label}, {
      Form("Jet R=%.1f", ana::JetRs[ir]),
      Form("|#eta^{#gamma}|<%.1f, |#eta^{jet}|<%.1f", ana::etacut, ana::etacut-ana::JetRs[ir]),
      Form("#Delta#phi>%.0f#pi/%.0f", ana::oppnum, ana::oppden)
    }, .18, .85, 16, p1->GetWh()/1.5);

  p2->cd();
  p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.2);
  p2->SetLeftMargin(.15);
  gPad->SetTicks(1,1);
  TH1F * frame2 = p2->DrawFrame(ana::ptBinsUsed[0], 0.90, ana::ptBinsUsed[nPtBinsUsed], 1.10);
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
  jestext.DrawLatex(.18,.28, Form("Data to MC JES = %.4f #pm %.4f", pa, paErr));

  c->SaveAs(pdfPath);
}

// One xJ-distribution comparison page, for a single photon-pT bin: the fixed MC
// reference, the raw (uncorrected) Data distribution, and the Data distribution at the
// study's best-fit jet-energy-scale (and, for the purity-corrected study, also
// background-subtracted) - all shape-normalized and divided by bin width for display,
// since ana::unfoldXjBins is non-uniform (same densityForDisplay convention as
// draw_insitu_xj.C/draw_purity_corrected.C). The purity-corrected histograms can go
// bin-by-bin negative (region C oversubtracting a noisy bin) - Integral() is only used
// to normalize when positive.
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
      Form("|#eta^{#gamma}|<%.1f, |#eta^{jet}|<%.1f", ana::etacut, ana::etacut-ana::JetRs[ir]),
      Form("#Delta#phi>%.0f#pi/%.0f", ana::oppnum, ana::oppden)
    }, .18, .85, 16, gPad->GetWh());

  c->SaveAs(pdfPath);
}

void grid_insitu(string systag = "nominal") {
  // Newly created histograms are not auto-registered to whatever TDirectory happens to
  // be gDirectory at construction time - without this, buildXjByPtBin/buildMCXjByPtBin's
  // fixed-name per-pT-bin histograms would auto-register into (and "Replacing existing
  // TH1" warn against) whichever radius subdirectory was left current by the *previous*
  // iteration's mkdir/cd below, since they're built before this iteration's own mkdir/cd
  // runs. Every actual save still goes through this file's explicit ->Write() calls
  // (unaffected by this setting), same precedent as grid_insitu_unfolded.C.
  TH1::AddDirectory(kFALSE);

  // One file/one PDF for the whole systag, all seven jet radii inside - opened/created
  // here, before the per-radius loop, instead of grid_insitu.C's old per-radius
  // filenames. Each radius's objects land in their own ana::rnames[ir] subdirectory of
  // fout (mkdir/cd'd right before that radius's own "Save" block below - a TTree binds
  // to whatever TDirectory is current at construction time, so this has to happen
  // before "results" is constructed, not just before its Write()); each radius's pages
  // become one more page in the same multi-page PDF via the standard "file.pdf["/
  // "file.pdf"/"file.pdf]" SaveAs bracket already used per-pT-bin below, just wrapped
  // one level higher.
  string pdfPathStr = Form("%s/grid_insitu_%s.pdf", insitu_pdf_dir, systag.c_str());
  TCanvas * c = new TCanvas("c","",700,700);
  c->SaveAs(Form("%s[", pdfPathStr.c_str()));

  string outfilename = Form("%s/grid_insitu_%s.root", insitu_output_dir, systag.c_str());
  TFile * fout = TFile::Open(outfilename.c_str(), "RECREATE");

  for (int ir = 0; ir < ana::nJetR; ir++) {

  string dataFile = insitu_utility::insituFilename(insitu_input_dir, "Data", "", systag);
  vector<DataEvent> dataA = insitu_utility::cacheDataEvents(dataFile.c_str(), 0, ir);
  vector<DataEvent> dataC = insitu_utility::cacheDataEvents(dataFile.c_str(), 2, ir);
  cout << "Cached Data events: region A=" << dataA.size() << " region C=" << dataC.size() << endl;

  // Low-xJ floor per used pT bin - same cut unfolder::check_pair applies at floorScale=1
  // before a reco jet enters hrecoxj/the response matrix (see src/insitu_utility.h's
  // lowXjFloor comment). Applied below to every mean(x_J)/shape computation (Data and
  // MC reference alike) so this scan excludes exactly the events the main unfolding
  // pipeline would exclude at the same jet radius.
  float lowXj[nPtBinsUsed];
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) lowXj[ipt] = insitu_utility::lowXjFloor(ir, ana::ptBinsUsed[ipt]);

  float refMean[nPtBinsUsed], refMeanErr[nPtBinsUsed];
  insitu_utility::referenceMeans({
      {insitu_utility::insituFilename(insitu_input_dir, "Photon5",  "pythia", systag), photon_scale[5]},
      {insitu_utility::insituFilename(insitu_input_dir, "Photon10", "pythia", systag), photon_scale[10]},
      {insitu_utility::insituFilename(insitu_input_dir, "Photon20", "pythia", systag), photon_scale[20]},
    }, 0, ir, refMean, refMeanErr, lowXj);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    cout << "MC reference <x_J> pt bin " << ipt << " [" << ana::ptBinsUsed[ipt] << "," << ana::ptBinsUsed[ipt+1]
         << "): " << refMean[ipt] << " +/- " << refMeanErr[ipt] << endl;
  }

  // Purity per photon-pT bin (region A and region C) - computed once, held fixed across
  // the whole pa scan. Error arrays are only needed by purityCorrectByPtBin's asymmetric
  // purity-uncertainty term below; the scalar-moment grid scan/computeCorrectedMeans use
  // only the central values (purity has always been held fixed, not scanned, here).
  float purity[nPtBinsUsed], purityErrLow[nPtBinsUsed], purityErrHigh[nPtBinsUsed];
  float purityC[nPtBinsUsed], purityCErrLow[nPtBinsUsed], purityCErrHigh[nPtBinsUsed];
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    purity[ipt]        = ana::getPurity(ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1], systag, ir);
    purityErrLow[ipt]  = ana::getPurityErrorLow(ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1], systag, ir);
    purityErrHigh[ipt] = ana::getPurityErrorHigh(ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1], systag, ir);
    purityC[ipt]        = ana::getPurityC(ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1], systag, ir);
    purityCErrLow[ipt]  = ana::getPurityCErrorLow(ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1], systag, ir);
    purityCErrHigh[ipt] = ana::getPurityCErrorHigh(ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1], systag, ir);
    cout << "Purity pt bin " << ipt << ": P_A=" << purity[ipt] << " P_C=" << purityC[ipt] << endl;
  }

  // -----------------------------
  // Grid scan: single overall jet-energy-scale factor pa, no pT-dependence - matches
  // run_grid.sh's gammajet-only mode (nb=1, pb=0). Scan window/step live in
  // insitu_utility.h (scanLow/scanHigh/scanN) so they're shared across every
  // grid_insitu*.C.
  // -----------------------------
  const int na = insitu_utility::scanN;
  const float lowa = insitu_utility::scanLow, higha = insitu_utility::scanHigh;

  TGraph * gchisqA    = new TGraph(na);
  TGraph * gchisqCorr = new TGraph(na);
  gchisqA->SetName("gchisq_regionA");
  gchisqA->SetTitle(";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});#chi^{2}");
  gchisqCorr->SetName("gchisq_puritycorrected");
  gchisqCorr->SetTitle(";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});#chi^{2}");

  float minchisqA = FLT_MAX, minpaA = 1;
  float minchisqCorr = FLT_MAX, minpaCorr = 1;
  int ibestA = 0, ibestCorr = 0;

  for (int ia = 0; ia < na; ia++) {
    float pa = lowa + ia*(higha-lowa)/na;

    vector<double> sumA(nPtBinsUsed,0), sumA2(nPtBinsUsed,0);
    vector<int> countA(nPtBinsUsed,0);
    for (auto & ev : dataA) {
      float x = (ev.jet_pt/pa)/ev.pho_pt;
      if (x < lowXj[ev.ptbin]) continue;
      sumA[ev.ptbin]  += x;
      sumA2[ev.ptbin] += x*x;
      countA[ev.ptbin]++;
    }
    vector<double> sumC(nPtBinsUsed,0), sumC2(nPtBinsUsed,0);
    vector<int> countC(nPtBinsUsed,0);
    for (auto & ev : dataC) {
      float x = (ev.jet_pt/pa)/ev.pho_pt;
      if (x < lowXj[ev.ptbin]) continue;
      sumC[ev.ptbin]  += x;
      sumC2[ev.ptbin] += x*x;
      countC[ev.ptbin]++;
    }

    float chisqA = 0, chisqCorr = 0;
    for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
      if (refMean[ipt] <= 0) continue;

      // Region A only (naive, background-contaminated fit).
      if (countA[ipt] > 0) {
        double mean = sumA[ipt]/countA[ipt];
        double var  = sumA2[ipt]/countA[ipt] - mean*mean;
        double err  = sqrt(std::max(var,0.)/countA[ipt]);
        double diff = 1 - mean/refMean[ipt];
        double errt = sqrt((err*err)/(refMean[ipt]*refMean[ipt])
            + mean*mean*refMeanErr[ipt]*refMeanErr[ipt]/pow(refMean[ipt],4));
        if (errt > 0) chisqA += diff*diff/(errt*errt);
      }

      // Purity-corrected (two-purity method, both regions scaled by the same pa) - see
      // unfold_utility::purityCorrectCoeffs for the coeffA/coeffC derivation.
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

    gchisqA->SetPoint(ia, pa, chisqA);
    gchisqCorr->SetPoint(ia, pa, chisqCorr);

    if (chisqA < minchisqA)       { minchisqA = chisqA;       minpaA = pa;       ibestA = ia; }
    if (chisqCorr < minchisqCorr) { minchisqCorr = chisqCorr; minpaCorr = pa; ibestCorr = ia; }
  }

  float errLowA, errHighA, errLowCorr, errHighCorr;
  insitu_utility::findError(gchisqA,    ibestA,    minchisqA,    errLowA,    errHighA);
  insitu_utility::findError(gchisqCorr, ibestCorr, minchisqCorr, errLowCorr, errHighCorr);

  cout << "\nFINAL RESULT (jet R=" << ana::JetRs[ir] << ", systag=" << systag << ")\n";
  cout << "Region A only:        p_a = " << minpaA
       << " +" << errHighA << "/-" << errLowA << " (chi2=" << minchisqA << ")" << endl;
  cout << "Purity-corrected:     p_a = " << minpaCorr
       << " +" << errHighCorr << "/-" << errLowCorr << " (chi2=" << minchisqCorr << ")" << endl;

  // -----------------------------
  // Build x_J histograms per photon-pT bin: the fixed MC reference, Data at pa=1 (raw)
  // and at each study's own best-fit pa (corrected) - same non-uniform binning as
  // draw_insitu_xj.C. Also used to build the inclusive (summed-over-pT-bin) histograms
  // saved to the output ROOT file below, for continuity with earlier versions of this
  // macro.
  // -----------------------------
  vector<pair<string,double>> mcSamples = {
    {insitu_utility::insituFilename(insitu_input_dir, "Photon5",  "pythia", systag), photon_scale[5]},
    {insitu_utility::insituFilename(insitu_input_dir, "Photon10", "pythia", systag), photon_scale[10]},
    {insitu_utility::insituFilename(insitu_input_dir, "Photon20", "pythia", systag), photon_scale[20]},
  };
  vector<TH1D*> hxjMC_pt       = insitu_utility::buildMCXjByPtBin(mcSamples, 0, ir, "hxjA_pythia", lowXj);
  vector<TH1D*> hxjA_raw_pt    = insitu_utility::buildXjByPtBin(dataA, 1.0,       nPtBinsUsed, "hxjA_data_raw", lowXj);
  vector<TH1D*> hxjA_corr_pt   = insitu_utility::buildXjByPtBin(dataA, minpaA,    nPtBinsUsed, "hxjA_data_corr", lowXj);
  vector<TH1D*> hxjC_raw_pt    = insitu_utility::buildXjByPtBin(dataC, 1.0,       nPtBinsUsed, "hxjC_data_raw", lowXj);
  vector<TH1D*> hxjA_atCorr_pt = insitu_utility::buildXjByPtBin(dataA, minpaCorr, nPtBinsUsed, "hxjA_data_atCorrScale", lowXj);
  vector<TH1D*> hxjC_atCorr_pt = insitu_utility::buildXjByPtBin(dataC, minpaCorr, nPtBinsUsed, "hxjC_data_atCorrScale", lowXj);
  // Purity-corrected (two-purity method), once raw (pa=1) and once at the
  // purity-corrected study's best-fit pa - purity is fixed either way (file header).
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
  // Mean(x_J) vs pT comparison plots - top: MC vs raw Data; bottom: raw ratio vs
  // corrected ratio (evaluated at each study's own best-fit pa). One page for Region A
  // alone, one page for the purity-corrected combination.
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
  // Save - into this radius's own subdirectory of the shared, once-opened fout (see
  // top of function). mkdir/cd has to happen before "results" (a TTree) is
  // constructed below, not just before its Write().
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

  TTree * wt = new TTree("results", "best-fit jet energy scale results");
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
