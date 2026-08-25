#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/insitu_utility.h"
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

// ana::findPtBin/etc. live in ana.cc, compiled into libgammajet_unfold.so - load it
// explicitly (see grid_insitu.C) so cling resolves the real compiled definitions.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Shape-chi2 variant of grid_insitu_jet12.C, exactly as grid_insitu_shapechi2.C is the
// shape-chi2 variant of grid_insitu.C: instead of picking the trial jet-energy-scale
// factor pa that makes Data Region A's mean(x_{J#gamma}) match the fixed Jet12-referenced
// MC mean, this picks pa by minimizing a bin-by-bin shape chi2 between Data's and the
// MC's shape-normalized xJ distributions (chi2 = sum of squared pulls between bin
// fractions - see referenceShape() below). Region A only, no purity correction and no
// Region C, exactly as in grid_insitu_jet12.C - this is the same non-purity-corrected
// Data-Region-A-vs-Jet12-Region-A consistency cross-check, just with the shape-sensitive
// fit criterion instead of the mean-matching one.
//
// Reuses grid_insitu_shapechi2.C's nXjBinsForChi2 low-count tail cutoff, coarseRebinPtBin
// pT-bin-2 rebinning, and 1/N errt floor verbatim (see that macro's comments for why each
// exists - debug_shapechi2_spike.C/debug_shapechi2_sawtooth.C) rather than inventing new
// bin-selection/uncertainty-combination logic for this cross-check (see
// gammajet_unfold/CLAUDE.md's bin-selection ground rule).

// insitu/ is split into inputs/ (the raw Data/Jet12 insitu ntuples, written by
// unfolder.h's production pipeline), output/ (this and the other grid_insitu*.C
// macros' own .root output), and pdfs/ (their .pdf output).
const char * insitu_input_dir  = "/home/samson72/sphnx/gammajet_unfold/insitu/inputs";
const char * insitu_output_dir = "/home/samson72/sphnx/gammajet_unfold/insitu/output";
const char * insitu_pdf_dir    = "/home/samson72/sphnx/gammajet_unfold/insitu/pdfs";

// Same restriction as grid_insitu.C/grid_insitu_jet12.C - only ana::ptBinsUsed
// (15-20, 20-25, 25-35 GeV).
const int nPtBinsUsed = ana::nPtBinsUsed;

// Same explicit, stated low-count-tail exclusion as grid_insitu_shapechi2.C - see that
// macro's comment (gammajet_unfold/CLAUDE.md's bin-selection ground rule).
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;

// Same pT-bin-2 coarse rebinning as grid_insitu_shapechi2.C, for the same low-Region-A-
// statistics reason (see that macro's comment) - Jet12's Region A sample size in this
// bin is comparable to grid_insitu.C's, so the same fix applies here.
const int coarseRebinPtBin = 2;
const int coarseGroupSize = 2;

// The MC reference sample - always "Jet12_long" now (the legacy "Jet12_full"
// cross-check has been retired - see grid_insitu_jet12.C's identical comment).
const string mcTrigger = "Jet12_long";

// Nominal cross-section weight for the Jet12 sample - see grid_insitu_jet12.C's
// identical comment for why this cancels out of every mean/shape computed below.
map<int,double> jet_scale = {{12,3.997e+06}};

// coarsenSum, coarsenQuadrature, referenceMeans, referenceShape, and
// computeRegionAMeans now live in src/insitu_utility.h/.cc (insitu_utility::
// namespace) - moved there after being found copy-pasted byte-for-byte across all six
// grid_insitu*.C macros (see that header's comment).

// One comparison page: top panel is mean(x_J) vs pT for MC and raw Data; bottom panel
// is the raw ratio (raw Data/MC) and the corrected ratio (best-fit-scaled Data/MC) -
// same layout as grid_insitu_jet12.C's drawJESPage(). The best-fit pa stamped here now
// comes from the shape chi2 scan below, not from matching this mean.
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
  jestext.DrawLatex(.18,.28, Form("Data to MC (Jet12) JES (shape #chi^{2}) = %.4f #pm %.4f", pa, paErr));

  c->SaveAs(pdfPath);
}

// One xJ-distribution comparison page, for a single photon-pT bin - unchanged from
// grid_insitu_jet12.C (still shape-normalized-plus-bin-width-density for display; the
// fit itself uses raw bin fractions, see referenceShape() above).
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
      Form("|#eta^{#gamma}|<%.1f, |#eta^{jet}|<%.1f", ana::etacut, ana::etacut-ana::JetRs[ir]),
      Form("#Delta#phi>%.0f#pi/%.0f", ana::oppnum, ana::oppden)
    }, .18, .85, 16, gPad->GetWh());

  c->SaveAs(pdfPath);
}

void grid_insitu_jet12_shapechi2(string systag = "nominal") {
  // See grid_insitu.C's identical comment: avoids per-pT-bin histograms auto-
  // registering into whichever radius subdirectory was left current by the previous
  // iteration's mkdir/cd.
  TH1::AddDirectory(kFALSE);

  // One file/one PDF for the whole systag, all seven jet radii inside - see
  // grid_insitu_jet12.C's identical comment.
  string pdfPathStr = Form("%s/grid_insitu_jet12_shapechi2_%s.pdf", insitu_pdf_dir, systag.c_str());
  TCanvas * c = new TCanvas("c","",700,700);
  c->SaveAs(Form("%s[", pdfPathStr.c_str()));

  string outfilename = Form("%s/grid_insitu_jet12_shapechi2_%s.root", insitu_output_dir, systag.c_str());
  TFile * fout = TFile::Open(outfilename.c_str(), "RECREATE");

  for (int ir = 0; ir < ana::nJetR; ir++) {

  const char * dataFile = Form("%s/Data_%s_insitu.root", insitu_input_dir, systag.c_str());
  vector<DataEvent> dataA = insitu_utility::cacheDataEvents(dataFile, 0, ir);
  cout << "Cached Data events: region A=" << dataA.size() << endl;

  // Low-xJ floor per used pT bin - same cut unfolder::check_pair applies at floorScale=1
  // before a reco jet enters hrecoxj/the response matrix (see src/insitu_utility.h's
  // lowXjFloor comment and grid_insitu_jet12.C's identical comment).
  float lowXj[nPtBinsUsed];
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) lowXj[ipt] = insitu_utility::lowXjFloor(ir, ana::ptBinsUsed[ipt]);

  vector<pair<string,double>> mcSamples = {
    {Form("%s/%s_pythia_%s_insitu.root", insitu_input_dir, mcTrigger.c_str(), systag.c_str()), jet_scale[12]},
  };

  float refMean[nPtBinsUsed], refMeanErr[nPtBinsUsed];
  insitu_utility::referenceMeans(mcSamples, 0, ir, refMean, refMeanErr, lowXj);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    cout << "Jet12 MC reference <x_J> pt bin " << ipt << " [" << ana::ptBinsUsed[ipt] << "," << ana::ptBinsUsed[ipt+1]
         << "): " << refMean[ipt] << " +/- " << refMeanErr[ipt] << endl;
  }

  // Fixed MC reference xJ SHAPE (bin fractions) - this is what the pa scan below
  // actually fits to, replacing refMean/refMeanErr as the fit target (same as
  // grid_insitu_shapechi2.C's relationship to grid_insitu.C).
  vector<vector<double>> refFrac, refFracErr;
  insitu_utility::referenceShape(mcSamples, 0, ir, refFrac, refFracErr, lowXj);

  // -----------------------------
  // Grid scan: single overall jet-energy-scale factor pa, no pT-dependence - same scan
  // window/step (insitu_utility.h's scanLow/scanHigh/scanN) as grid_insitu_jet12.C, but
  // the chi2 minimized at each point is now a shape chi2 (sum of squared pulls between
  // Data's and MC's shape-normalized xJ bin fractions), not the old mean(x_J)-matching
  // chi2 - same fit criterion as grid_insitu_shapechi2.C's Region-A block, Region A only
  // (no purity correction, no Region C - same restriction as grid_insitu_jet12.C).
  // -----------------------------
  const int na = insitu_utility::scanN;
  const float lowa = insitu_utility::scanLow, higha = insitu_utility::scanHigh;

  TGraph * gchisq = new TGraph(na);
  gchisq->SetName("gchisq_regionA_jet12ref");
  gchisq->SetTitle(";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});Shape #chi^{2}");

  float minchisq = FLT_MAX, minpa = 1;
  int ibest = 0;

  for (int ia = 0; ia < na; ia++) {
    float pa = lowa + ia*(higha-lowa)/na;

    // Per-(pT bin, xJ bin) raw (unweighted) Data counts at this trial pa - same as
    // grid_insitu_shapechi2.C's scan loop.
    vector<vector<double>> countA(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
    for (auto & ev : dataA) {
      float x = (ev.jet_pt/pa)/ev.pho_pt;
      if (x < lowXj[ev.ptbin]) continue;
      int ixj = ana::findUnfoldXjBin(x);
      if (ixj < 0 || ixj >= ana::nUnfoldXjBins) continue;
      countA[ev.ptbin][ixj] += 1;
    }

    float chisq = 0;
    for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
      double NA = 0;
      for (double n : countA[ipt]) NA += n;
      if (NA <= 0) continue;

      // pT bin 2 fits on coarsened bins (see coarseRebinPtBin above); every other pT
      // bin uses the native fine binning - same as grid_insitu_shapechi2.C.
      vector<double> useA, useRefFrac, useRefFracErr;
      int nBins;
      if (ipt == coarseRebinPtBin) {
        int startBin = ana::findUnfoldXjBin(lowXj[ipt]);
        useA          = insitu_utility::coarsenSum(countA[ipt], startBin, nXjBinsForChi2, coarseGroupSize);
        useRefFrac    = insitu_utility::coarsenSum(refFrac[ipt], startBin, nXjBinsForChi2, coarseGroupSize);
        useRefFracErr = insitu_utility::coarsenQuadrature(refFracErr[ipt], startBin, nXjBinsForChi2, coarseGroupSize);
        nBins = (int)useA.size();
      } else {
        useA = countA[ipt];
        useRefFrac = refFrac[ipt]; useRefFracErr = refFracErr[ipt];
        nBins = nXjBinsForChi2;
      }

      // Pull-squared per xJ bin between Data's raw bin fraction and the fixed MC
      // reference fraction, summed over the first nBins bins - same 1/NA errt floor as
      // grid_insitu_shapechi2.C's Region-A block (see that macro's comment for why).
      double errFloor = 1.0/NA;
      for (int ib = 0; ib < nBins; ib++) {
        double fData = useA[ib]/NA;
        double errData = sqrt(useA[ib])/NA;
        double errt = sqrt(errData*errData + useRefFracErr[ib]*useRefFracErr[ib]);
        errt = std::max(errt, errFloor);
        double diff = fData - useRefFrac[ib];
        chisq += diff*diff/(errt*errt);
      }
    }

    gchisq->SetPoint(ia, pa, chisq);
    if (chisq < minchisq) { minchisq = chisq; minpa = pa; ibest = ia; }
  }

  float errLow, errHigh;
  insitu_utility::findError(gchisq, ibest, minchisq, errLow, errHigh);

  cout << "\nFINAL RESULT (jet R=" << ana::JetRs[ir] << ", systag=" << systag << ", shape chi2, non-purity-corrected, "
       << mcTrigger << " reference)\n";
  cout << "Region A (Data) vs. Region A (" << mcTrigger << " MC):  p_a = " << minpa
       << " +" << errHigh << "/-" << errLow << " (chi2=" << minchisq << ")" << endl;

  // -----------------------------
  // Build x_J histograms per photon-pT bin: the fixed Jet12 MC reference, Data at
  // pa=1 (raw) and at the scan's best-fit pa (corrected) - same non-uniform binning
  // as grid_insitu_jet12.C.
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
  // Mean(x_J) vs pT comparison plot - display only now (the fit no longer targets this
  // mean) - top: MC vs raw Data; bottom: raw ratio vs corrected ratio (evaluated at the
  // scan's shape-chi2 best-fit pa).
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
  // Save - see grid_insitu.C's identical comment.
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

  TTree * wt = new TTree("results", "best-fit jet energy scale result (shape chi2, non-purity-corrected, Jet12 reference)");
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
