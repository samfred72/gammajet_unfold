#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
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

// See the R__LOAD_LIBRARY comment in grid_insitu.C - same reasoning applies here.
R__LOAD_LIBRARY(libgammajet_unfold.so);
R__LOAD_LIBRARY(libRooUnfold.so);

// Unfolded extension of grid_insitu.C's purity-corrected in-situ JES study.
//
// grid_insitu.C finds the single overall jet-energy-scale factor pa (jet_pt_corrected =
// jet_pt/pa) that makes Data's RECO-level mean(x_{J#gamma}) match Pythia8's RECO-level
// mean, both for Region A alone and for the purity-corrected A-minus-background
// combination - see that file's header for the full rationale. This macro asks the
// physically sharper version of the same question: instead of comparing two reco-level
// means (which both carry the same detector smearing and so can agree even if the
// underlying truth-level scale is off), it UNFOLDS the purity-corrected Data spectrum
// through the nominal photon+jet response matrix and compares the unfolded mean(x_J) to
// the fixed Pythia8 TRUTH-level mean, per photon-pT bin. This is the same background
// subtraction as grid_insitu.C's "purity-corrected" branch, plus the same
// unfold-before-compare step already used (at fixed nominal JES only) by
// drawing/draw_purity_corrected.C - here it's repeated at every trial pa in the scan.
//
// Region A and Region C jets are scaled by the same trial pa at every grid point, purity
// is computed once and held fixed (unaffected by jet energy scale), and the response
// matrix itself is built at the nominal (uncorrected) MC jet energy scale and held fixed
// across the whole scan - only the measured (Data) input changes with pa. This mirrors
// exactly how grid_insitu.C holds its own MC reference fixed while scanning Data alone.

// insitu/ is split into inputs/ (the raw Data insitu ntuple, written by unfolder.h's
// production pipeline), output/ (this and the other grid_insitu*.C macros' own .root
// output), and pdfs/ (their .pdf output).
const char * insitu_input_dir  = "/home/samson72/sphnx/gammajet_unfold/insitu/inputs";
const char * insitu_output_dir = "/home/samson72/sphnx/gammajet_unfold/insitu/output";
const char * insitu_pdf_dir    = "/home/samson72/sphnx/gammajet_unfold/insitu/pdfs";

const int nPtBinsUsed = ana::nPtBinsUsed;

// Bayesian-unfolding iteration count - same "best-iteration scan result" choice
// drawing/draw_purity_corrected.C uses (see that file's niterate comment /
// pdfs/purity_corrected_iterations_nominal.pdf), reused here rather than re-deriving it,
// since it comes from the same response matrix and the same purity-corrected input.
const int niterate = 2;

// struct DataEvent and cacheDataEvents now live in src/insitu_utility.h/.cc - unlike
// grid_insitu.C's default (restrictToUsed=true) usage, this file calls
// insitu_utility::cacheDataEvents(..., false) below to keep every ana::ptBins bin
// (0..nPtBins-1, including the low-pT migration buffer and high-pT overflow bins), not
// just the nPtBinsUsed reported ones - the response matrix's flattened (pT,xJ) measured
// vector needs a complete input for cross-pT-bin migration to unfold correctly, the same
// reason drawing/draw_purity_corrected.C's buildFullyCorrected() purity-corrects all
// ana::nPtBins slices instead of just the used ones.

// buildXjByPtBin and purityCorrectByPtBin now live in src/insitu_utility.h/.cc
// (insitu_utility:: namespace) - moved there after being found copy-pasted (differing
// only in bin count / zero-vs-real purity-error arrays) across all six
// grid_insitu*.C macros (see that header's comment). This macro always calls them
// ana::nPtBins-sized with zero-filled purity-error arrays, since purity is held fixed
// across the whole pa scan (file header) and the asymmetric purity-uncertainty term
// isn't needed here.

// Purity-correct (all nPtBins), reflatten into the response matrix's dimensionality,
// unfold through the fixed nominal-JES response, and return mean(x_J)/error per USED
// photon-pT bin (index 0..nPtBinsUsed-1, offset from ana::firstUsedPtBin) - the unfolded
// analogue of grid_insitu.C's computeCorrectedMeans(). If unfoldedOut is non-null, it is
// filled with clones of the per-used-pT-bin unfolded x_J histograms (caller owns them).
void computeUnfoldedMeans(const vector<DataEvent> & dataA, const vector<DataEvent> & dataC, float pa,
    const float purity[], const float purityC[], RooUnfoldResponse * response, TH1D * respRecoTemplate,
    float mean[], float err[], const float lowXj[], vector<TH1D*> * unfoldedOut = nullptr) {
  vector<TH1D*> hA    = insitu_utility::buildXjByPtBin(dataA, pa, ana::nPtBins, "hUnfA_tmp", lowXj);
  vector<TH1D*> hC    = insitu_utility::buildXjByPtBin(dataC, pa, ana::nPtBins, "hUnfC_tmp", lowXj);
  // Purity is held fixed across the whole pa scan (file header), so the asymmetric
  // purity-uncertainty term isn't needed here - zero error arrays.
  float zeroErr[ana::nPtBins] = {0};
  vector<TH1D*> hCorr = insitu_utility::purityCorrectByPtBin(hA, hC, ana::nPtBins,
      purity, zeroErr, zeroErr, purityC, zeroErr, zeroErr, "hUnfCorr_tmp");

  TH1D * flatCorrected = (TH1D*)respRecoTemplate->Clone("flatCorrected_tmp");
  flatCorrected->Reset("ICES");
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) unfold_utility::reflattenXj(hCorr[ipt], ipt, flatCorrected);

  TH1D * flatUnfolded = unfold_utility::unfoldOnce(response, flatCorrected, niterate, "flatUnfolded_tmp");

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

// findError, meanGraph, ratioGraph, drawSPhenixLabel now live in src/insitu_utility.h/.cc.

// One comparison page: top panel is unfolded mean(x_J) vs pT for Truth (fixed) and raw
// (pa=1) unfolded Data; bottom panel is the raw ratio and the corrected ratio (evaluated
// at the scan's best-fit pa) - the latter should sit flat at 1 by construction. Same
// layout as grid_insitu.C's drawJESPage, with "Data/MC" -> "Unfolded/Truth".
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
      Form("|#eta^{#gamma}|<%.1f, |#eta^{jet}|<%.1f", ana::etacut, ana::etacut-ana::JetRs[ir]),
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
  jestext.DrawLatex(.18,.28, Form("Data to MC JES (unfolded vs. truth) = %.4f #pm %.4f", pa, paErr));

  c->SaveAs(pdfPath);
}

// x_J shape comparison for one photon-pT bin: fixed truth reference, unfolded Data at
// raw (pa=1) JES, and unfolded Data at the scan's best-fit pa - all shape-normalized
// (density, then unit-area) since reco/unfolded/truth sit at different absolute scales
// from reconstruction efficiency, same convention as
// drawing/draw_purity_corrected.C's page 2.
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
      Form("|#eta^{#gamma}|<%.1f, |#eta^{jet}|<%.1f", ana::etacut, ana::etacut-ana::JetRs[ir]),
      Form("#Delta#phi>%.0f#pi/%.0f", ana::oppnum, ana::oppden)
    }, .18, .85, 16, gPad->GetWh());

  c->SaveAs(pdfPath);
}

void grid_insitu_unfolded(string systag = "nominal", int na = insitu_utility::scanN) {
  // Newly created histograms are not registered to any TDirectory, so the ~20
  // temporaries allocated per pa grid point (computeUnfoldedMeans, called na times) don't
  // pile up in gROOT's object list or collide on name across iterations - they're freed
  // by their own explicit `delete` calls instead. This also means the per-pT-bin
  // histograms that ARE meant to be saved below don't auto-register into whichever
  // radius subdirectory was left current by the previous iteration's mkdir/cd (see
  // grid_insitu.C's identical comment) - both are covered by this one call.
  TH1::AddDirectory(kFALSE);

  // Response-matrix source (see below) and output file/PDF are shared across every
  // radius - constructed/opened once here, before the per-radius loop, instead of
  // per-radius as before.
  drawer d("pythia", systag);

  string pdfPathStr = Form("%s/grid_insitu_unfolded_%s.pdf", insitu_pdf_dir, systag.c_str());
  TCanvas * c = new TCanvas("c","",700,700);
  c->SaveAs(Form("%s[", pdfPathStr.c_str()));

  string outfilename = Form("%s/grid_insitu_unfolded_%s.root", insitu_output_dir, systag.c_str());
  TFile * fout = TFile::Open(outfilename.c_str(), "RECREATE");

  for (int ir = 0; ir < ana::nJetR; ir++) {

  string dataFile = insitu_utility::insituFilename(insitu_input_dir, "Data", "", systag);
  vector<DataEvent> dataA = insitu_utility::cacheDataEvents(dataFile.c_str(), 0, ir, false);
  vector<DataEvent> dataC = insitu_utility::cacheDataEvents(dataFile.c_str(), 2, ir, false);
  cout << "Cached Data events: region A=" << dataA.size() << " region C=" << dataC.size() << endl;

  // Low-xJ floor per ana::ptBins bin (all nPtBins, since buildXjByPtBin fills every one
  // of them, not just the used bins - see the cacheDataEvents comment above) - same cut
  // unfolder::check_pair applies at floorScale=1 before a reco jet enters hrecoxj/the
  // response matrix (see src/insitu_utility.h's lowXjFloor comment).
  float lowXj[ana::nPtBins];
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) lowXj[ipt] = insitu_utility::lowXjFloor(ir, ana::ptBins[ipt]);

  // Purity per ana::ptBins bin (all nPtBins, not just the used ones - see
  // cacheDataEvents comment) - computed once, held fixed across the whole pa scan, same
  // as grid_insitu.C.
  float purity[ana::nPtBins], purityC[ana::nPtBins];
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    purity[ipt]  = ana::getPurity(ana::ptBins[ipt], ana::ptBins[ipt+1], systag, ir);
    purityC[ipt] = ana::getPurityC(ana::ptBins[ipt], ana::ptBins[ipt+1], systag, ir);
    cout << "Purity pt bin " << ipt << " [" << ana::ptBins[ipt] << "," << ana::ptBins[ipt+1]
         << "): P_A=" << purity[ipt] << " P_C=" << purityC[ipt] << endl;
  }

  // Response matrix: the physically meaningful photon+jet MC response (cross-section-
  // weighted combination of Photon5/10/20, type=1/isample=-1), built at the nominal
  // (uncorrected) MC jet energy scale and held fixed across the whole pa scan - only the
  // trial-pa-rescaled Data measured spectrum changes below. Same source/convention as
  // drawing/draw_purity_corrected.C's response matrix. (drawer d itself is shared/
  // hoisted above the radius loop - see top of function.)
  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i", ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i", ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i", ir), 1);
  RooUnfoldResponse * response = new RooUnfoldResponse(respRecoTemplate, respTruthTemplate, respMatrix2D);

  // Fixed Pythia8 gamma+jet TRUTH-level mean(x_J) per used photon-pT bin - the reference
  // the unfolded Data mean is compared against below (never rescaled: JES is a
  // reco-level detector effect, truth is untouched by it).
  float truthMean[nPtBinsUsed], truthMeanErr[nPtBinsUsed];
  for (int k = 0; k < nPtBinsUsed; k++) {
    int ipt = ana::firstUsedPtBin + k;
    TH1D * hTruth = unfold_utility::unflattenXj(respTruthTemplate, ipt, "hTruthMean_tmp");
    truthMean[k] = hTruth->GetMean();
    truthMeanErr[k] = hTruth->GetMeanError();
    cout << "Truth <x_J> pt bin " << k << " [" << ana::ptBinsUsed[k] << "," << ana::ptBinsUsed[k+1]
         << "): " << truthMean[k] << " +/- " << truthMeanErr[k] << endl;
    delete hTruth;
  }

  // -----------------------------
  // Grid scan: single overall jet-energy-scale factor pa, no pT-dependence - same scan
  // window (insitu_utility.h's scanLow/scanHigh) as grid_insitu.C. At each pa,
  // purity-correct region A/C (all ana::nPtBins slices), unfold through the fixed
  // response above, and chi2 the unfolded mean(x_J) (used bins only) against the fixed
  // truth mean.
  // -----------------------------
  const float lowa = insitu_utility::scanLow, higha = insitu_utility::scanHigh;

  TGraph * gchisqUnfold = new TGraph(na);
  gchisqUnfold->SetName("gchisq_unfolded");
  gchisqUnfold->SetTitle(";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});#chi^{2}");

  float minchisqUnfold = FLT_MAX, minpaUnfold = 1;
  int ibestUnfold = 0;

  for (int ia = 0; ia < na; ia++) {
    float pa = lowa + ia*(higha-lowa)/na;

    float mean[nPtBinsUsed], err[nPtBinsUsed];
    computeUnfoldedMeans(dataA, dataC, pa, purity, purityC, response, respRecoTemplate, mean, err, lowXj);

    float chisq = 0;
    for (int k = 0; k < nPtBinsUsed; k++) {
      if (truthMean[k] <= 0) continue;
      double diff = 1 - mean[k]/truthMean[k];
      double errt = sqrt((err[k]*err[k])/(truthMean[k]*truthMean[k])
          + mean[k]*mean[k]*truthMeanErr[k]*truthMeanErr[k]/pow(truthMean[k],4));
      if (errt > 0) chisq += diff*diff/(errt*errt);
    }

    gchisqUnfold->SetPoint(ia, pa, chisq);
    if (chisq < minchisqUnfold) { minchisqUnfold = chisq; minpaUnfold = pa; ibestUnfold = ia; }

    if (ia % 100 == 0) cout << "  scan " << ia << "/" << na << ": pa=" << pa << " chi2=" << chisq << endl;
  }

  float errLowUnfold, errHighUnfold;
  insitu_utility::findError(gchisqUnfold, ibestUnfold, minchisqUnfold, errLowUnfold, errHighUnfold);

  cout << "\nFINAL RESULT (jet R=" << ana::JetRs[ir] << ", systag=" << systag << ", unfolded vs. truth)\n";
  cout << "Purity-corrected + unfolded:  p_a = " << minpaUnfold
       << " +" << errHighUnfold << "/-" << errLowUnfold << " (chi2=" << minchisqUnfold << ")" << endl;

  // -----------------------------
  // Build final comparison histograms/graphs at pa=1 (raw) and pa=minpaUnfold
  // (best-fit), and truth - one page for the mean-vs-pT summary, one per used pT bin
  // for the xJ shape comparison.
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
  // Save - see grid_insitu.C's identical comment.
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

  TTree * wt = new TTree("results", "best-fit jet energy scale result (unfolded vs. truth)");
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
