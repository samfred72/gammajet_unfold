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

// See the R__LOAD_LIBRARY comment in grid_insitu.C - same reasoning applies here.
R__LOAD_LIBRARY(libgammajet_unfold.so);
R__LOAD_LIBRARY(libRooUnfold.so);

// Shape-chi2 variant of grid_insitu_unfolded.C, itself the unfolded extension of
// grid_insitu.C - see those two files' headers for the full rationale of purity-
// correcting and unfolding before comparing to Pythia8 TRUTH-level x_{J#gamma}.
// grid_insitu_unfolded.C picks the trial jet-energy-scale factor pa that makes the
// UNFOLDED Data mean(x_J) match the fixed truth-level mean - a single number per pT
// bin, blind to any shape difference that happens to preserve the mean (exactly the
// same limitation grid_insitu_shapechi2.C addresses at reco level). This version
// instead picks pa by minimizing a bin-by-bin shape chi2 between the unfolded Data xJ
// distribution and the fixed Pythia8 TRUTH-level xJ distribution, each shape-normalized
// to unit area (bin fraction, not bin-width density - same reasoning as
// grid_insitu_shapechi2.C's referenceShape()), summed over xJ bins (dropping the
// low-stat tail) and over used pT bins. Region A and Region C are still purity-
// corrected and unfolded through the fixed nominal-JES response exactly as in
// grid_insitu_unfolded.C - only the fit criterion changes, so the two best-fit pa
// values can be compared directly.

// insitu/ is split into inputs/ (the raw Data insitu ntuple, written by unfolder.h's
// production pipeline), output/ (this and the other grid_insitu*.C macros' own .root
// output), and pdfs/ (their .pdf output).
const char * insitu_input_dir  = ana::path("insitu/inputs");
const char * insitu_output_dir = ana::path("insitu/output");
const char * insitu_pdf_dir    = ana::path("insitu/pdfs");

const int nPtBinsUsed = ana::nPtBinsUsed;

// Bayesian-unfolding iteration count - same choice as grid_insitu_unfolded.C /
// drawing/draw_purity_corrected.C.
const int niterate = 2;

// The last 3 xJ bins in each pT bin have very low counts, so a per-bin shape chi2
// computed against them is dominated by their noise rather than genuine shape
// agreement - excluded here, same explicit-threshold precedent as
// drawing/draw_purity_corrected.C's nXjBinsForChi2 and grid_insitu_shapechi2.C's own
// (see gammajet_unfold/CLAUDE.md's bin-selection ground rule).
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;

// struct DataEvent and cacheDataEvents live in src/insitu_utility.h/.cc - this file
// calls insitu_utility::cacheDataEvents(..., false) below to keep every ana::ptBins bin
// (not just the nPtBinsUsed reported ones), same reason as grid_insitu_unfolded.C: the
// response matrix's flattened (pT,xJ) measured vector needs a complete input for
// cross-pT-bin migration to unfold correctly.

// buildXjByPtBin and purityCorrectByPtBin now live in src/insitu_utility.h/.cc
// (insitu_utility:: namespace) - see grid_insitu_unfolded.C's identical comment.

// Purity-correct (all nPtBins), reflatten, unfold through the fixed nominal-JES
// response, and return mean(x_J)/error per used photon-pT bin - unchanged from
// grid_insitu_unfolded.C. Kept only for the mean(x_J) display page now; the fit
// criterion below uses computeUnfoldedShapeChi2 instead.
void computeUnfoldedMeans(const vector<DataEvent> & dataA, const vector<DataEvent> & dataC, float pa,
    const float purity[], const float purityC[], RooUnfoldResponse * response, TH1D * respRecoTemplate,
    float mean[], float err[], const float lowXj[], vector<TH1D*> * unfoldedOut = nullptr) {
  vector<TH1D*> hA    = insitu_utility::buildXjByPtBin(dataA, pa, ana::nPtBins, "hUnfA_tmp", lowXj);
  vector<TH1D*> hC    = insitu_utility::buildXjByPtBin(dataC, pa, ana::nPtBins, "hUnfC_tmp", lowXj);
  // Purity is held fixed across the whole pa scan (file header) - zero error arrays,
  // see grid_insitu_unfolded.C's identical comment.
  float zeroErr[ana::nPtBins] = {0};
  vector<TH1D*> hCorr = insitu_utility::purityCorrectByPtBin(hA, hC, ana::nPtBins,
      purity, zeroErr, zeroErr, purityC, zeroErr, zeroErr, "hUnfCorr_tmp");

  TH1D * flatCorrected = (TH1D*)respRecoTemplate->Clone("flatCorrected_tmp");
  flatCorrected->Reset("ICES");
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) unfold_utility::reflattenXj(hCorr[ipt], ipt, flatCorrected);

  // includeSystematics=false - see computeUnfoldedShapeChi2's identical comment below;
  // this copy is only used for the (twice-per-radius, not per-scan-point) display page.
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

// Same internal pipeline as computeUnfoldedMeans (purity-correct, reflatten, unfold),
// but instead of extracting mean(x_J), computes a shape chi2 between the unfolded xJ
// bin fractions and the fixed Pythia8 TRUTH-level bin fractions (truthFrac/
// truthFracErr, built once outside the pa scan - see the caller), summed over the
// first nXjBinsForChi2 xJ bins and over all used pT bins. This is the fit criterion
// the pa scan below actually minimizes.
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

  // includeSystematics=false: this is the actual fit criterion below, called once per
  // pa-scan point (na=insitu_utility::scanN x ana::nJetR calls total) - the expensive
  // response-matrix-statistics covariance term (see unfold_utility.h's includeSystematics
  // comment, ~3s/call) made the full scan take hours. errUnf below still reflects the
  // (cheap, always-on) data-statistics covariance term, just not the response-matrix
  // contribution - acceptable here since this method is a cross-check against
  // grid_insitu.C's plain purity-corrected fit (the one draw_jes_summary.C actually
  // sources ana::jesNominal from), not the headline result.
  TH1D * flatUnfolded = unfold_utility::unfoldOnce(response, flatCorrected, niterate, "flatUnfoldedShape_tmp", false);

  float chisq = 0;
  for (int k = 0; k < nPtBinsUsed; k++) {
    int ipt = ana::firstUsedPtBin + k;
    TH1D * hU = unfold_utility::unflattenXj(flatUnfolded, ipt, "hUnfoldedShapePt_tmp");
    double N = hU->Integral();
    if (N > 0) {
      // errt floored at 1/N - same fix and same reason as grid_insitu_shapechi2.C's
      // scan loop (see its comment): a bin that unfolds to ~0 content can still have
      // errUnf~0, which understates the real uncertainty and lets one near-empty bin
      // dominate the chi2 (see debug_shapechi2_spike.C, written against the Region-A
      // version of this same pathology).
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

// findError, meanGraph, ratioGraph, drawSPhenixLabel live in src/insitu_utility.h/.cc.

// One comparison page: top panel is unfolded mean(x_J) vs pT for Truth (fixed) and raw
// (pa=1) unfolded Data; bottom panel is the raw ratio and the corrected ratio (evaluated
// at the scan's shape-chi2 best-fit pa) - unchanged from grid_insitu_unfolded.C, display
// only now (the fit no longer targets this mean, so the corrected ratio isn't
// guaranteed to sit exactly at 1 - a useful cross-check of how much shape-matching
// still brings the mean into line).
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
  jestext.DrawLatex(.18,.28, Form("Data to MC JES (unfolded vs. truth, shape #chi^{2}) = %.4f #pm %.4f", pa, paErr));

  c->SaveAs(pdfPath);
}

// x_J shape comparison for one photon-pT bin - unchanged from grid_insitu_unfolded.C.
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

void grid_insitu_unfolded_shapechi2(string systag = "nominal", int na = insitu_utility::scanN) {
  // Newly created histograms are not registered to any TDirectory, so the temporaries
  // allocated per pa grid point (computeUnfoldedShapeChi2, called na times) don't pile
  // up in gROOT's object list or collide on name across iterations - and, per
  // grid_insitu.C's identical comment, this also keeps the per-pT-bin histograms that
  // ARE meant to be saved from auto-registering into whichever radius subdirectory was
  // left current by the previous iteration's mkdir/cd below.
  TH1::AddDirectory(kFALSE);

  // Response-matrix source and output file/PDF are shared across every radius -
  // constructed/opened once here, before the per-radius loop.
  drawer d("pythia", systag);

  string pdfPathStr = Form("%s/grid_insitu_unfolded_shapechi2_%s.pdf", insitu_pdf_dir, systag.c_str());
  TCanvas * c = new TCanvas("c","",700,700);
  c->SaveAs(Form("%s[", pdfPathStr.c_str()));

  string outfilename = Form("%s/grid_insitu_unfolded_shapechi2_%s.root", insitu_output_dir, systag.c_str());
  TFile * fout = TFile::Open(outfilename.c_str(), "RECREATE");

  for (int ir = 0; ir < ana::nJetR; ir++) {

  string dataFile = insitu_utility::insituFilename(insitu_input_dir, "Data", "", systag);
  vector<DataEvent> dataA = insitu_utility::cacheDataEvents(dataFile.c_str(), 0, ir, false);
  vector<DataEvent> dataC = insitu_utility::cacheDataEvents(dataFile.c_str(), 2, ir, false);
  cout << "Cached Data events: region A=" << dataA.size() << " region C=" << dataC.size() << endl;

  // Low-xJ floor per ana::ptBins bin (all nPtBins, since buildXjByPtBin fills every one
  // of them) - same cut unfolder::check_pair applies at floorScale=1 before a reco jet
  // enters hrecoxj/the response matrix (see src/insitu_utility.h's lowXjFloor comment).
  float lowXj[ana::nPtBins];
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) lowXj[ipt] = insitu_utility::lowXjFloor(ir, ana::ptBins[ipt]);

  // Purity per ana::ptBins bin (all nPtBins) - computed once, held fixed across the
  // whole pa scan, same as grid_insitu_unfolded.C.
  float purity[ana::nPtBins], purityC[ana::nPtBins];
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    purity[ipt]  = ana::getPurity(ana::ptBins[ipt], ana::ptBins[ipt+1], systag, ir);
    purityC[ipt] = ana::getPurityC(ana::ptBins[ipt], ana::ptBins[ipt+1], systag, ir);
    cout << "Purity pt bin " << ipt << " [" << ana::ptBins[ipt] << "," << ana::ptBins[ipt+1]
         << "): P_A=" << purity[ipt] << " P_C=" << purityC[ipt] << endl;
  }

  // Response matrix: same construction as grid_insitu_unfolded.C. (drawer d is shared/
  // hoisted above the radius loop - see top of function.)
  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i", ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i", ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i", ir), 1);
  RooUnfoldResponse * response = new RooUnfoldResponse(respRecoTemplate, respTruthTemplate, respMatrix2D);

  // Fixed Pythia8 gamma+jet TRUTH-level mean(x_J) per used photon-pT bin - kept only
  // for the mean(x_J) display page now (see truthFrac/truthFracErr below for the
  // actual fit target).
  float truthMean[nPtBinsUsed], truthMeanErr[nPtBinsUsed];
  // Fixed Pythia8 gamma+jet TRUTH-level xJ SHAPE (bin fraction) and its per-bin error,
  // per used photon-pT bin - the shape-chi2 analogue of truthMean/truthMeanErr, and
  // the actual fit target below. respTruthTemplate already carries the correct
  // (weighted-MC) Sumw2 bin errors, so no separate sumw/sumw2 accumulation is needed
  // here (unlike grid_insitu_shapechi2.C's referenceShape(), which has to build this
  // from raw insitutree events instead of an existing histogram) - per-bin fraction
  // error uses the same treat-N-as-fixed simplification as everywhere else in this
  // codebase (err(f_i) = binError_i/N).
  vector<vector<double>> truthFrac(nPtBinsUsed), truthFracErr(nPtBinsUsed);
  for (int k = 0; k < nPtBinsUsed; k++) {
    int ipt = ana::firstUsedPtBin + k;
    TH1D * hTruth = unfold_utility::unflattenXj(respTruthTemplate, ipt, "hTruthShape_tmp");
    truthMean[k] = hTruth->GetMean();
    truthMeanErr[k] = hTruth->GetMeanError();
    double N = hTruth->Integral();
    truthFrac[k].assign(ana::nUnfoldXjBins, 0.);
    truthFracErr[k].assign(ana::nUnfoldXjBins, 0.);
    if (N > 0) {
      for (int ixj = 1; ixj <= ana::nUnfoldXjBins; ixj++) {
        truthFrac[k][ixj-1]    = hTruth->GetBinContent(ixj)/N;
        truthFracErr[k][ixj-1] = hTruth->GetBinError(ixj)/N;
      }
    }
    cout << "Truth <x_J> pt bin " << k << " [" << ana::ptBinsUsed[k] << "," << ana::ptBinsUsed[k+1]
         << "): " << truthMean[k] << " +/- " << truthMeanErr[k] << endl;
    delete hTruth;
  }

  // -----------------------------
  // Grid scan: single overall jet-energy-scale factor pa, no pT-dependence - same scan
  // window (insitu_utility.h's scanLow/scanHigh) as grid_insitu_unfolded.C. At each pa,
  // purity-correct region A/C (all ana::nPtBins slices), unfold through the fixed
  // response above, and shape-chi2 the unfolded xJ bin fractions (used bins only,
  // dropping the low-stat tail) against the fixed truth shape.
  // -----------------------------
  const float lowa = insitu_utility::scanLow, higha = insitu_utility::scanHigh;

  TGraph * gchisqUnfold = new TGraph(na);
  gchisqUnfold->SetName("gchisq_unfolded");
  gchisqUnfold->SetTitle(";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});Shape #chi^{2}");

  float minchisqUnfold = FLT_MAX, minpaUnfold = 1;
  int ibestUnfold = 0;

  for (int ia = 0; ia < na; ia++) {
    float pa = lowa + ia*(higha-lowa)/na;

    float chisq = computeUnfoldedShapeChi2(dataA, dataC, pa, purity, purityC, response, respRecoTemplate,
        truthFrac, truthFracErr, lowXj);

    gchisqUnfold->SetPoint(ia, pa, chisq);
    if (chisq < minchisqUnfold) { minchisqUnfold = chisq; minpaUnfold = pa; ibestUnfold = ia; }

    if (ia % 100 == 0) cout << "  scan " << ia << "/" << na << ": pa=" << pa << " chi2=" << chisq << endl;
  }

  float errLowUnfold, errHighUnfold;
  insitu_utility::findError(gchisqUnfold, ibestUnfold, minchisqUnfold, errLowUnfold, errHighUnfold);

  cout << "\nFINAL RESULT (jet R=" << ana::JetRs[ir] << ", systag=" << systag << ", unfolded vs. truth, shape chi2)\n";
  cout << "Purity-corrected + unfolded:  p_a = " << minpaUnfold
       << " +" << errHighUnfold << "/-" << errLowUnfold << " (chi2=" << minchisqUnfold << ")" << endl;

  // -----------------------------
  // Build final comparison histograms/graphs at pa=1 (raw) and pa=minpaUnfold
  // (best-fit), and truth - unchanged from grid_insitu_unfolded.C, using the
  // (still-available) computeUnfoldedMeans for display purposes.
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

  TTree * wt = new TTree("results", "best-fit jet energy scale result (unfolded vs. truth, shape chi2)");
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
