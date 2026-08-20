#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
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

const char * insitu_dir = "/home/samson72/sphnx/gammajet_unfold/insitu";

const int nPtBinsUsed = ana::nPtBinsUsed;

// The insitu trees are R=0.4 only (unfolder.cc fills insitu_tree under an explicit
// `ir == 2` gate) and the response matrix retrieved below must be the R=0.4 one to
// correspond to the same jets - unlike grid_insitu.C's ir function parameter (which only
// ever makes sense at 2 for the same reason but is passed through anyway), ir is fixed
// here since nothing else in this file can meaningfully vary it.
const int ir = 2;

// Bayesian-unfolding iteration count - same "best-iteration scan result" choice
// drawing/draw_purity_corrected.C uses (see that file's niterate comment /
// pdfs/purity_corrected_iterations_nominal.pdf), reused here rather than re-deriving it,
// since it comes from the same response matrix and the same purity-corrected input.
const int niterate = 2;

struct DataEvent { float pho_pt, jet_pt; int ptbin; };

// Unlike grid_insitu.C's cacheDataEvents, this keeps every ana::ptBins bin (0..nPtBins-1,
// including the low-pT migration buffer and high-pT overflow bins), not just the
// nPtBinsUsed reported ones - the response matrix's flattened (pT,xJ) measured vector
// needs a complete input for cross-pT-bin migration to unfold correctly, the same reason
// drawing/draw_purity_corrected.C's buildFullyCorrected() purity-corrects all
// ana::nPtBins slices instead of just the used ones.
vector<DataEvent> cacheDataEvents(const char * filename, int abcdSelect) {
  vector<DataEvent> events;
  TFile * f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie()) {
    cout << "WARNING: could not open " << filename << endl;
    return events;
  }
  TTree * t = (TTree*)f->Get("insitutree");
  Float_t pho_pt, jet_pt;
  Int_t abcd;
  t->SetBranchAddress("pho_pt", &pho_pt);
  t->SetBranchAddress("jet_pt", &jet_pt);
  t->SetBranchAddress("abcd", &abcd);
  Long64_t nentries = t->GetEntries();
  for (Long64_t e = 0; e < nentries; e++) {
    t->GetEntry(e);
    if (abcd != abcdSelect) continue;
    int ipt = ana::findPtBin(pho_pt);
    if (ipt < 0) continue;
    events.push_back({pho_pt, jet_pt, ipt});
  }
  f->Close();
  return events;
}

// x_J histogram per ana::ptBins bin (all nPtBins, ana::unfoldXjBins binning - same
// per-pT-bin layout unfold_utility::unflattenXj/reflattenXj use), at a trial
// jet-energy-scale factor pa.
vector<TH1D*> buildXjByPtBin(const vector<DataEvent> & data, float pa, const char * prefix) {
  vector<TH1D*> h(ana::nPtBins);
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    h[ipt] = new TH1D(Form("%s_pt%d", prefix, ipt), ";x_{J#gamma};Counts", ana::nUnfoldXjBins, ana::unfoldXjBins);
  }
  for (auto & ev : data) {
    h[ev.ptbin]->Fill((ev.jet_pt/pa)/ev.pho_pt);
  }
  return h;
}

// Purity-correct region A/C per ana::ptBins bin: A - (1-P)*(N_A/N_C)*C, same formula as
// grid_insitu.C's purityCorrectByPtBin / drawing/draw_purity_corrected.C's
// purityCorrectP (minus that file's separate purity-uncertainty term, which
// grid_insitu.C's pa scan already omits too - purity is held fixed across the scan, see
// file header). Falls back to raw region A (scale=0) when region C is empty in a pT
// slice - only the unused buffer bins are ever at risk of this.
vector<TH1D*> purityCorrectByPtBin(const vector<TH1D*> & hA, const vector<TH1D*> & hC,
    const float purity[], const char * prefix) {
  vector<TH1D*> h(ana::nPtBins);
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    h[ipt] = (TH1D*)hA[ipt]->Clone(Form("%s_pt%d", prefix, ipt));
    float NA = hA[ipt]->Integral(), NC = hC[ipt]->Integral();
    float scale = (NC > 0) ? (1-purity[ipt])*(NA/NC) : 0;
    for (int ib = 1; ib <= h[ipt]->GetNbinsX(); ib++) {
      float a = hA[ipt]->GetBinContent(ib), ae = hA[ipt]->GetBinError(ib);
      float c = hC[ipt]->GetBinContent(ib), ce = hC[ipt]->GetBinError(ib);
      h[ipt]->SetBinContent(ib, a - scale*c);
      h[ipt]->SetBinError(ib, sqrt(ae*ae + pow(scale*ce,2)));
    }
  }
  return h;
}

// Purity-correct (all nPtBins), reflatten into the response matrix's dimensionality,
// unfold through the fixed nominal-JES response, and return mean(x_J)/error per USED
// photon-pT bin (index 0..nPtBinsUsed-1, offset from ana::firstUsedPtBin) - the unfolded
// analogue of grid_insitu.C's computeCorrectedMeans(). If unfoldedOut is non-null, it is
// filled with clones of the per-used-pT-bin unfolded x_J histograms (caller owns them).
void computeUnfoldedMeans(const vector<DataEvent> & dataA, const vector<DataEvent> & dataC, float pa,
    const float purity[], RooUnfoldResponse * response, TH1D * respRecoTemplate,
    float mean[], float err[], vector<TH1D*> * unfoldedOut = nullptr) {
  vector<TH1D*> hA    = buildXjByPtBin(dataA, pa, "hUnfA_tmp");
  vector<TH1D*> hC    = buildXjByPtBin(dataC, pa, "hUnfC_tmp");
  vector<TH1D*> hCorr = purityCorrectByPtBin(hA, hC, purity, "hUnfCorr_tmp");

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

// Scans outward from the minimum on a chi2-vs-pa graph for the two Delta-chi2=1 points -
// identical to grid_insitu.C's findError().
void findError(TGraph * g, int ibest, float minchisq, float & errLow, float & errHigh) {
  double xbest, ytmp;
  g->GetPoint(ibest, xbest, ytmp);
  errLow = xbest - g->GetX()[0];
  errHigh = g->GetX()[g->GetN()-1] - xbest;
  double x, y;
  for (int i = ibest; i >= 0; i--) {
    g->GetPoint(i, x, y);
    if (y - minchisq > 1.0) { errLow = xbest - x; break; }
  }
  for (int i = ibest; i < g->GetN(); i++) {
    g->GetPoint(i, x, y);
    if (y - minchisq > 1.0) { errHigh = x - xbest; break; }
  }
}

TGraphErrors * meanGraph(const float mean[], const float err[], const char * name) {
  TGraphErrors * g = new TGraphErrors(nPtBinsUsed);
  g->SetName(name);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    float lo = ana::ptBinsUsed[ipt], hi = ana::ptBinsUsed[ipt+1];
    g->SetPoint(ipt, (lo+hi)/2.0, mean[ipt]);
    g->SetPointError(ipt, (hi-lo)/2.0, err[ipt]);
  }
  return g;
}

TGraphErrors * ratioGraph(const float meanNum[], const float errNum[],
    const float meanDen[], const float errDen[], const char * name) {
  TGraphErrors * g = new TGraphErrors(nPtBinsUsed);
  g->SetName(name);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    float lo = ana::ptBinsUsed[ipt], hi = ana::ptBinsUsed[ipt+1];
    float num = meanNum[ipt], den = meanDen[ipt];
    if (num <= 0 || den <= 0) { g->SetPoint(ipt, (lo+hi)/2.0, 0); g->SetPointError(ipt, (hi-lo)/2.0, 0); continue; }
    float ratio = num/den;
    float err = ratio*sqrt(pow(errNum[ipt]/num,2) + pow(errDen[ipt]/den,2));
    g->SetPoint(ipt, (lo+hi)/2.0, ratio);
    g->SetPointError(ipt, (hi-lo)/2.0, err);
  }
  return g;
}

// Same sPHENIX label block as grid_insitu.C's drawSPhenixLabel - see that file's comment.
void drawSPhenixLabel(vector<string> samples, vector<string> features, float drawx, float drawy, int fontsize, float csize) {
  float titlescale = 1.25;
  float subtitlescale = 1.25;
  float ydiff = fontsize * 0.0017 * 700.0/csize;
  auto drawOne = [&](const char * text, float xp, float yp, int size) {
    TLatex * tex = new TLatex(xp, yp, text);
    tex->SetTextFont(43);
    tex->SetTextSize(size);
    tex->SetTextColor(kBlack);
    tex->SetLineWidth(1);
    tex->SetNDC();
    tex->Draw();
  };
  drawOne("#bf{#it{sPHENIX}} #kern[0.5]{Internal}", drawx, drawy, (int)(fontsize*titlescale));
  for (unsigned i = 0; i < samples.size(); i++) {
    drawOne(samples[i].c_str(), drawx, drawy-ydiff*subtitlescale*(i+1), (int)(fontsize*subtitlescale));
  }
  for (unsigned i = 0; i < features.size(); i++) {
    drawOne(features[i].c_str(), drawx, drawy-ydiff*subtitlescale*samples.size()-ydiff*(i+1)*subtitlescale, fontsize);
  }
}

// One comparison page: top panel is unfolded mean(x_J) vs pT for Truth (fixed) and raw
// (pa=1) unfolded Data; bottom panel is the raw ratio and the corrected ratio (evaluated
// at the scan's best-fit pa) - the latter should sit flat at 1 by construction. Same
// layout as grid_insitu.C's drawJESPage, with "Data/MC" -> "Unfolded/Truth".
void drawJESPage(TCanvas * c, const char * pdfPath, const char * label,
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
  drawSPhenixLabel({label}, {
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
void drawXjPage(TCanvas * c, const char * pdfPath, const char * label, float ptlow, float pthigh,
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

  drawSPhenixLabel({label, Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ptlow, pthigh)}, {
      Form("Jet R=%.1f", ana::JetRs[ir]),
      Form("p_{T}^{jet} > %.0f GeV", ana::jet_calib_pt_cut[ir]),
      Form("|#eta^{#gamma}|<%.1f, |#eta^{jet}|<%.1f", ana::etacut, ana::etacut-ana::JetRs[ir]),
      Form("#Delta#phi>%.0f#pi/%.0f", ana::oppnum, ana::oppden)
    }, .18, .85, 16, gPad->GetWh());

  c->SaveAs(pdfPath);
}

void grid_insitu_unfolded(string systag = "nominal", int na = 1000) {
  // Newly created histograms are not registered to any TDirectory, so the ~20
  // temporaries allocated per pa grid point (computeUnfoldedMeans, called na times) don't
  // pile up in gROOT's object list or collide on name across iterations - they're freed
  // by their own explicit `delete` calls instead.
  TH1::AddDirectory(kFALSE);

  const char * dataFile = Form("%s/Data_%s_insitu.root", insitu_dir, systag.c_str());
  vector<DataEvent> dataA = cacheDataEvents(dataFile, 0);
  vector<DataEvent> dataC = cacheDataEvents(dataFile, 2);
  cout << "Cached Data events: region A=" << dataA.size() << " region C=" << dataC.size() << endl;

  // Purity per ana::ptBins bin (all nPtBins, not just the used ones - see
  // cacheDataEvents comment) - computed once, held fixed across the whole pa scan, same
  // as grid_insitu.C.
  float purity[ana::nPtBins];
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    purity[ipt] = ana::getPurity(ana::ptBins[ipt], ana::ptBins[ipt+1], systag);
    cout << "Purity pt bin " << ipt << " [" << ana::ptBins[ipt] << "," << ana::ptBins[ipt+1] << "): " << purity[ipt] << endl;
  }

  // Response matrix: the physically meaningful photon+jet MC response (cross-section-
  // weighted combination of Photon5/10/20, type=1/isample=-1), built at the nominal
  // (uncorrected) MC jet energy scale and held fixed across the whole pa scan - only the
  // trial-pa-rescaled Data measured spectrum changes below. Same source/convention as
  // drawing/draw_purity_corrected.C's response matrix.
  drawer d("pythia", systag);
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
  // Grid scan: single overall jet-energy-scale factor pa, no pT-dependence - same
  // [0.95,1.05] range as grid_insitu.C. At each pa, purity-correct region A/C (all
  // ana::nPtBins slices), unfold through the fixed response above, and chi2 the
  // unfolded mean(x_J) (used bins only) against the fixed truth mean.
  // -----------------------------
  const float lowa = 0.95, higha = 1.05;

  TGraph * gchisqUnfold = new TGraph(na);
  gchisqUnfold->SetName("gchisq_unfolded");
  gchisqUnfold->SetTitle(";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});#chi^{2}");

  float minchisqUnfold = FLT_MAX, minpaUnfold = 1;
  int ibestUnfold = 0;

  for (int ia = 0; ia < na; ia++) {
    float pa = lowa + ia*(higha-lowa)/na;

    float mean[nPtBinsUsed], err[nPtBinsUsed];
    computeUnfoldedMeans(dataA, dataC, pa, purity, response, respRecoTemplate, mean, err);

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
  findError(gchisqUnfold, ibestUnfold, minchisqUnfold, errLowUnfold, errHighUnfold);

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
  computeUnfoldedMeans(dataA, dataC, 1.0,         purity, response, respRecoTemplate, rawMean,  rawErr,  &hUnfoldRaw);
  computeUnfoldedMeans(dataA, dataC, minpaUnfold, purity, response, respRecoTemplate, bestMean, bestErr, &hUnfoldBest);

  vector<TH1D*> hTruthPt(nPtBinsUsed);
  for (int k = 0; k < nPtBinsUsed; k++) {
    int ipt = ana::firstUsedPtBin + k;
    hTruthPt[k] = unfold_utility::unflattenXj(respTruthTemplate, ipt, Form("hxjtruth_pt%d", k));
  }

  TGraphErrors * gTruth      = meanGraph(truthMean, truthMeanErr, "gMeanTruth");
  TGraphErrors * gUnfoldRaw  = meanGraph(rawMean,  rawErr,  "gMeanUnfold_raw");
  TGraphErrors * gRatioRaw   = ratioGraph(rawMean,  rawErr,  truthMean, truthMeanErr, "gRatio_unfold_raw");
  TGraphErrors * gRatioCorr  = ratioGraph(bestMean, bestErr, truthMean, truthMeanErr, "gRatio_unfold_corrected");

  const char * pdfPath = Form("%s/grid_insitu_unfolded_%s.pdf", insitu_dir, systag.c_str());
  TCanvas * c = new TCanvas("c","",700,700);
  c->SaveAs(Form("%s[", pdfPath));
  drawJESPage(c, pdfPath, "Purity-corrected, unfolded", gTruth, gUnfoldRaw, gRatioRaw, gRatioCorr,
      minpaUnfold, errLowUnfold, errHighUnfold);
  for (int k = 0; k < nPtBinsUsed; k++) {
    drawXjPage(c, pdfPath, "Purity-corrected, unfolded", ana::ptBinsUsed[k], ana::ptBinsUsed[k+1],
        hTruthPt[k], hUnfoldRaw[k], hUnfoldBest[k]);
  }
  c->SaveAs(Form("%s]", pdfPath));
  cout << "Wrote " << pdfPath << endl;

  // -----------------------------
  // Save
  // -----------------------------
  const char * outfilename = Form("%s/grid_insitu_unfolded_%s.root", insitu_dir, systag.c_str());
  TFile * fout = TFile::Open(outfilename, "RECREATE");
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
  fout->Close();
  cout << "Wrote " << outfilename << endl;
}
