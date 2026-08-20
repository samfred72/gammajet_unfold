#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
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

const char * insitu_dir = "/home/samson72/sphnx/gammajet_unfold/insitu";

// Only ana::ptBinsUsed (15-20, 20-25, 25-35 GeV) is used for every calculation and
// plot below - both the low-pT migration-only buffer bin (13-15 GeV, ana::ptBins[0])
// and the high-pT overflow bin (35-100 GeV) are excluded, since neither is a reported
// physics bin (see ana.h's ptBins/ptBinsUsed/firstUsedPtBin comment) and the top one
// also has too few Data events for a meaningful in-situ point.
const int nPtBinsUsed = ana::nPtBinsUsed;

// Cross-section weights for combining the Photon5/10/20 MC samples - same numbers as
// drawer.h's scalemap[isphoton=1][sample] for sim="pythia".
map<int,double> photon_scale = {{5,146359.3},{10,6944.675},{20,130.4461}};

struct DataEvent { float pho_pt, jet_pt; int ptbin; };

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
    if (ipt < ana::firstUsedPtBin || ipt >= ana::firstUsedPtBin + nPtBinsUsed) continue;
    ipt -= ana::firstUsedPtBin;
    events.push_back({pho_pt, jet_pt, ipt});
  }
  f->Close();
  return events;
}

// Weighted mean/error of x=jet_pt/pho_pt per photon-pT bin, combining several
// (filename,weight) MC samples. This is the fixed reference the grid scan compares
// against - MC's own jet energy scale is the "truth" here, so these events are never
// rescaled below.
void referenceMeans(const vector<pair<string,double>> & samples, int abcdSelect,
    float refMean[], float refMeanErr[]) {
  vector<double> sumw(nPtBinsUsed,0), sumw2(nPtBinsUsed,0), sumwx(nPtBinsUsed,0), sumwx2(nPtBinsUsed,0);
  for (auto & s : samples) {
    TFile * f = TFile::Open(s.first.c_str(), "READ");
    if (!f || f->IsZombie()) {
      cout << "WARNING: could not open " << s.first << endl;
      continue;
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
      if (ipt < ana::firstUsedPtBin || ipt >= ana::firstUsedPtBin + nPtBinsUsed) continue;
      ipt -= ana::firstUsedPtBin;
      double x = jet_pt/pho_pt;
      double w = s.second;
      sumw[ipt]   += w;
      sumw2[ipt]  += w*w;
      sumwx[ipt]  += w*x;
      sumwx2[ipt] += w*x*x;
    }
    f->Close();
  }
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    if (sumw[ipt] <= 0) { refMean[ipt] = 0; refMeanErr[ipt] = 0; continue; }
    double mean = sumwx[ipt]/sumw[ipt];
    double var  = sumwx2[ipt]/sumw[ipt] - mean*mean;
    double neff = sumw[ipt]*sumw[ipt]/sumw2[ipt]; // Kish effective sample size
    refMean[ipt] = mean;
    refMeanErr[ipt] = sqrt(std::max(var,0.)/neff);
  }
}

// Scans outward from the minimum on a chi2-vs-pa graph for the two points where
// chi2 first crosses minchisq+1 (68% CL for one parameter) - the same Delta-chi2=1
// convention grid_insitu.C uses for its own (!domulti) single-scale-factor case.
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

// Region-A-only mean(x_J)/error per photon-pT bin, at a given trial jet-energy-scale
// factor pa - same computation as the inner loop of the grid scan above, factored out
// so it can also be evaluated at pa=1 (raw, uncorrected) and at the scan's winning pa
// for the mean-vs-pT comparison plot below.
void computeRegionAMeans(const vector<DataEvent> & dataA, float pa, float mean[], float err[]) {
  vector<double> sum(nPtBinsUsed,0), sum2(nPtBinsUsed,0);
  vector<int> count(nPtBinsUsed,0);
  for (auto & ev : dataA) {
    float x = (ev.jet_pt/pa)/ev.pho_pt;
    sum[ev.ptbin]  += x;
    sum2[ev.ptbin] += x*x;
    count[ev.ptbin]++;
  }
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    if (count[ipt] == 0) { mean[ipt] = 0; err[ipt] = 0; continue; }
    double m   = sum[ipt]/count[ipt];
    double var = sum2[ipt]/count[ipt] - m*m;
    mean[ipt] = m;
    err[ipt]  = sqrt(std::max(var,0.)/count[ipt]);
  }
}

// Purity-corrected mean(x_J)/error per photon-pT bin, at a given trial pa - both region
// A and region C are scaled by the same pa, purity[] is fixed (see file header).
void computeCorrectedMeans(const vector<DataEvent> & dataA, const vector<DataEvent> & dataC,
    float pa, const float purity[], float mean[], float err[]) {
  vector<double> sumA(nPtBinsUsed,0), sumA2(nPtBinsUsed,0);
  vector<int> countA(nPtBinsUsed,0);
  vector<double> sumC(nPtBinsUsed,0), sumC2(nPtBinsUsed,0);
  vector<int> countC(nPtBinsUsed,0);
  for (auto & ev : dataA) {
    float x = (ev.jet_pt/pa)/ev.pho_pt;
    sumA[ev.ptbin] += x; sumA2[ev.ptbin] += x*x; countA[ev.ptbin]++;
  }
  for (auto & ev : dataC) {
    float x = (ev.jet_pt/pa)/ev.pho_pt;
    sumC[ev.ptbin] += x; sumC2[ev.ptbin] += x*x; countC[ev.ptbin]++;
  }
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    if (countA[ipt] == 0 || countC[ipt] == 0) { mean[ipt] = 0; err[ipt] = 0; continue; }
    double NA = countA[ipt], NC = countC[ipt];
    double scale = (1-purity[ipt])*(NA/NC);
    double sumXcorr  = sumA[ipt]  - scale*sumC[ipt];
    double sumX2corr = sumA2[ipt] - scale*sumC2[ipt];
    double Ncorr = NA - scale*NC;
    if (Ncorr <= 0) { mean[ipt] = 0; err[ipt] = 0; continue; }
    double m   = sumXcorr/Ncorr;
    double var = sumX2corr/Ncorr - m*m;
    mean[ipt] = m;
    err[ipt]  = sqrt(std::max(var,0.)/Ncorr);
  }
}

// Mean(x_J) vs. photon pT, one point per ana::ptBins bin (x error = half bin width).
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

// Ratio of two mean(x_J) arrays (e.g. Data/MC) vs. photon pT, errors combined assuming
// the numerator and denominator are independent.
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

// sPHENIX label block: bold-italic "sPHENIX Internal" title, then one line per sample,
// then one line per feature - same text/font convention and stacking formula as
// drawer::drawAll() (see src/drawer.cc), reimplemented locally so this self-contained
// macro doesn't have to construct a full drawer (which opens a batch of unrelated
// unfolding-output files it has no other use for).
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

// x_J histogram per photon-pT bin (ana::unfoldXjBins binning) from a cached Data
// sample, at a given trial jet-energy-scale factor pa.
vector<TH1D*> buildXjByPtBin(const vector<DataEvent> & data, float pa, const char * prefix) {
  vector<TH1D*> h(nPtBinsUsed);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    h[ipt] = new TH1D(Form("%s_pt%d", prefix, ipt), ";x_{J#gamma};Counts", ana::nUnfoldXjBins, ana::unfoldXjBins);
  }
  for (auto & ev : data) {
    h[ev.ptbin]->Fill((ev.jet_pt/pa)/ev.pho_pt);
  }
  return h;
}

// x_J histogram per photon-pT bin for the fixed (never rescaled) cross-section-weighted
// Pythia8 MC reference - same samples/weights as referenceMeans() above.
vector<TH1D*> buildMCXjByPtBin(const vector<pair<string,double>> & samples, const char * prefix) {
  vector<TH1D*> h(nPtBinsUsed);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    h[ipt] = new TH1D(Form("%s_pt%d", prefix, ipt), ";x_{J#gamma};Counts", ana::nUnfoldXjBins, ana::unfoldXjBins);
  }
  for (auto & s : samples) {
    TFile * f = TFile::Open(s.first.c_str(), "READ");
    if (!f || f->IsZombie()) continue;
    TTree * t = (TTree*)f->Get("insitutree");
    Float_t pho_pt, jet_pt; Int_t abcd;
    t->SetBranchAddress("pho_pt", &pho_pt);
    t->SetBranchAddress("jet_pt", &jet_pt);
    t->SetBranchAddress("abcd", &abcd);
    Long64_t nentries = t->GetEntries();
    for (Long64_t e = 0; e < nentries; e++) {
      t->GetEntry(e);
      if (abcd != 0) continue;
      int ipt = ana::findPtBin(pho_pt);
      if (ipt < ana::firstUsedPtBin || ipt >= ana::firstUsedPtBin + nPtBinsUsed) continue;
      ipt -= ana::firstUsedPtBin;
      h[ipt]->Fill(jet_pt/pho_pt, s.second);
    }
    f->Close();
  }
  return h;
}

// Purity-correct region A/C per photon-pT bin: A - (1-P)*(N_A/N_C)*C, same formula as
// draw_insitu_xj.C's purityCorrectP(), applied bin-by-bin. hA/hC must already be built
// at the same jet-energy-scale factor (see buildXjByPtBin above).
vector<TH1D*> purityCorrectByPtBin(const vector<TH1D*> & hA, const vector<TH1D*> & hC,
    const float purity[], const char * prefix) {
  vector<TH1D*> h(nPtBinsUsed);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    h[ipt] = (TH1D*)hA[ipt]->Clone(Form("%s_pt%d", prefix, ipt));
    float NA = hA[ipt]->Integral(), NC = hC[ipt]->Integral();
    float scale = (NC > 0) ? (1-purity[ipt])*(NA/NC) : 0;
    for (int ib = 1; ib <= h[ipt]->GetNbinsX(); ib++) {
      float content = hA[ipt]->GetBinContent(ib) - scale*hC[ipt]->GetBinContent(ib);
      h[ipt]->SetBinContent(ib, content);
    }
  }
  return h;
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

  drawSPhenixLabel({label, Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ptlow, pthigh)}, {
      Form("Jet R=%.1f", ana::JetRs[ir]),
      Form("p_{T}^{jet} > %.0f GeV", ana::jet_calib_pt_cut[ir]),
      Form("|#eta^{#gamma}|<%.1f, |#eta^{jet}|<%.1f", ana::etacut, ana::etacut-ana::JetRs[ir]),
      Form("#Delta#phi>%.0f#pi/%.0f", ana::oppnum, ana::oppden)
    }, .18, .85, 16, gPad->GetWh());

  c->SaveAs(pdfPath);
}

void grid_insitu(string systag = "nominal", int ir = 2) {

  const char * dataFile = Form("%s/Data_%s_insitu.root", insitu_dir, systag.c_str());
  vector<DataEvent> dataA = cacheDataEvents(dataFile, 0);
  vector<DataEvent> dataC = cacheDataEvents(dataFile, 2);
  cout << "Cached Data events: region A=" << dataA.size() << " region C=" << dataC.size() << endl;

  float refMean[nPtBinsUsed], refMeanErr[nPtBinsUsed];
  referenceMeans({
      {Form("%s/Photon5_pythia_%s_insitu.root",  insitu_dir, systag.c_str()), photon_scale[5]},
      {Form("%s/Photon10_pythia_%s_insitu.root", insitu_dir, systag.c_str()), photon_scale[10]},
      {Form("%s/Photon20_pythia_%s_insitu.root", insitu_dir, systag.c_str()), photon_scale[20]},
    }, 0, refMean, refMeanErr);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    cout << "MC reference <x_J> pt bin " << ipt << " [" << ana::ptBinsUsed[ipt] << "," << ana::ptBinsUsed[ipt+1]
         << "): " << refMean[ipt] << " +/- " << refMeanErr[ipt] << endl;
  }

  // Purity per photon-pT bin - computed once, held fixed across the whole pa scan.
  float purity[nPtBinsUsed];
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    purity[ipt] = ana::getPurity(ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1], systag);
    cout << "Purity pt bin " << ipt << ": " << purity[ipt] << endl;
  }

  // -----------------------------
  // Grid scan: single overall jet-energy-scale factor pa, no pT-dependence - matches
  // run_grid.sh's gammajet-only mode (na=1000 over [0.95,1.05], nb=1, pb=0).
  // -----------------------------
  const int na = 1000;
  const float lowa = 0.95, higha = 1.05;

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
      sumA[ev.ptbin]  += x;
      sumA2[ev.ptbin] += x*x;
      countA[ev.ptbin]++;
    }
    vector<double> sumC(nPtBinsUsed,0), sumC2(nPtBinsUsed,0);
    vector<int> countC(nPtBinsUsed,0);
    for (auto & ev : dataC) {
      float x = (ev.jet_pt/pa)/ev.pho_pt;
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

      // Purity-corrected: A - (1-P)*(N_A/N_C)*C, both scaled by the same pa.
      if (countA[ipt] > 0 && countC[ipt] > 0) {
        double NA = countA[ipt], NC = countC[ipt];
        double scale = (1-purity[ipt])*(NA/NC);
        double sumXcorr  = sumA[ipt]  - scale*sumC[ipt];
        double sumX2corr = sumA2[ipt] - scale*sumC2[ipt];
        double Ncorr = NA - scale*NC;
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
  findError(gchisqA,    ibestA,    minchisqA,    errLowA,    errHighA);
  findError(gchisqCorr, ibestCorr, minchisqCorr, errLowCorr, errHighCorr);

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
    {Form("%s/Photon5_pythia_%s_insitu.root",  insitu_dir, systag.c_str()), photon_scale[5]},
    {Form("%s/Photon10_pythia_%s_insitu.root", insitu_dir, systag.c_str()), photon_scale[10]},
    {Form("%s/Photon20_pythia_%s_insitu.root", insitu_dir, systag.c_str()), photon_scale[20]},
  };
  vector<TH1D*> hxjMC_pt       = buildMCXjByPtBin(mcSamples, "hxjA_pythia");
  vector<TH1D*> hxjA_raw_pt    = buildXjByPtBin(dataA, 1.0,       "hxjA_data_raw");
  vector<TH1D*> hxjA_corr_pt   = buildXjByPtBin(dataA, minpaA,    "hxjA_data_corr");
  vector<TH1D*> hxjC_raw_pt    = buildXjByPtBin(dataC, 1.0,       "hxjC_data_raw");
  vector<TH1D*> hxjA_atCorr_pt = buildXjByPtBin(dataA, minpaCorr, "hxjA_data_atCorrScale");
  vector<TH1D*> hxjC_atCorr_pt = buildXjByPtBin(dataC, minpaCorr, "hxjC_data_atCorrScale");
  // Purity-corrected (A - (1-P)*(N_A/N_C)*C), once raw (pa=1) and once at the
  // purity-corrected study's best-fit pa - purity is fixed either way (file header).
  vector<TH1D*> hxjcorr_raw_pt  = purityCorrectByPtBin(hxjA_raw_pt,    hxjC_raw_pt,    purity, "hxjcorrected_data_raw");
  vector<TH1D*> hxjcorr_best_pt = purityCorrectByPtBin(hxjA_atCorr_pt, hxjC_atCorr_pt, purity, "hxjcorrected_data_bestscale");

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
  computeRegionAMeans(dataA, 1.0,    rawMeanA,  rawErrA);
  computeRegionAMeans(dataA, minpaA, corrMeanA, corrErrA);

  float rawMeanCorr[nPtBinsUsed], rawErrCorr[nPtBinsUsed];
  float bestMeanCorr[nPtBinsUsed], bestErrCorr[nPtBinsUsed];
  computeCorrectedMeans(dataA, dataC, 1.0,        purity, rawMeanCorr,  rawErrCorr);
  computeCorrectedMeans(dataA, dataC, minpaCorr,  purity, bestMeanCorr, bestErrCorr);

  TGraphErrors * gMC              = meanGraph(refMean, refMeanErr, "gMeanMC");
  TGraphErrors * gDataRaw_A       = meanGraph(rawMeanA,  rawErrA,  "gMeanData_regionA_raw");
  TGraphErrors * gRatioRaw_A      = ratioGraph(rawMeanA,  rawErrA,  refMean, refMeanErr, "gRatio_regionA_raw");
  TGraphErrors * gRatioCorr_A     = ratioGraph(corrMeanA, corrErrA, refMean, refMeanErr, "gRatio_regionA_corrected");

  TGraphErrors * gDataRaw_Corr    = meanGraph(rawMeanCorr,  rawErrCorr,  "gMeanData_puritycorrected_raw");
  TGraphErrors * gRatioRaw_Corr   = ratioGraph(rawMeanCorr,  rawErrCorr,  refMean, refMeanErr, "gRatio_puritycorrected_raw");
  TGraphErrors * gRatioCorr_Corr  = ratioGraph(bestMeanCorr, bestErrCorr, refMean, refMeanErr, "gRatio_puritycorrected_corrected");

  const char * pdfPath = Form("%s/grid_insitu_%s.pdf", insitu_dir, systag.c_str());
  TCanvas * c = new TCanvas("c","",700,700);
  c->SaveAs(Form("%s[", pdfPath));
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
  c->SaveAs(Form("%s]", pdfPath));
  cout << "Wrote " << pdfPath << endl;

  // -----------------------------
  // Save
  // -----------------------------
  const char * outfilename = Form("%s/grid_insitu_%s.root", insitu_dir, systag.c_str());
  TFile * fout = TFile::Open(outfilename, "RECREATE");
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
  fout->Close();
  cout << "Wrote " << outfilename << endl;
}
