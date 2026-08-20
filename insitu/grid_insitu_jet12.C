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

// ana::findPtBin/etc. live in ana.cc, compiled into libgammajet_unfold.so - load it
// explicitly (see grid_insitu.C) so cling resolves the real compiled definitions.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Non-purity-corrected in-situ JES cross-check: same grid-scan machinery as
// grid_insitu.C, but comparing the *un-purity-corrected* Data Region A (no A-C
// background subtraction) against the *un-purity-corrected* Region A of a single
// QCD-dijet-triggered Pythia8 MC sample (Jet12_full - the 10x-higher-statistics
// towerntup production kept in the purity_check project; see src/treeuser.h's
// trigger=="Jet12_full" special case), instead of grid_insitu.C's purity-corrected
// Data vs. real prompt-photon Pythia8 gamma+jet MC (Photon5+10+20).
//
// This is a consistency check on the primary (purity-corrected, Photon-MC-referenced)
// result: does Data's naive, background-contaminated Region A line up with a QCD MC
// sample's own naive Region A (itself mostly fake-photon-triggered dijet background)
// under the same single-scale-factor fit? Neither side is ABCD-subtracted here, so
// there is no purity[] and no Region C in this macro at all - see draw_insitu_xj_jet12.C
// for the corresponding non-purity-corrected xJ-shape comparison.
//
// Reuses grid_insitu.C's Delta-chi2=1 error convention, ana::ptBinsUsed restriction
// (only the 15-20, 20-25, 25-35 GeV reported bins - see ana.h's ptBins/ptBinsUsed/
// firstUsedPtBin comment; drops both the 13-15 GeV migration buffer and the 35-100 GeV
// overflow, the latter for low Data statistics), and mean(x_J)-based chi2 definition
// unchanged - no new uncertainty-combination or bin-selection logic is introduced here
// (see gammajet_unfold/CLAUDE.md).

const char * insitu_dir = "/home/samson72/sphnx/gammajet_unfold/insitu";

// Same restriction as grid_insitu.C - only ana::ptBinsUsed (15-20, 20-25, 25-35 GeV),
// dropping both the 13-15 GeV migration buffer bin and the 35-100 GeV overflow bin
// (the latter for low Data statistics).
const int nPtBinsUsed = ana::nPtBinsUsed;

// Nominal cross-section weight for the Jet12 sample - same number as drawer.h's
// scalemap[isphoton=0][12] for sim="pythia". Jet12_full is a higher-statistics copy
// of the same underlying production (see src/treeuser.h), not a separate cross
// section, so this weight is a documentation/consistency convention only: with a
// single MC sample as the reference, it cancels out of every mean(x_J) computed
// below since all Jet12_full events share it.
map<int,double> jet_scale = {{12,3.997e+06}};

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
// (filename,weight) MC samples - same as grid_insitu.C's referenceMeans(), reused
// here with a single-entry sample list (Jet12 pythia, abcd=0/"Region A"). This is
// the fixed reference the grid scan compares against; these events are never
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

// Same Delta-chi2=1 (68% CL, one parameter) scan as grid_insitu.C's findError().
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
// factor pa - same as grid_insitu.C's computeRegionAMeans().
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

// sPHENIX label block - same convention as grid_insitu.C's drawSPhenixLabel(), see
// there for why this is reimplemented locally rather than going through drawer::drawAll().
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
// same layout as grid_insitu.C's drawJESPage(), single study (no purity-corrected page).
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
  jestext.DrawLatex(.18,.28, Form("Data to MC (Jet12) JES = %.4f #pm %.4f", pa, paErr));

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
// Jet12 pythia reference - same samples/weights as referenceMeans() above.
vector<TH1D*> buildMCXjByPtBin(const vector<pair<string,double>> & samples, int abcdSelect, const char * prefix) {
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
      if (abcd != abcdSelect) continue;
      int ipt = ana::findPtBin(pho_pt);
      if (ipt < ana::firstUsedPtBin || ipt >= ana::firstUsedPtBin + nPtBinsUsed) continue;
      ipt -= ana::firstUsedPtBin;
      h[ipt]->Fill(jet_pt/pho_pt, s.second);
    }
    f->Close();
  }
  return h;
}

// One xJ-distribution comparison page, for a single photon-pT bin: the fixed Jet12 MC
// reference, the raw (uncorrected) Data distribution, and the Data distribution at the
// study's best-fit jet-energy-scale - all shape-normalized and divided by bin width for
// display (ana::unfoldXjBins is non-uniform), same densityForDisplay convention as
// draw_insitu_xj_jet12.C/grid_insitu.C.
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

  drawSPhenixLabel({label, Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ptlow, pthigh)}, {
      Form("Jet R=%.1f", ana::JetRs[ir]),
      Form("p_{T}^{jet} > %.0f GeV", ana::jet_calib_pt_cut[ir]),
      Form("|#eta^{#gamma}|<%.1f, |#eta^{jet}|<%.1f", ana::etacut, ana::etacut-ana::JetRs[ir]),
      Form("#Delta#phi>%.0f#pi/%.0f", ana::oppnum, ana::oppden)
    }, .18, .85, 16, gPad->GetWh());

  c->SaveAs(pdfPath);
}

void grid_insitu_jet12(string systag = "nominal", int ir = 2) {

  const char * dataFile = Form("%s/Data_%s_insitu.root", insitu_dir, systag.c_str());
  vector<DataEvent> dataA = cacheDataEvents(dataFile, 0);
  cout << "Cached Data events: region A=" << dataA.size() << endl;

  vector<pair<string,double>> mcSamples = {
    {Form("%s/Jet12_full_pythia_%s_insitu.root", insitu_dir, systag.c_str()), jet_scale[12]},
  };

  float refMean[nPtBinsUsed], refMeanErr[nPtBinsUsed];
  referenceMeans(mcSamples, 0, refMean, refMeanErr);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    cout << "Jet12 MC reference <x_J> pt bin " << ipt << " [" << ana::ptBinsUsed[ipt] << "," << ana::ptBinsUsed[ipt+1]
         << "): " << refMean[ipt] << " +/- " << refMeanErr[ipt] << endl;
  }

  // -----------------------------
  // Grid scan: single overall jet-energy-scale factor pa, no pT-dependence - same
  // scan range/step as grid_insitu.C's gammajet-only mode.
  // -----------------------------
  const int na = 1000;
  const float lowa = 0.95, higha = 1.05;

  TGraph * gchisq = new TGraph(na);
  gchisq->SetName("gchisq_regionA_jet12ref");
  gchisq->SetTitle(";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});#chi^{2}");

  float minchisq = FLT_MAX, minpa = 1;
  int ibest = 0;

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

    float chisq = 0;
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

    gchisq->SetPoint(ia, pa, chisq);
    if (chisq < minchisq) { minchisq = chisq; minpa = pa; ibest = ia; }
  }

  float errLow, errHigh;
  findError(gchisq, ibest, minchisq, errLow, errHigh);

  cout << "\nFINAL RESULT (jet R=" << ana::JetRs[ir] << ", systag=" << systag << ", non-purity-corrected, Jet12 reference)\n";
  cout << "Region A (Data) vs. Region A (Jet12 MC):  p_a = " << minpa
       << " +" << errHigh << "/-" << errLow << " (chi2=" << minchisq << ")" << endl;

  // -----------------------------
  // Build x_J histograms per photon-pT bin: the fixed Jet12 MC reference, Data at
  // pa=1 (raw) and at the scan's best-fit pa (corrected) - same non-uniform binning
  // as grid_insitu.C.
  // -----------------------------
  vector<TH1D*> hxjMC_pt    = buildMCXjByPtBin(mcSamples, 0, "hxjA_jet12");
  vector<TH1D*> hxjA_raw_pt  = buildXjByPtBin(dataA, 1.0,   "hxjA_data_raw");
  vector<TH1D*> hxjA_corr_pt = buildXjByPtBin(dataA, minpa, "hxjA_data_corr");

  auto sumPtBins = [&](const vector<TH1D*> & h, const char * name) {
    TH1D * hsum = (TH1D*)h[0]->Clone(name);
    for (int ipt = 1; ipt < nPtBinsUsed; ipt++) hsum->Add(h[ipt]);
    return hsum;
  };
  TH1D * hxjA_jet12          = sumPtBins(hxjMC_pt,     "hxjA_jet12");
  TH1D * hxjA_data_raw       = sumPtBins(hxjA_raw_pt,  "hxjA_data_raw");
  TH1D * hxjA_data_bestscale = sumPtBins(hxjA_corr_pt, "hxjA_data_bestscale");

  // -----------------------------
  // Mean(x_J) vs pT comparison plot - top: MC vs raw Data; bottom: raw ratio vs
  // corrected ratio (evaluated at the scan's best-fit pa). Single page - no
  // purity-corrected branch, unlike grid_insitu.C.
  // -----------------------------
  gStyle->SetOptStat(0);

  float rawMeanA[nPtBinsUsed], rawErrA[nPtBinsUsed];
  float corrMeanA[nPtBinsUsed], corrErrA[nPtBinsUsed];
  computeRegionAMeans(dataA, 1.0,   rawMeanA,  rawErrA);
  computeRegionAMeans(dataA, minpa, corrMeanA, corrErrA);

  TGraphErrors * gMC          = meanGraph(refMean, refMeanErr, "gMeanMC_jet12");
  TGraphErrors * gDataRaw     = meanGraph(rawMeanA,  rawErrA,  "gMeanData_regionA_raw");
  TGraphErrors * gRatioRaw    = ratioGraph(rawMeanA,  rawErrA,  refMean, refMeanErr, "gRatio_regionA_raw");
  TGraphErrors * gRatioCorr   = ratioGraph(corrMeanA, corrErrA, refMean, refMeanErr, "gRatio_regionA_corrected");

  const char * pdfPath = Form("%s/grid_insitu_jet12_%s.pdf", insitu_dir, systag.c_str());
  TCanvas * c = new TCanvas("c","",700,700);
  c->SaveAs(Form("%s[", pdfPath));
  drawJESPage(c, pdfPath, "Region A, non-purity-corrected (Jet12 ref.)", ir, gMC, gDataRaw, gRatioRaw, gRatioCorr, minpa, errLow, errHigh);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    drawXjPage(c, pdfPath, "Region A, non-purity-corrected (Jet12 ref.)", ir, ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1],
        hxjMC_pt[ipt], hxjA_raw_pt[ipt], hxjA_corr_pt[ipt]);
  }
  c->SaveAs(Form("%s]", pdfPath));
  cout << "Wrote " << pdfPath << endl;

  // -----------------------------
  // Save
  // -----------------------------
  const char * outfilename = Form("%s/grid_insitu_jet12_%s.root", insitu_dir, systag.c_str());
  TFile * fout = TFile::Open(outfilename, "RECREATE");
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

  TTree * wt = new TTree("results", "best-fit jet energy scale result (non-purity-corrected, Jet12 reference)");
  float wpa = minpa, wchisq = minchisq, werrLow = errLow, werrHigh = errHigh;
  wt->Branch("pa_regionA_jet12ref", &wpa);
  wt->Branch("chisq_regionA_jet12ref", &wchisq);
  wt->Branch("errLow_regionA_jet12ref", &werrLow);
  wt->Branch("errHigh_regionA_jet12ref", &werrHigh);
  wt->Fill();
  wt->Write();
  fout->Close();
  cout << "Wrote " << outfilename << endl;
}
