#include "../src/ana.h"
#include "../src/insitu_utility.h"
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TStyle.h"
using namespace std;

R__LOAD_LIBRARY(libgammajet_unfold.so);

// Digs into the sawtooth pattern grid_insitu.C (shape method)'s Region-A shape-chi2 curve
// picked up after the low-xJ floor (insitu_utility::lowXjFloor) was added: as pa scans,
// x_J(pa) = jet_pt/pa/pho_pt shifts continuously, so each event crosses the FIXED
// lowXj[ptbin] threshold at its own pa_cross = jet_pt/(lowXj[ptbin]*pho_pt), leaving the
// sample permanently as pa increases past that point (x_J is monotonically decreasing in
// pa, so this is a one-way exit, not a flicker). This script recomputes NA(pa) per pT bin
// (the shared normalization denominator for every bin fraction) alongside chi2A(pa) over
// the same scan grid as production, to see whether chi2 jumps line up with clusters of
// events exiting together rather than one at a time.
// insitu/ is split into inputs/ (the raw Data/Photon insitu ntuples, written by
// unfolder.h's production pipeline) and pdfs/ (this debug macro's own .pdf output,
// no .root output).
const char * insitu_input_dir = ana::path("insitu/inputs");
const char * insitu_pdf_dir   = ana::path("insitu/pdfs");
const int nPtBinsUsed = ana::nPtBinsUsed;
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;

map<int,double> photon_scale = {{5,146359.3},{10,6944.675},{20,130.4461}};

void referenceShape(const vector<pair<string,double>> & samples, int abcdSelect, int ir,
    vector<vector<double>> & refFrac, vector<vector<double>> & refFracErr, const float lowXj[]) {
  vector<vector<double>> sumw(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
  vector<vector<double>> sumw2(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
  for (auto & s : samples) {
    TFile * f = TFile::Open(s.first.c_str(), "READ");
    if (!f || f->IsZombie()) { cout << "WARNING: could not open " << s.first << endl; continue; }
    TTree * t = (TTree*)f->Get("insitutree");
    Float_t pho_pt, jet_pt, mcWeight;
    Int_t abcd, evIr;
    t->SetBranchAddress("pho_pt", &pho_pt);
    t->SetBranchAddress("jet_pt", &jet_pt);
    t->SetBranchAddress("abcd", &abcd);
    t->SetBranchAddress("weight", &mcWeight);
    t->SetBranchAddress("ir", &evIr);
    Long64_t nentries = t->GetEntries();
    for (Long64_t e = 0; e < nentries; e++) {
      t->GetEntry(e);
      if (abcd != abcdSelect) continue;
      if (evIr != ir) continue;
      int ipt = ana::findPtBin(pho_pt);
      if (ipt < ana::firstUsedPtBin || ipt >= ana::firstUsedPtBin + nPtBinsUsed) continue;
      ipt -= ana::firstUsedPtBin;
      double x = jet_pt/pho_pt;
      if (x < lowXj[ipt]) continue;
      int ixj = ana::findUnfoldXjBin(x);
      if (ixj < 0 || ixj >= ana::nUnfoldXjBins) continue;
      double w = s.second*mcWeight;
      sumw[ipt][ixj]  += w;
      sumw2[ipt][ixj] += w*w;
    }
    f->Close();
  }
  refFrac.assign(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
  refFracErr.assign(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    double N = 0;
    for (double w : sumw[ipt]) N += w;
    if (N <= 0) continue;
    for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
      refFrac[ipt][ixj]    = sumw[ipt][ixj]/N;
      refFracErr[ipt][ixj] = sqrt(sumw2[ipt][ixj])/N;
    }
  }
}

void debug_shapechi2_sawtooth(string systag = "nominal", int ir = 2) {
  gStyle->SetOptStat(0);

  string dataFile = insitu_utility::insituFilename(insitu_input_dir, "Data", "", systag);
  vector<DataEvent> dataA = insitu_utility::cacheDataEvents(dataFile.c_str(), 0, ir);
  cout << "Cached Data Region A events: " << dataA.size() << endl;

  float lowXj[nPtBinsUsed];
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    lowXj[ipt] = insitu_utility::lowXjFloor(ir, ana::ptBinsUsed[ipt]);
    cout << "lowXj[" << ipt << "] (pT " << ana::ptBinsUsed[ipt] << "-" << ana::ptBinsUsed[ipt+1]
         << " GeV) = " << lowXj[ipt] << endl;
  }

  // pa_cross per event: the pa above which this event's x_J(pa) drops below its pT
  // bin's floor and it exits the sample. Events are sorted here purely so the printed
  // "cluster" report below reads in scan order - the scan loop itself doesn't need this.
  struct Cross { float pa_cross; int ipt; };
  vector<Cross> crossings;
  for (auto & ev : dataA) {
    float pa_cross = ev.jet_pt/(lowXj[ev.ptbin]*ev.pho_pt);
    if (pa_cross >= insitu_utility::scanLow && pa_cross < insitu_utility::scanHigh) {
      crossings.push_back({pa_cross, ev.ptbin});
    }
  }
  sort(crossings.begin(), crossings.end(), [](const Cross & a, const Cross & b) { return a.pa_cross < b.pa_cross; });
  cout << "\nEvents with a floor-crossing pa inside the scan window: " << crossings.size()
       << " / " << dataA.size() << " total Region A events" << endl;

  vector<pair<string,double>> mcSamples = {
    {insitu_utility::insituFilename(insitu_input_dir, "Photon5",  "pythia", systag), photon_scale[5]},
    {insitu_utility::insituFilename(insitu_input_dir, "Photon10", "pythia", systag), photon_scale[10]},
    {insitu_utility::insituFilename(insitu_input_dir, "Photon20", "pythia", systag), photon_scale[20]},
  };
  vector<vector<double>> refFrac, refFracErr;
  referenceShape(mcSamples, 0, ir, refFrac, refFracErr, lowXj);

  const int na = insitu_utility::scanN;
  const float lowa = insitu_utility::scanLow, higha = insitu_utility::scanHigh;

  TGraph * gchisq = new TGraph(na);
  gchisq->SetName("gchisqA_debug");
  vector<TGraph*> gNA(nPtBinsUsed);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    gNA[ipt] = new TGraph(na);
    gNA[ipt]->SetName(Form("gNA_pt%d", ipt));
  }

  vector<float> chisqArr(na);
  for (int ia = 0; ia < na; ia++) {
    float pa = lowa + ia*(higha-lowa)/na;

    vector<vector<double>> countA(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
    for (auto & ev : dataA) {
      float x = (ev.jet_pt/pa)/ev.pho_pt;
      if (x < lowXj[ev.ptbin]) continue;
      int ixj = ana::findUnfoldXjBin(x);
      if (ixj < 0 || ixj >= ana::nUnfoldXjBins) continue;
      countA[ev.ptbin][ixj] += 1;
    }

    float chisqA = 0;
    for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
      double NA = 0;
      for (double n : countA[ipt]) NA += n;
      gNA[ipt]->SetPoint(ia, pa, NA);
      if (NA <= 0) continue;
      double errFloor = 1.0/NA;
      for (int ixj = 0; ixj < nXjBinsForChi2; ixj++) {
        double fData = countA[ipt][ixj]/NA;
        double errData = sqrt(countA[ipt][ixj])/NA;
        double errt = sqrt(errData*errData + refFracErr[ipt][ixj]*refFracErr[ipt][ixj]);
        errt = std::max(errt, errFloor);
        double diff = fData - refFrac[ipt][ixj];
        chisqA += diff*diff/(errt*errt);
      }
    }
    gchisq->SetPoint(ia, pa, chisqA);
    chisqArr[ia] = chisqA;
  }

  // Find the biggest upward jumps in chi2 between consecutive scan points, and report
  // how many events cross the floor in that same pa step, in which pT bin(s).
  vector<pair<float,int>> jumps; // (delta chi2, ia)
  for (int ia = 1; ia < na; ia++) jumps.push_back({chisqArr[ia] - chisqArr[ia-1], ia});
  sort(jumps.begin(), jumps.end(), [](const pair<float,int> & a, const pair<float,int> & b) { return a.first > b.first; });

  cout << "\nTop 10 upward chi2 jumps (Region A):" << endl;
  for (int k = 0; k < 10 && k < (int)jumps.size(); k++) {
    int ia = jumps[k].second;
    float pa_lo = lowa + (ia-1)*(higha-lowa)/na;
    float pa_hi = lowa + ia*(higha-lowa)/na;
    int nCross = 0;
    map<int,int> crossByPt;
    for (auto & c : crossings) {
      if (c.pa_cross >= pa_lo && c.pa_cross < pa_hi) { nCross++; crossByPt[c.ipt]++; }
    }
    cout << "  pa=" << pa_lo << "->" << pa_hi << "  dchi2=" << jumps[k].first
         << "  chi2=" << chisqArr[ia-1] << "->" << chisqArr[ia]
         << "  events crossing floor in this step=" << nCross;
    for (auto & pr : crossByPt) cout << " [ptbin" << pr.first << ":" << pr.second << "]";
    cout << endl;
  }

  // Also report the widest gaps between consecutive sorted crossing pa's, per pT bin -
  // a large gap means a long pa stretch with a completely fixed NA (no exits), which
  // would show up as one of the smoother "tooth" segments between cliffs.
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    vector<float> pas;
    for (auto & c : crossings) if (c.ipt == ipt) pas.push_back(c.pa_cross);
    cout << "\nptbin " << ipt << ": " << pas.size() << " floor-crossings in-window";
    if (!pas.empty()) cout << ", first=" << pas.front() << " last=" << pas.back();
    cout << endl;
  }

  TCanvas * c = new TCanvas("c", "", 1000, 900);
  c->Divide(1, 2);
  c->cd(1);
  gPad->SetLeftMargin(.12);
  gchisq->SetTitle(";p_{a};Shape #chi^{2} (Region A)");
  gchisq->SetMarkerStyle(kFullDotSmall);
  gchisq->Draw("AP");

  c->cd(2);
  gPad->SetLeftMargin(.12);
  int colors[3] = {kBlue+1, kRed+1, kGreen+2};
  TLegend * leg = new TLegend(.7, .7, .92, .9);
  leg->SetBorderSize(0);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    gNA[ipt]->SetLineColor(colors[ipt]);
    gNA[ipt]->SetMarkerColor(colors[ipt]);
    gNA[ipt]->SetMarkerStyle(kFullDotSmall);
    gNA[ipt]->SetTitle(";p_{a};N_{A}(p_{a}) surviving the floor");
    gNA[ipt]->Draw(ipt == 0 ? "AP" : "P SAME");
    leg->AddEntry(gNA[ipt], Form("pT bin %d [%.0f,%.0f)", ipt, ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1]), "p");
  }
  leg->Draw();

  const char * pdfPath = Form("%s/debug_shapechi2_sawtooth_%s_%s.pdf", insitu_pdf_dir, systag.c_str(), ana::rnames[ir]);
  c->SaveAs(pdfPath);
  cout << "\nWrote " << pdfPath << endl;
}
