#include "../src/ana.h"
#include "../src/insitu_utility.h"
#include <string>
#include <vector>
#include <cfloat>
#include <cmath>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TParameter.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
using namespace std;

// ana::/insitu_utility:: are implemented in ana.cc/insitu_utility.cc, compiled into
// libgammajet_unfold.so - load it explicitly (see insitu/grid_insitu.C's identical
// comment) so cling resolves the real compiled definitions.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// In-situ JES closure test, stage 2: reads the per-sample insitu-style trees
// make_closure_trees.C wrote (one "data"/"sim" pair per Photon5/10/20 pythia sample - a
// "data" half with jet_pt scaled by a known injected factor, and a "sim" half left at
// nominal reco scale), stitches the three samples together by cross-section weight
// (same weights make_closure_trees.C already baked into each row's "weight" branch),
// and runs the SAME Region-A mean(x_J) chi2 scan insitu/grid_insitu.C uses on real Data
// - scanning a single overall jet-energy-scale factor pa (jet_pt_corrected =
// jet_pt/pa) to find the value that makes the combined "data" half's mean(x_J) match
// the combined "sim" half's, per photon-pT bin. Since the "data" half's true best-fit
// pa is exactly the injected scale by construction, this is the closure check: does the
// fitting method recover what was put in?
//
// R=0.4 (ir=2) only - the primary/nominal jet radius used throughout this project (e.g.
// insitu/draw_insitu_xj.C's fixed ir=2, unfolder.cc's jesNominal[2] baked into
// production) - not the full 7-radius scan grid_insitu.C runs for the real in-situ
// study. Region A only, no purity/background-subtraction correction (unlike
// grid_insitu.C's second, purity-corrected method) - this closure tests the
// scale-fitting machinery itself, not the ABCD purity correction, so there's no need to
// reproduce that half of grid_insitu.C here.
//
// Deliberately does NOT reuse insitu_utility::cacheDataEvents/computeRegionAMeans/
// referenceMeans: cacheDataEvents/computeRegionAMeans are unweighted (correct for real
// Data, which always has weight=1), but BOTH halves here are stitched, cross-section-
// weighted MC (see make_closure_trees.C) - the "data" half needs the same weighted
// mean/error treatment referenceMeans gives the MC reference in grid_insitu.C, not the
// unweighted one grid_insitu.C's actual Region-A Data scan uses. The local
// cacheClosureEvents/weightedRegionAMeans below are that one weighted implementation,
// reused symmetrically for both halves, and cacheClosureEvents is called once per
// sample file below and concatenated - the actual stitch, mirroring how
// insitu_utility::referenceMeans loops a (filename,weight) list, just with the weight
// already baked into each row instead of passed alongside the filename.

const char * closure_input_dir  = ana::path("insitu_closure/inputs");
const char * closure_output_dir = ana::path("insitu_closure/output");
const char * closure_pdf_dir    = ana::path("insitu_closure/pdfs");

const int nPtBinsUsed = ana::nPtBinsUsed;

// R=0.4 - see the file-header comment for why this closure test is single-radius.
const int ir = 2;

// Same three stitched samples make_closure_trees.C's photon_scale map covers - see that
// file's comment (pythia-only; Photon5 has no herwig equivalent in this pipeline).
const vector<string> closureSamples = {"Photon5", "Photon10", "Photon20"};

// This closure's own pa scan window - deliberately NOT insitu_utility::scanLow/scanHigh
// ([0.90,1.00), tuned to Data's actual, known-to-be-small JES gap): make_closure_trees.C
// injects a scale anywhere across [0.9,1.1] by default, so the scan here needs its own,
// wider window to be able to find it, with some margin on both sides to still resolve
// the chi2 minimum's uncertainty band even for an injected scale near an edge.
const float paLow = 0.80, paHigh = 1.20;
const int na = 4000;

struct ClosureEvent { float pho_pt, jet_pt, weight; int ptbin; };

// Same event selection/re-indexing convention as insitu_utility::cacheDataEvents, plus
// the "weight" branch that function deliberately skips (see its own comment - real
// Data's weight is always 1, so it's a no-op there; not so for our stitched-MC halves).
vector<ClosureEvent> cacheClosureEvents(const char * filename, int abcdSelect, int ir) {
  vector<ClosureEvent> events;
  TFile * f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie()) {
    cout << "WARNING: could not open " << filename << endl;
    return events;
  }
  TTree * t = (TTree*)f->Get("insitutree");
  Float_t pho_pt, jet_pt, weight;
  Int_t abcd, evIr;
  t->SetBranchAddress("pho_pt", &pho_pt);
  t->SetBranchAddress("jet_pt", &jet_pt);
  t->SetBranchAddress("abcd", &abcd);
  t->SetBranchAddress("weight", &weight);
  t->SetBranchAddress("ir", &evIr);
  Long64_t nentries = t->GetEntries();
  for (Long64_t e = 0; e < nentries; e++) {
    t->GetEntry(e);
    if (abcd != abcdSelect) continue;
    if (evIr != ir) continue;
    int ipt = ana::findPtBin(pho_pt);
    if (ipt < ana::firstUsedPtBin || ipt >= ana::firstUsedPtBin + ana::nPtBinsUsed) continue;
    ipt -= ana::firstUsedPtBin;
    events.push_back({pho_pt, jet_pt, weight, ipt});
  }
  f->Close();
  return events;
}

// Weighted mean(x_J)/error per used photon-pT bin at a given trial pa - same
// weighted-mean/Kish-effective-sample-size error formula as
// insitu_utility::referenceMeans, just operating on an already-cached event list (so it
// can be called once per pa inside the grid scan below, like computeRegionAMeans, not
// just once as referenceMeans is for a fixed reference).
void weightedRegionAMeans(const vector<ClosureEvent> & events, float pa, float mean[], float err[], const float lowXj[]) {
  vector<double> sumw(nPtBinsUsed,0), sumw2(nPtBinsUsed,0), sumwx(nPtBinsUsed,0), sumwx2(nPtBinsUsed,0);
  for (auto & ev : events) {
    double x = (ev.jet_pt/pa)/ev.pho_pt;
    if (x < lowXj[ev.ptbin]) continue;
    double w = ev.weight;
    sumw[ev.ptbin]   += w;
    sumw2[ev.ptbin]  += w*w;
    sumwx[ev.ptbin]  += w*x;
    sumwx2[ev.ptbin] += w*x*x;
  }
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    if (sumw[ipt] <= 0) { mean[ipt] = 0; err[ipt] = 0; continue; }
    double m    = sumwx[ipt]/sumw[ipt];
    double var  = sumwx2[ipt]/sumw[ipt] - m*m;
    double neff = sumw[ipt]*sumw[ipt]/sumw2[ipt]; // Kish effective sample size
    mean[ipt] = m;
    err[ipt]  = sqrt(std::max(var,0.)/neff);
  }
}

// One comparison page: top panel is mean(x_J) vs photon pT for the "sim" reference and
// the "data" half at its raw (pa=1, uncorrected) scale; bottom panel is the raw ratio
// (uncorrected data/sim) and the corrected ratio (data/sim at the recovered best-fit
// pa) - same two-panel layout/convention as insitu/grid_insitu.C's drawJESPage, adapted
// to this closure test's single radius and known-injected-vs-recovered pa instead of
// Data's actual (unknown) JES gap.
void drawMeanPage(TCanvas * c, const char * pdfPath,
    TGraphErrors * gSim, TGraphErrors * gDataRaw, TGraphErrors * gRatioRaw, TGraphErrors * gRatioCorr,
    double injectedScale, float pa, float paErrLow, float paErrHigh) {
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
  gSim->SetLineColor(kMagenta+1);
  gSim->SetMarkerColor(kMagenta+1);
  gSim->SetMarkerStyle(21);
  gSim->SetLineWidth(2);
  gSim->Draw("p same");
  gDataRaw->SetLineColor(kBlue);
  gDataRaw->SetMarkerColor(kBlue);
  gDataRaw->SetMarkerStyle(20);
  gDataRaw->SetLineWidth(2);
  gDataRaw->Draw("p same");
  TLegend * l1 = new TLegend(.55,.1,.85,.3);
  l1->SetLineWidth(0);
  l1->AddEntry(gSim, "\"Sim\" (unscaled reference)");
  l1->AddEntry(gDataRaw, "\"Data\" (injected scale, uncorrected)");
  l1->Draw();
  insitu_utility::drawSPhenixLabel({"In-situ JES closure test"}, {
      Form("Jet R=%.1f", ana::JetRs[ir]),
      Form("Injected p_{a} = %.4f", injectedScale)
    }, .18, .85, 16, p1->GetWh()/1.5);

  p2->cd();
  p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.2);
  p2->SetLeftMargin(.15);
  gPad->SetTicks(1,1);
  TH1F * frame2 = p2->DrawFrame(ana::ptBinsUsed[0], 0.90, ana::ptBinsUsed[nPtBinsUsed], 1.10);
  frame2->GetYaxis()->SetTitle("Data/Sim");
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
  l2->AddEntry(gRatioRaw,  "Uncorrected ratio");
  l2->AddEntry(gRatioCorr, "Corrected ratio");
  l2->Draw();

  float paErr = (paErrLow+paErrHigh)/2.0;
  TLatex jestext;
  jestext.SetNDC();
  jestext.SetTextColor(kRed);
  jestext.DrawLatex(.18,.28, Form("Recovered p_{a} = %.4f #pm %.4f", pa, paErr));

  c->SaveAs(pdfPath);
}

void fit_closure() {
  TH1::AddDirectory(kFALSE);

  // One file pair per stitched sample - make_closure_trees.C is now run once per
  // sample (see its own header comment / run_closure.sh) so the three can be launched
  // as parallel ROOT processes instead of one process looping all three sequentially.
  vector<string> dataFiles, simFiles;
  for (const string & s : closureSamples) {
    dataFiles.push_back(Form("%s/ClosureData_%s_insitu.root", closure_input_dir, s.c_str()));
    simFiles.push_back(Form("%s/ClosureSim_%s_insitu.root", closure_input_dir, s.c_str()));
  }

  // injectedScale is written identically into every sample's data file (see
  // make_closure_trees.C) - read it back from whichever one opens first, and sanity-
  // check the rest agree (a mismatch would mean run_closure.sh's three parallel
  // invocations were NOT given the same injected value, which breaks the whole
  // premise of combining them).
  double injectedScale = -1;
  for (const string & df : dataFiles) {
    TFile * ftruth = TFile::Open(df.c_str(), "READ");
    if (!ftruth || ftruth->IsZombie()) {
      cout << "fit_closure: could not open " << df << " - run make_closure_trees.C first." << endl;
      return;
    }
    TParameter<double> * pInjected = (TParameter<double>*)ftruth->Get("injectedScale");
    double thisScale = pInjected ? pInjected->GetVal() : -1;
    ftruth->Close();
    if (thisScale < 0) {
      cout << "fit_closure: no injectedScale found in " << df << " - was it written by make_closure_trees.C?" << endl;
      return;
    }
    if (injectedScale < 0) injectedScale = thisScale;
    else if (fabs(thisScale - injectedScale) > 1e-6) {
      cout << "WARNING: " << df << " has injectedScale=" << thisScale
           << ", inconsistent with " << injectedScale << " from an earlier sample - "
           << "were all three make_closure_trees.C runs given the same injected scale?" << endl;
    }
  }
  cout << "Injected closure JES scale (truth): " << injectedScale << endl;

  // Concatenate all three samples' Region-A, R=0.4 events - the actual stitch. Each row
  // already carries its sample's cross-section weight (see make_closure_trees.C), so a
  // plain concatenation (not a weighted merge) is correct: weightedRegionAMeans below
  // sums per-row weights itself.
  vector<ClosureEvent> dataA, simA;
  for (size_t is = 0; is < dataFiles.size(); is++) {
    vector<ClosureEvent> d = cacheClosureEvents(dataFiles[is].c_str(), 0, ir);
    vector<ClosureEvent> s = cacheClosureEvents(simFiles[is].c_str(), 0, ir);
    dataA.insert(dataA.end(), d.begin(), d.end());
    simA.insert(simA.end(), s.begin(), s.end());
  }
  cout << "R=" << ana::JetRs[ir] << ": \"data\" A=" << dataA.size() << " \"sim\" A=" << simA.size() << endl;
  if (dataA.empty() || simA.empty()) {
    cout << "fit_closure: empty sample at R=" << ana::JetRs[ir] << " - aborting." << endl;
    return;
  }

  // Low-xJ floor per used pT bin - same cut make_closure_trees.C's closure_check_pair
  // (and unfolder::check_pair) applies at floorScale=1, so this scan excludes exactly
  // the same low-xJ events the main pipeline would at the same jet radius.
  float lowXj[nPtBinsUsed];
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) lowXj[ipt] = insitu_utility::lowXjFloor(ir, ana::ptBinsUsed[ipt]);

  // Reference mean(x_J): the closure "sim" half's own cross-section-weighted mean at
  // pa=1 - it is never rescaled (see make_closure_trees.C), so this plays the role
  // grid_insitu.C's referenceMeans() plays for real Data's fixed MC reference.
  float refMean[nPtBinsUsed], refMeanErr[nPtBinsUsed];
  weightedRegionAMeans(simA, 1.0, refMean, refMeanErr, lowXj);

  TGraph * gchisq = new TGraph(na);
  gchisq->SetName(Form("gchisq_%s", ana::rnames[ir]));
  gchisq->SetTitle(";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});#chi^{2}");

  float minchisq = FLT_MAX, minpa = 1;
  int ibest = 0;
  for (int ia = 0; ia < na; ia++) {
    float pa = paLow + ia*(paHigh-paLow)/na;
    float mean[nPtBinsUsed], err[nPtBinsUsed];
    weightedRegionAMeans(dataA, pa, mean, err, lowXj);
    float chisq = 0;
    for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
      if (refMean[ipt] <= 0 || mean[ipt] <= 0) continue;
      float diff = 1 - mean[ipt]/refMean[ipt];
      float errt = sqrt((err[ipt]*err[ipt])/(refMean[ipt]*refMean[ipt])
          + mean[ipt]*mean[ipt]*refMeanErr[ipt]*refMeanErr[ipt]/pow(refMean[ipt],4));
      if (errt > 0) chisq += diff*diff/(errt*errt);
    }
    gchisq->SetPoint(ia, pa, chisq);
    if (chisq < minchisq) { minchisq = chisq; minpa = pa; ibest = ia; }
  }
  float errLow, errHigh;
  insitu_utility::findError(gchisq, ibest, minchisq, errLow, errHigh);
  float errAvg = (errLow+errHigh)/2.0;
  float pull = (errAvg > 0 ? (minpa-injectedScale)/errAvg : 0);

  cout << "recovered p_a = " << minpa << " +" << errHigh << "/-" << errLow
       << "  (chi2=" << minchisq << ", injected=" << injectedScale << ", pull=" << pull << ")" << endl;

  // Raw (pa=1, i.e. exactly the injected scale, uncorrected) and best-fit-corrected
  // (pa=minpa) "data" means, recomputed here rather than reused from inside the scan
  // loop above (minpa generally doesn't land exactly on a loop iteration's own pa) -
  // same convention as grid_insitu.C's own post-scan mean recomputation.
  float rawMean[nPtBinsUsed], rawErr[nPtBinsUsed];
  float corrMean[nPtBinsUsed], corrErr[nPtBinsUsed];
  weightedRegionAMeans(dataA, 1.0,   rawMean,  rawErr,  lowXj);
  weightedRegionAMeans(dataA, minpa, corrMean, corrErr, lowXj);

  TGraphErrors * gSim       = insitu_utility::meanGraph(refMean, refMeanErr, "gMeanSim");
  TGraphErrors * gDataRaw   = insitu_utility::meanGraph(rawMean, rawErr, "gMeanData_raw");
  TGraphErrors * gRatioRaw  = insitu_utility::ratioGraph(rawMean,  rawErr,  refMean, refMeanErr, "gRatio_uncorrected");
  TGraphErrors * gRatioCorr = insitu_utility::ratioGraph(corrMean, corrErr, refMean, refMeanErr, "gRatio_corrected");

  string pdfPathStr = Form("%s/fit_closure.pdf", closure_pdf_dir);
  TCanvas * c = new TCanvas("c","",700,700);
  c->SaveAs(Form("%s[", pdfPathStr.c_str()));

  gPad->SetLeftMargin(.15);
  gPad->SetTicks(1,1);
  gchisq->Draw("AL");
  double ymax = gchisq->GetHistogram()->GetMaximum();
  TLine * vline = new TLine(injectedScale, 0, injectedScale, ymax);
  vline->SetLineColor(kRed);
  vline->SetLineStyle(9);
  vline->SetLineWidth(2);
  vline->Draw("same");
  insitu_utility::drawSPhenixLabel({"In-situ JES closure test", Form("Jet R=%.1f", ana::JetRs[ir])}, {
      Form("Injected p_{a} = %.4f (red)", injectedScale),
      Form("Recovered p_{a} = %.4f +%.4f/-%.4f", minpa, errHigh, errLow)
    }, .18, .85, 16, c->GetWh());
  c->SaveAs(pdfPathStr.c_str());

  drawMeanPage(c, pdfPathStr.c_str(), gSim, gDataRaw, gRatioRaw, gRatioCorr, injectedScale, minpa, errLow, errHigh);

  c->SaveAs(Form("%s]", pdfPathStr.c_str()));

  string outfilename = Form("%s/fit_closure.root", closure_output_dir);
  TFile * fout = TFile::Open(outfilename.c_str(), "RECREATE");
  gchisq->Write();
  gSim->Write();
  gDataRaw->Write();
  gRatioRaw->Write();
  gRatioCorr->Write();

  TTree * wt = new TTree("results", "closure test result (R=0.4)");
  double w_injected = injectedScale;
  int w_ir = ir; float w_R = ana::JetRs[ir], w_pa = minpa, w_errLow = errLow, w_errHigh = errHigh, w_pull = pull;
  wt->Branch("injectedScale", &w_injected);
  wt->Branch("ir", &w_ir);
  wt->Branch("JetR", &w_R);
  wt->Branch("pa", &w_pa);
  wt->Branch("errLow", &w_errLow);
  wt->Branch("errHigh", &w_errHigh);
  wt->Branch("pull", &w_pull);
  wt->Fill();
  wt->Write();
  fout->Close();

  cout << "\nWrote " << pdfPathStr << endl;
  cout << "Wrote " << outfilename << endl;
}
