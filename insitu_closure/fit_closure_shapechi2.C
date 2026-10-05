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

// Shape-chi2 variant of fit_closure.C, same relationship as insitu/grid_insitu.C to
// insitu/grid_insitu_shapechi2.C: fit_closure.C picks the trial jet-energy-scale factor
// pa that makes the "data" half's mean(x_J) match the "sim" half's - a single number per
// pT bin, blind to any shape difference that happens to preserve the mean. This version
// instead picks pa by minimizing a bin-by-bin shape chi2 between the "data" half's and
// "sim" half's xJ distributions, each shape-normalized to unit area (bin fraction, not
// bin-width density), summed over xJ bins (and pT bins) as chi2 = sum_i [(f_i^data -
// f_i^sim)/sigma_i]^2 - a genuinely different, shape-sensitive fit criterion. Everything
// else (R=0.4 only, Region A only, the stitched Photon5+10+20 combination, the pa scan
// window/step) is unchanged from fit_closure.C, so the two best-fit pa values can be
// compared directly to see how much the fit criterion itself matters - e.g. whether a
// shape fit recovers the injected scale better or worse than mean-matching in the
// low-pT bins, where the low-xJ floor truncates a large fraction of the distribution
// (see the mean version's own low-pT-bin discussion).
//
// Deliberately does NOT reuse insitu_utility::referenceShape: like fit_closure.C's own
// weightedRegionAMeans (vs. cacheDataEvents/computeRegionAMeans), referenceShape reads
// straight from a file once for a fixed reference - here BOTH halves are stitched,
// cross-section-weighted MC, and the "data" half's shape has to be recomputed at every
// trial pa from an already-cached in-memory event list, not re-read from disk per grid
// point. weightedRegionAShape below is that one weighted, re-callable implementation,
// used for both the fixed "sim" reference (pa=1, called once) and the "data" scan
// (called once per trial pa).

const char * closure_input_dir  = ana::path("insitu_closure/inputs");
const char * closure_output_dir = ana::path("insitu_closure/output");
const char * closure_pdf_dir    = ana::path("insitu_closure/pdfs");

const int nPtBinsUsed = ana::nPtBinsUsed;

// R=0.4 - see fit_closure.C's header comment for why this closure test is single-radius.
const int ir = 2;

// Same three stitched samples make_closure_trees.C's photon_scale map covers.
const vector<string> closureSamples = {"Photon5", "Photon10", "Photon20"};

// Same pa scan window as fit_closure.C - see that file's comment.
const float paLow = 0.80, paHigh = 1.20;
const int na = 4000;

// The last 3 xJ bins in each pT bin have very low counts, so a per-bin shape chi2
// computed against them is dominated by their noise rather than genuine shape
// agreement - same explicit-threshold precedent as insitu/grid_insitu_shapechi2.C and
// drawing/draw_covariance_chi2.C:60 (gammajet_unfold/CLAUDE.md's bin-selection ground
// rule: exclude low-count bins only with a stated, explicit threshold).
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;

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

// Weighted mean(x_J)/error per used photon-pT bin at a given trial pa - unchanged from
// fit_closure.C, kept here for the display-only mean/ratio comparison page (see that
// file's comment: the fit criterion below is shape, not mean, but plotting the mean
// anyway shows how much shape-matching still brings it into line - same cross-check
// insitu/grid_insitu_shapechi2.C keeps its mean panel for).
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

// Weighted xJ SHAPE (bin fraction, not density) and its per-bin error, per used
// photon-pT bin, at a given trial pa - the shape-chi2 analogue of weightedRegionAMeans
// above / the weighted counterpart of insitu_utility::referenceShape (same
// sqrt(sum w^2)/N bin-fraction-error convention). `neff[ipt]` is that pT bin's Kish
// effective sample size (sum(w)^2/sum(w^2), same formula weightedRegionAMeans uses) -
// the weighted generalization of a raw event count, needed below for a per-bin error
// floor: a bin fraction built from Neff effectively-independent events can't be known
// finer than 1/Neff, same reasoning as (and same failure mode without it as)
// grid_insitu_shapechi2.C's 1/NA floor on its raw, unweighted Data counts.
void weightedRegionAShape(const vector<ClosureEvent> & events, float pa,
    vector<vector<double>> & frac, vector<vector<double>> & fracErr, vector<double> & neff, const float lowXj[]) {
  vector<vector<double>> sumw(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
  vector<vector<double>> sumw2(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
  for (auto & ev : events) {
    double x = (ev.jet_pt/pa)/ev.pho_pt;
    if (x < lowXj[ev.ptbin]) continue;
    int ixj = ana::findUnfoldXjBin(x);
    if (ixj < 0 || ixj >= ana::nUnfoldXjBins) continue;
    double w = ev.weight;
    sumw[ev.ptbin][ixj]  += w;
    sumw2[ev.ptbin][ixj] += w*w;
  }
  frac.assign(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
  fracErr.assign(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
  neff.assign(nPtBinsUsed, 0.);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    double N = 0, N2 = 0;
    for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) { N += sumw[ipt][ixj]; N2 += sumw2[ipt][ixj]; }
    if (N <= 0) continue;
    neff[ipt] = N*N/N2;
    for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
      frac[ipt][ixj]    = sumw[ipt][ixj]/N;
      fracErr[ipt][ixj] = sqrt(sumw2[ipt][ixj])/N;
    }
  }
}

// Same two-panel mean/ratio comparison page as fit_closure.C's drawMeanPage, duplicated
// here (rather than shared) per this project's precedent of self-contained ROOT macros
// (e.g. insitu/grid_insitu.C and grid_insitu_shapechi2.C each carry their own copy of
// drawJESPage/drawXjPage) - only the axis label changes (this study's pa comes from a
// shape fit, not a mean fit).
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
  insitu_utility::drawSPhenixLabel({"In-situ JES closure test (shape #chi^{2})"}, {
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
  l2->AddEntry(gRatioCorr, "Corrected ratio (shape #chi^{2} best fit)");
  l2->Draw();

  float paErr = (paErrLow+paErrHigh)/2.0;
  TLatex jestext;
  jestext.SetNDC();
  jestext.SetTextColor(kRed);
  jestext.DrawLatex(.18,.28, Form("Recovered p_{a} (shape #chi^{2}) = %.4f #pm %.4f", pa, paErr));

  c->SaveAs(pdfPath);
}

void fit_closure_shapechi2() {
  TH1::AddDirectory(kFALSE);

  vector<string> dataFiles, simFiles;
  for (const string & s : closureSamples) {
    dataFiles.push_back(Form("%s/ClosureData_%s_insitu.root", closure_input_dir, s.c_str()));
    simFiles.push_back(Form("%s/ClosureSim_%s_insitu.root", closure_input_dir, s.c_str()));
  }

  // injectedScale is written identically into every sample's data file (see
  // make_closure_trees.C) - read it back from whichever one opens first, and sanity-
  // check the rest agree, same as fit_closure.C.
  double injectedScale = -1;
  for (const string & df : dataFiles) {
    TFile * ftruth = TFile::Open(df.c_str(), "READ");
    if (!ftruth || ftruth->IsZombie()) {
      cout << "fit_closure_shapechi2: could not open " << df << " - run make_closure_trees.C first." << endl;
      return;
    }
    TParameter<double> * pInjected = (TParameter<double>*)ftruth->Get("injectedScale");
    double thisScale = pInjected ? pInjected->GetVal() : -1;
    ftruth->Close();
    if (thisScale < 0) {
      cout << "fit_closure_shapechi2: no injectedScale found in " << df << " - was it written by make_closure_trees.C?" << endl;
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

  // Concatenate all three samples' Region-A, R=0.4 events - the actual stitch.
  vector<ClosureEvent> dataA, simA;
  for (size_t is = 0; is < dataFiles.size(); is++) {
    vector<ClosureEvent> d = cacheClosureEvents(dataFiles[is].c_str(), 0, ir);
    vector<ClosureEvent> s = cacheClosureEvents(simFiles[is].c_str(), 0, ir);
    dataA.insert(dataA.end(), d.begin(), d.end());
    simA.insert(simA.end(), s.begin(), s.end());
  }
  cout << "R=" << ana::JetRs[ir] << ": \"data\" A=" << dataA.size() << " \"sim\" A=" << simA.size() << endl;
  if (dataA.empty() || simA.empty()) {
    cout << "fit_closure_shapechi2: empty sample at R=" << ana::JetRs[ir] << " - aborting." << endl;
    return;
  }

  float lowXj[nPtBinsUsed];
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) lowXj[ipt] = insitu_utility::lowXjFloor(ir, ana::ptBinsUsed[ipt]);

  // Fixed "sim" reference shape (bin fractions) - never rescaled, computed once. This
  // is what the pa scan below actually fits to, replacing fit_closure.C's refMean/
  // refMeanErr as the fit target.
  vector<vector<double>> refFrac, refFracErr; vector<double> refNeff;
  weightedRegionAShape(simA, 1.0, refFrac, refFracErr, refNeff, lowXj);

  // -----------------------------
  // Grid scan: same window/step as fit_closure.C, but the chi2 minimized at each point
  // is now a shape chi2 (sum of squared pulls between "data"'s and "sim"'s
  // shape-normalized xJ bin fractions), not the mean(x_J)-matching chi2.
  // -----------------------------
  TGraph * gchisq = new TGraph(na);
  gchisq->SetName(Form("gchisq_shape_%s", ana::rnames[ir]));
  gchisq->SetTitle(";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});Shape #chi^{2}");

  float minchisq = FLT_MAX, minpa = 1;
  int ibest = 0;
  for (int ia = 0; ia < na; ia++) {
    float pa = paLow + ia*(paHigh-paLow)/na;
    vector<vector<double>> dataFrac, dataFracErr; vector<double> dataNeff;
    weightedRegionAShape(dataA, pa, dataFrac, dataFracErr, dataNeff, lowXj);

    float chisq = 0;
    for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
      if (dataNeff[ipt] <= 0) continue;
      // Error floor - see weightedRegionAShape's comment: a bin fraction built from
      // this pT bin's Neff effectively-independent "data" events can't be known finer
      // than 1/Neff, regardless of what a single (possibly zero-count) bin's own
      // sqrt(sum w^2)/N formula says.
      double errFloor = 1.0/dataNeff[ipt];
      for (int ixj = 0; ixj < nXjBinsForChi2; ixj++) {
        double diff = dataFrac[ipt][ixj] - refFrac[ipt][ixj];
        double errt = sqrt(dataFracErr[ipt][ixj]*dataFracErr[ipt][ixj] + refFracErr[ipt][ixj]*refFracErr[ipt][ixj]);
        errt = std::max(errt, errFloor);
        chisq += diff*diff/(errt*errt);
      }
    }
    gchisq->SetPoint(ia, pa, chisq);
    if (chisq < minchisq) { minchisq = chisq; minpa = pa; ibest = ia; }
  }
  float errLow, errHigh;
  insitu_utility::findError(gchisq, ibest, minchisq, errLow, errHigh);
  float errAvg = (errLow+errHigh)/2.0;
  float pull = (errAvg > 0 ? (minpa-injectedScale)/errAvg : 0);

  cout << "recovered p_a (shape chi2) = " << minpa << " +" << errHigh << "/-" << errLow
       << "  (chi2=" << minchisq << ", injected=" << injectedScale << ", pull=" << pull << ")" << endl;

  // Raw (pa=1) and best-fit-corrected (pa=minpa) "data" means - display only (the fit
  // criterion above is shape, not mean), same convention as
  // insitu/grid_insitu_shapechi2.C's own mean panel.
  float refMean[nPtBinsUsed], refMeanErr[nPtBinsUsed];
  weightedRegionAMeans(simA, 1.0, refMean, refMeanErr, lowXj);
  float rawMean[nPtBinsUsed], rawErr[nPtBinsUsed];
  float corrMean[nPtBinsUsed], corrErr[nPtBinsUsed];
  weightedRegionAMeans(dataA, 1.0,   rawMean,  rawErr,  lowXj);
  weightedRegionAMeans(dataA, minpa, corrMean, corrErr, lowXj);

  TGraphErrors * gSim       = insitu_utility::meanGraph(refMean, refMeanErr, "gMeanSim");
  TGraphErrors * gDataRaw   = insitu_utility::meanGraph(rawMean, rawErr, "gMeanData_raw");
  TGraphErrors * gRatioRaw  = insitu_utility::ratioGraph(rawMean,  rawErr,  refMean, refMeanErr, "gRatio_uncorrected");
  TGraphErrors * gRatioCorr = insitu_utility::ratioGraph(corrMean, corrErr, refMean, refMeanErr, "gRatio_corrected");

  string pdfPathStr = Form("%s/fit_closure_shapechi2.pdf", closure_pdf_dir);
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
  insitu_utility::drawSPhenixLabel({"In-situ JES closure test (shape #chi^{2})", Form("Jet R=%.1f", ana::JetRs[ir])}, {
      Form("Injected p_{a} = %.4f (red)", injectedScale),
      Form("Recovered p_{a} = %.4f +%.4f/-%.4f", minpa, errHigh, errLow)
    }, .18, .85, 16, c->GetWh());
  c->SaveAs(pdfPathStr.c_str());

  drawMeanPage(c, pdfPathStr.c_str(), gSim, gDataRaw, gRatioRaw, gRatioCorr, injectedScale, minpa, errLow, errHigh);

  c->SaveAs(Form("%s]", pdfPathStr.c_str()));

  string outfilename = Form("%s/fit_closure_shapechi2.root", closure_output_dir);
  TFile * fout = TFile::Open(outfilename.c_str(), "RECREATE");
  gchisq->Write();
  gSim->Write();
  gDataRaw->Write();
  gRatioRaw->Write();
  gRatioCorr->Write();

  TTree * wt = new TTree("results", "closure test result, shape chi2 (R=0.4)");
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
