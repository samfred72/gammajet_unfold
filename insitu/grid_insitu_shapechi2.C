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
// load it explicitly (see grid_insitu.C) so cling resolves the real compiled
// definitions instead of misbinding against the sibling gammajet project's own
// ana/drawer classes on the same library path.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Shape-chi2 variant of grid_insitu.C's in-situ JES study. grid_insitu.C picks the
// trial jet-energy-scale factor pa that makes Data's mean(x_{J#gamma}) match the fixed
// Pythia8 reference mean - a single number per pT bin, blind to any shape difference
// that happens to preserve the mean. This version instead picks pa by minimizing a
// bin-by-bin shape chi2 between Data's and the fixed MC's xJ distributions, each
// shape-normalized to unit area (bin fraction, not bin-width density - see
// referenceShape() below for why), summed over xJ bins (and pT bins) as a standard
// sum-of-pulls-squared, chi2 = sum_i [(f_i^data - f_i^MC)/sigma_i]^2. This is a
// genuinely different, shape-sensitive fit criterion, not just a display change -
// everything else (Region A alone vs. the purity-corrected A-C combination, the pa scan
// range in insitu_utility.h, the Delta-chi2=1 error convention, and all display pages)
// is unchanged from grid_insitu.C, so the two best-fit pa values can be compared
// directly to see how much the fit criterion itself matters.
//
// Region A and Region C jets are scaled by the *same* trial pa at every grid point,
// exactly as in grid_insitu.C. The purity P is computed once, outside the pa loop, and
// held fixed across the whole scan (same reasoning as grid_insitu.C's file header).

// insitu/ is split into inputs/ (the raw Data/Photon insitu ntuples, written by
// unfolder.h's production pipeline), output/ (this and the other grid_insitu*.C
// macros' own .root output), and pdfs/ (their .pdf output).
const char * insitu_input_dir  = "/home/samson72/sphnx/gammajet_unfold/insitu/inputs";
const char * insitu_output_dir = "/home/samson72/sphnx/gammajet_unfold/insitu/output";
const char * insitu_pdf_dir    = "/home/samson72/sphnx/gammajet_unfold/insitu/pdfs";

// Only ana::ptBinsUsed (15-20, 20-25, 25-35 GeV) is used for every calculation and
// plot below - same restriction as grid_insitu.C (see ana.h's
// ptBins/ptBinsUsed/firstUsedPtBin comment).
const int nPtBinsUsed = ana::nPtBinsUsed;

// The last 3 xJ bins in each pT bin have very low counts, so a per-bin shape chi2
// computed against them is dominated by their noise rather than genuine shape
// agreement - excluded here, same explicit-threshold precedent as
// drawing/draw_purity_corrected.C's nXjBinsForChi2 and
// drawing/draw_covariance_chi2.C:60 (see gammajet_unfold/CLAUDE.md's bin-selection
// ground rule: exclude low-count bins only with a stated, explicit threshold).
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;

// pT bin 2 (25-35 GeV) has only ~100 raw Region-A events (NA), so the native 0.1-wide
// ana::unfoldXjBins granularity gives too few events per bin within the fit range
// (low-xJ floor to nXjBinsForChi2) for a stable per-bin shape-chi2 pull - a single event
// migrating between adjacent bins as pa scans can swing that bin's pull by several
// sigma and dominate the total chi2 (see debug_shapechi2_spike.C /
// debug_shapechi2_sawtooth.C, which trace exactly this). For the shape-chi2 comparison
// only (not the mean-based fit in grid_insitu.C, not the stored/displayed histograms
// below), pT bin 2's xJ bins are merged pairwise, from its own low-xJ floor bin through
// the shared nXjBinsForChi2 tail cutoff - explicit, stated rebinning per this project's
// bin-selection ground rule (see draw_covariance_chi2.C:60). Every other pT bin keeps
// the native binning.
const int coarseRebinPtBin = 2;
const int coarseGroupSize = 2;

// coarsenSum/coarsenQuadrature, referenceMeans, referenceShape, computeRegionAMeans,
// computeCorrectedMeans, buildXjByPtBin, buildMCXjByPtBin, and purityCorrectByPtBin now
// live in src/insitu_utility.h/.cc (insitu_utility:: namespace) - moved there after
// being found copy-pasted byte-for-byte across all six grid_insitu*.C macros (see that
// header's comment).

// Cross-section weights for combining the Photon5/10/20 MC samples - same numbers as
// drawer.h's scalemap[isphoton=1][sample] for sim="pythia".
map<int,double> photon_scale = {{5,146359.3},{10,6944.675},{20,130.4461}};

// One comparison page: top panel is mean(x_J) vs pT for MC and raw Data; bottom panel
// is the raw ratio (raw Data/MC) and the corrected ratio (best-fit-scaled Data/MC) -
// unchanged from grid_insitu.C. The best-fit pa stamped here now comes from the shape
// chi2 scan below, not from matching this mean - so the corrected ratio is no longer
// guaranteed to sit exactly at 1 by construction (a useful cross-check in itself: how
// much does shape-matching still bring the mean into line?).
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
      Form("|#eta^{#gamma}|<%.1f, |#eta^{jet}|<%.1f", ana::photonEtaCut, ana::etacut-ana::JetRs[ir]),
      Form("#Delta#phi>%.0f#pi/%.0f", ana::oppnum, ana::oppden)
    }, .18, .85, 16, p1->GetWh()/1.5);

  p2->cd();
  p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.2);
  p2->SetLeftMargin(.15);
  gPad->SetTicks(1,1);
  TH1F * frame2 = p2->DrawFrame(ana::ptBinsUsed[0], 0.80, ana::ptBinsUsed[nPtBinsUsed], 1.10);
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
  jestext.DrawLatex(.18,.28, Form("Data to MC JES (shape #chi^{2}) = %.4f #pm %.4f", pa, paErr));

  c->SaveAs(pdfPath);
}

// One xJ-distribution comparison page, for a single photon-pT bin - unchanged from
// grid_insitu.C (still shape-normalized-plus-bin-width-density for display, per this
// project's plotting convention; the fit itself uses raw bin fractions, see
// referenceShape() above for why that's a separate choice from this display).
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
      Form("|#eta^{#gamma}|<%.1f, |#eta^{jet}|<%.1f", ana::photonEtaCut, ana::etacut-ana::JetRs[ir]),
      Form("#Delta#phi>%.0f#pi/%.0f", ana::oppnum, ana::oppden)
    }, .18, .85, 16, gPad->GetWh());

  c->SaveAs(pdfPath);
}

void grid_insitu_shapechi2(string systag = "nominal") {
  // See grid_insitu.C's identical comment: avoids per-pT-bin histograms auto-
  // registering into whichever radius subdirectory was left current by the previous
  // iteration's mkdir/cd.
  TH1::AddDirectory(kFALSE);

  // One file/one PDF for the whole systag, all seven jet radii inside - see
  // grid_insitu.C's identical comment for the mkdir/cd-before-TTree-construction detail.
  string pdfPathStr = Form("%s/grid_insitu_shapechi2_%s.pdf", insitu_pdf_dir, systag.c_str());
  TCanvas * c = new TCanvas("c","",700,700);
  c->SaveAs(Form("%s[", pdfPathStr.c_str()));

  string outfilename = Form("%s/grid_insitu_shapechi2_%s.root", insitu_output_dir, systag.c_str());
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

  vector<pair<string,double>> mcSamples = {
    {insitu_utility::insituFilename(insitu_input_dir, "Photon5",  "pythia", systag), photon_scale[5]},
    {insitu_utility::insituFilename(insitu_input_dir, "Photon10", "pythia", systag), photon_scale[10]},
    {insitu_utility::insituFilename(insitu_input_dir, "Photon20", "pythia", systag), photon_scale[20]},
  };

  float refMean[nPtBinsUsed], refMeanErr[nPtBinsUsed];
  insitu_utility::referenceMeans(mcSamples, 0, ir, refMean, refMeanErr, lowXj);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    cout << "MC reference <x_J> pt bin " << ipt << " [" << ana::ptBinsUsed[ipt] << "," << ana::ptBinsUsed[ipt+1]
         << "): " << refMean[ipt] << " +/- " << refMeanErr[ipt] << endl;
  }

  // Fixed MC reference xJ SHAPE (bin fractions) - this is what the pa scan below
  // actually fits to, replacing grid_insitu.C's refMean/refMeanErr as the fit target.
  vector<vector<double>> refFrac, refFracErr;
  insitu_utility::referenceShape(mcSamples, 0, ir, refFrac, refFracErr, lowXj);

  // Purity per photon-pT bin (region A and region C) - computed once, held fixed across
  // the whole pa scan, same as grid_insitu.C.
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
  // Grid scan: single overall jet-energy-scale factor pa, no pT-dependence - same scan
  // window/step (insitu_utility.h's scanLow/scanHigh/scanN) as grid_insitu.C, but the
  // chi2 minimized at each point is now a shape chi2 (sum of squared pulls between
  // Data's and MC's shape-normalized xJ bin fractions), not the old mean(x_J)-matching
  // chi2.
  // -----------------------------
  const int na = insitu_utility::scanN;
  const float lowa = insitu_utility::scanLow, higha = insitu_utility::scanHigh;

  TGraph * gchisqA    = new TGraph(na);
  TGraph * gchisqCorr = new TGraph(na);
  gchisqA->SetName("gchisq_regionA");
  gchisqA->SetTitle(";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});Shape #chi^{2}");
  gchisqCorr->SetName("gchisq_puritycorrected");
  gchisqCorr->SetTitle(";p_{a} (jet_{pt,corrected} = jet_{pt}/p_{a});Shape #chi^{2}");

  float minchisqA = FLT_MAX, minpaA = 1;
  float minchisqCorr = FLT_MAX, minpaCorr = 1;
  int ibestA = 0, ibestCorr = 0;

  for (int ia = 0; ia < na; ia++) {
    float pa = lowa + ia*(higha-lowa)/na;

    // Per-(pT bin, xJ bin) raw (unweighted) Data counts at this trial pa - the shape-chi2
    // analogue of the scalar sumA/sumA2/countA accumulators in grid_insitu.C's own loop.
    vector<vector<double>> countA(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
    vector<vector<double>> countC(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
    for (auto & ev : dataA) {
      float x = (ev.jet_pt/pa)/ev.pho_pt;
      if (x < lowXj[ev.ptbin]) continue;
      int ixj = ana::findUnfoldXjBin(x);
      if (ixj < 0 || ixj >= ana::nUnfoldXjBins) continue;
      countA[ev.ptbin][ixj] += 1;
    }
    for (auto & ev : dataC) {
      float x = (ev.jet_pt/pa)/ev.pho_pt;
      if (x < lowXj[ev.ptbin]) continue;
      int ixj = ana::findUnfoldXjBin(x);
      if (ixj < 0 || ixj >= ana::nUnfoldXjBins) continue;
      countC[ev.ptbin][ixj] += 1;
    }

    float chisqA = 0, chisqCorr = 0;
    for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
      double NA = 0;
      for (double n : countA[ipt]) NA += n;

      // pT bin 2 fits on coarsened bins (see coarseRebinPtBin above); every other pT
      // bin uses the native fine binning, i.e. useA/useRefFrac/useRefFracErr are just
      // the original per-bin arrays and nBins is nXjBinsForChi2 (unchanged behavior).
      vector<double> useA, useC, useRefFrac, useRefFracErr;
      int nBins;
      if (ipt == coarseRebinPtBin) {
        int startBin = ana::findUnfoldXjBin(lowXj[ipt]);
        useA         = insitu_utility::coarsenSum(countA[ipt], startBin, nXjBinsForChi2, coarseGroupSize);
        useC         = insitu_utility::coarsenSum(countC[ipt], startBin, nXjBinsForChi2, coarseGroupSize);
        useRefFrac   = insitu_utility::coarsenSum(refFrac[ipt], startBin, nXjBinsForChi2, coarseGroupSize);
        useRefFracErr = insitu_utility::coarsenQuadrature(refFracErr[ipt], startBin, nXjBinsForChi2, coarseGroupSize);
        nBins = (int)useA.size();
      } else {
        useA = countA[ipt]; useC = countC[ipt];
        useRefFrac = refFrac[ipt]; useRefFracErr = refFracErr[ipt];
        nBins = nXjBinsForChi2;
      }

      // Region A only (naive, background-contaminated fit): pull-squared per xJ bin
      // between Data's raw bin fraction and the fixed MC reference fraction, summed
      // over the first nBins bins (low-stat tail dropped, see nXjBinsForChi2 above).
      //
      // errt is floored at 1/NA: sqrt(countA)/NA is the usual Poisson error on a bin
      // fraction, but it goes to exactly 0 when countA=0, understating what a zero-count
      // observation actually leaves open - a Poisson process with a real, nonzero rate
      // routinely produces a handful of zero-count bins (e.g. rate~9 in the 25-35 GeV
      // pT slice's sparsest bins - Region A there has only NA~100 events total), so an
      // exactly-zero error is never justified. 1/NA is the coarsest resolvable step in
      // a bin fraction built from NA raw counts, so no bin can claim to be known finer
      // than that regardless of what the naive sqrt(count) formula says. Without this,
      // a single zero-count bin sitting next to a well-populated MC reference bin can
      // produce a pull of -70+ from one event's worth of statistical noise (see
      // debug_shapechi2_spike.C).
      if (NA > 0) {
        double errFloor = 1.0/NA;
        for (int ib = 0; ib < nBins; ib++) {
          double fData = useA[ib]/NA;
          double errData = sqrt(useA[ib])/NA;
          double errt = sqrt(errData*errData + useRefFracErr[ib]*useRefFracErr[ib]);
          errt = std::max(errt, errFloor);
          double diff = fData - useRefFrac[ib];
          chisqA += diff*diff/(errt*errt);
        }
      }

      // Purity-corrected (two-purity method, both regions scaled by the same pa) - same
      // coeffA/coeffC as grid_insitu.C's computeCorrectedMeans(), now applied to the
      // per-xJ-bin counts instead of Sum(x_J)/Sum(x_J^2) (still linear, see
      // unfold_utility::purityCorrectCoeffs). Stat error on the corrected fraction comes
      // from A/C counting stats alone (Poisson on the raw counts before the linear
      // combination) - purity's own uncertainty is held fixed across this scan, same
      // simplification as grid_insitu.C's mean-based chi2.
      double NC = 0;
      for (double n : countC[ipt]) NC += n;
      if (NA > 0 && NC > 0) {
        float coeffA, coeffC;
        unfold_utility::purityCorrectCoeffs(purity[ipt], purityC[ipt], NA, NC, coeffA, coeffC);
        double Ncorr = coeffA*NA - coeffC*NC;
        if (Ncorr > 0) {
          // Same 1/Ncorr error floor as the Region A block above, and for the same
          // reason - a zero-count bin in either region can drive errCorr to 0 even
          // though it says nothing more than "at most 1/Ncorr of the corrected sample".
          double errFloor = 1.0/Ncorr;
          for (int ib = 0; ib < nBins; ib++) {
            double corrCount = coeffA*useA[ib] - coeffC*useC[ib];
            double fCorr = corrCount/Ncorr;
            double errCorrCount = sqrt(coeffA*coeffA*useA[ib] + coeffC*coeffC*useC[ib]);
            double errCorr = errCorrCount/Ncorr;
            double errt = sqrt(errCorr*errCorr + useRefFracErr[ib]*useRefFracErr[ib]);
            errt = std::max(errt, errFloor);
            double diff = fCorr - useRefFrac[ib];
            chisqCorr += diff*diff/(errt*errt);
          }
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

  cout << "\nFINAL RESULT (jet R=" << ana::JetRs[ir] << ", systag=" << systag << ", shape chi2)\n";
  cout << "Region A only:        p_a = " << minpaA
       << " +" << errHighA << "/-" << errLowA << " (chi2=" << minchisqA << ")" << endl;
  cout << "Purity-corrected:     p_a = " << minpaCorr
       << " +" << errHighCorr << "/-" << errLowCorr << " (chi2=" << minchisqCorr << ")" << endl;

  // -----------------------------
  // Build x_J histograms per photon-pT bin: the fixed MC reference, Data at pa=1 (raw)
  // and at each study's own best-fit pa (corrected) - unchanged from grid_insitu.C.
  // -----------------------------
  vector<TH1D*> hxjMC_pt       = insitu_utility::buildMCXjByPtBin(mcSamples, 0, ir, "hxjA_pythia", lowXj);
  vector<TH1D*> hxjA_raw_pt    = insitu_utility::buildXjByPtBin(dataA, 1.0,       nPtBinsUsed, "hxjA_data_raw", lowXj);
  vector<TH1D*> hxjA_corr_pt   = insitu_utility::buildXjByPtBin(dataA, minpaA,    nPtBinsUsed, "hxjA_data_corr", lowXj);
  vector<TH1D*> hxjC_raw_pt    = insitu_utility::buildXjByPtBin(dataC, 1.0,       nPtBinsUsed, "hxjC_data_raw", lowXj);
  vector<TH1D*> hxjA_atCorr_pt = insitu_utility::buildXjByPtBin(dataA, minpaCorr, nPtBinsUsed, "hxjA_data_atCorrScale", lowXj);
  vector<TH1D*> hxjC_atCorr_pt = insitu_utility::buildXjByPtBin(dataC, minpaCorr, nPtBinsUsed, "hxjC_data_atCorrScale", lowXj);
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
  // Mean(x_J) vs pT comparison plots - display only now (the fit no longer targets
  // this mean) - top: MC vs raw Data; bottom: raw ratio vs corrected ratio (evaluated
  // at each study's own shape-chi2 best-fit pa).
  // -----------------------------
  gStyle->SetOptStat(0);

  float rawMeanA[nPtBinsUsed], rawErrA[nPtBinsUsed];
  float corrMeanA[nPtBinsUsed], corrErrA[nPtBinsUsed];
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
  // Save - see grid_insitu.C's identical comment.
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

  TTree * wt = new TTree("results", "best-fit jet energy scale results (shape chi2)");
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
