#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/unfolder.h"
#include "../src/unfold_utility.h"
#include "../src/reweight_utility.h"
// The original gammajet project builds its OWN, differently-laid-out drawer/ana classes
// into /home/samson72/root/lib/libgammajet.so, sitting on the same library search path
// as this project's libgammajet_unfold.so. Without forcing which one loads first, ROOT's
// implicit symbol autoload can bind drawer/ana calls to the wrong (mismatched-layout)
// library and segfault - explicit load removes the ambiguity.
R__LOAD_LIBRARY(libgammajet_unfold.so);
R__LOAD_LIBRARY(libRooUnfold.so); // needed because we instantiate unfolder directly below

// Diagnostic: how much does the reconstructed jet's EM fraction (jet_emfrac[ir]) matter
// for the measured x_{J#gamma} shape, and how much of that is just Data/MC not agreeing
// on the emfrac distribution itself?
//
// Part 1 (pages 1-2): splits paired events into ana::emfracBins ({0, 0.5, 0.8, 1} -
// ana.h:175, already defined but never actually used anywhere in the codebase until this
// macro) and overlays the resulting x_J shapes, once for Data and once for the combined
// Photon5+10+20 pythia MC.
//
// Part 2 (pages 3-4): derives a per-emfrac-bin reweighting factor
// w(emfrac) = normalized_Data(emfrac) / normalized_MC(emfrac) (a finer, ana::emfracBins-
// independent binning - see fineEmfracBins below) from the same event samples, applies it
// as an extra per-event weight on top of MC's existing weight when refilling the x_J
// histogram, and compares MC's x_J shape before vs. after this reweighting (with Data
// shown alongside as the target the reweighting is meant to approach).
//
// Part 3 (page 5): <x_J> Data/MC ratio vs. photon p_T, one curve per ana::emfracBins slice
// plus one for the inclusive/"Total" sample, across every ana::ptBinsUsed bin - unlike
// Parts 1-2's single fixed bin, this scans the full reported pT range to see whether the
// emfrac dependence seen there is pT-dependent.
//
// Reuses unfolder::check_pair/check_keep_MC directly (by constructing a throwaway
// unfolder instance per sample and never calling fill_matrix/unfold/end on it) rather
// than reimplementing the pairing/eta/dphi/xJ-floor cuts by hand, so this can't silently
// drift from the production selection. Reconstructs the same nominal-systag reco-level
// quantities fill_matrix() would (unfolder.cc:149-232) for MC: the always-on nominal EMR smear
// on cluster pT and the truth-derived jet_pt_smear_truth[ir] jet pT. Data deliberately
// does NOT get unfolder.cc's per-radius JES correction (ana::jesNominal[ir]) here - uses
// raw jet_pt_calib[ir] directly, on explicit request, matching
// temporary_study/draw_xj_data_purity.C's convention rather than the production pipeline's.
//
// Purity correction: Data (never MC - it's signal-only pythia with no ABCD background to
// subtract) is purity-corrected via the same two-purity background subtraction as
// temporary_study/draw_xj_data_purity.C and insitu/grid_insitu.C - Region A
// (findabcdBin==0) events combined with their Region C (findabcdBin==2: good iso, bad
// bdt) siblings, via unfold_utility::purityCorrectCoeffs/purityCorrect, with purity
// P_A/P_C (+ErrorLow/ErrorHigh) from ana::getPurity/getPurityC at systag="nominal", ir=2 -
// i.e. the purity_nominal.root inputs, not re-derived here. Parts 1-2 purity-correct per
// emfrac slice (assuming P_A/P_C - a function of photon p_T only - doesn't itself depend
// on jet emfrac); Part 3 purity-corrects per (p_T bin, emfrac category) directly on the
// exact unbinned mean, not a binned histogram mean.
//
// Deliberately simplified relative to the real analysis in one remaining way, noted on
// every plot: R=0.4 (ir=2) only, matching this project's standalone-diagnostic convention
// (draw_response_matrix.C, draw_refolding.C, etc. all hardcode the same radius).
//
// MC Photon5/10/20 combination uses the exact same per-sample cross-section scale
// factors as the real pipeline (drawer::getScale(true, sample) - reads the same
// `scalemap` drawer::combineMC uses), applied on top of the same per-event vz/cluster-pT
// reweighting (src/reweight_utility.h) unfolder.cc itself applies - not re-derived or
// approximated. The Part 2 emfrac reweighting is a further, additional per-event weight
// on top of both of those, not a replacement for either.

const int ir = 2; // R=0.4, nominal - see header comment
const set<int> photonSampleCodes = {5, 10, 20};

// The single pT bin Parts 1-2 focus on - one of ana::ptBinsUsed's 7 reported bins
// (15-20 GeV), used on explicit request as a representative mid-range bin. Part 3 instead
// scans every ana::ptBinsUsed bin.
const double singleBinPtLo = ana::ptBins[5];
const double singleBinPtHi = ana::ptBins[6];

// Draws one color-matched "<varLabel> = mean" line per histogram, NDC coordinates,
// starting at (x,yTop) and stepping down by dy per line - meant to sit directly below an
// already-drawn TLegend in the same pad, one line per legend entry in the same order/
// color, matching the mean-labeling convention already used elsewhere in this project
// (fig:insituxjcomp / insitu/draw_insitu_xj.C). Means are read directly off each raw
// (weighted, not yet shape-normalized or bin-width-divided) histogram passed in -
// TH1::GetMean() is unaffected by an overall Scale(), but WOULD be distorted by
// densityForDisplay's per-bin division by (non-uniform) bin width, so this must be
// called with the pre-display histograms, not their _disp clones.
// width should match the legend's own width (x2-x1) these lines sit below, so the
// opaque backing box lines up with it rather than looking like a mismatched patch.
void drawMeans(const vector<TH1D*> & hMeanSrc, const vector<int> & colors,
    const char * varLabel, double x, double yTop, double dy = 0.045, double width = 0.35) {
  int n = (int)hMeanSrc.size();
  // Opaque backing box: the curve/marker content these labels sit "below the legend"
  // over varies bin-to-bin and page-to-page (peak position shifts with the plotted
  // quantity), so no fixed NDC position can be guaranteed clear of it - a solid box
  // drawn first, then text on top, keeps the labels legible regardless.
  TPave * bg = new TPave(x - 0.02, yTop - (n - 1) * dy - 0.04, x - 0.02 + width, yTop + 0.025, 0, "NDC");
  bg->SetFillColor(kWhite);
  bg->SetFillStyle(1001);
  bg->SetBorderSize(0);
  bg->Draw();
  for (int i = 0; i < n; i++) {
    TLatex * t = new TLatex(x, yTop - i * dy, Form("<%s> = %.3f", varLabel, hMeanSrc[i]->GetMean()));
    t->SetNDC();
    t->SetTextColor(colors[i]);
    t->SetTextSize(0.035);
    t->Draw();
  }
}

// One passing (paired, Region A or C, per check_pair) event's emfrac, x_J, photon pT, the
// mcWeight*crossSectionScale weight it already carries (1 for Data), and which ABCD region
// (0=A, 2=C) it landed in - collected once per sample across the full ana::ptBinsUsed
// range and reused for Parts 1-3 (filtered/binned differently by each), instead of looping
// the tree more than once per sample.
struct EmfracEvent { float emfrac; float xj; float weight; float phopt; int abcdRegion; };

// Loops one (trigger,sim) sample once and returns one EmfracEvent per event passing the
// same cuts drawEmfracPage's caller used to rely on fillEmfracXj for (every
// ana::ptBinsUsed-covered photon pT, R=0.4, check_pair). extraScale is the cross-section
// combination factor for MC (drawer::getScale(true, sample)), 1 for Data. includeRegionC
// additionally keeps Region C (findabcdBin==2) events alongside Region A - only ever
// passed true for Data, which needs its Region C sibling to purity-correct; MC is
// signal-only pythia with no ABCD background to subtract, so it's Region A only.
vector<EmfracEvent> collectEvents(string trigger, string sim, bool isMCsample, double extraScale,
    bool includeRegionC = false) {
  // A brace-init-list literal passed directly as this constructor's vector<string>
  // argument (e.g. unfolder uf(trigger, sim, {"nominal"});) crashes Cling at parse time
  // (reproduced and isolated in-session) - pass a pre-built vector instead.
  vector<string> nominalSystag = {"nominal"};
  unfolder uf(trigger, sim, nominalSystag);
  Reweighter rw;
  TRandom3 rnd(12345);
  vector<EmfracEvent> events;

  Long64_t nentries = uf.t->GetEntriesFast();
  cout << "  " << trigger << " (" << sim << "): " << nentries << " entries" << endl;
  for (Long64_t e = 0; e < nentries; e++) {
    uf.t->GetEntry(e);
    if (fabs(uf.vz) > ana::vzcut) continue;

    if (isMCsample) {
      vector<bool> keepMC = uf.check_keep_MC(uf.truth_cluster_pt, uf.cluster_pt,
          uf.truth_jet_pt, uf.jet_pt_smear_truth, trigger);
      if (!keepMC[ir]) continue;
    }

    float mcWeight = isMCsample ? rw.GetWeight(uf.vz, uf.cluster_pt) : 1.0f;

    // Nominal-systag reco-level photon: MC carries the always-on nominal EMR smear
    // (ana::emResolutionSigma nominal, PPG12 prescription - see ana.h);
    // Data is used as-is. Not the same random sequence unfolder.cc itself draws (its own
    // TRandom member is private, and bit-for-bit reproducibility isn't needed for a
    // standalone diagnostic), just the same smearing recipe.
    float recoClusterPt = uf.cluster_pt;
    if (isMCsample) recoClusterPt += rnd.Gaus(0, ana::emResolutionSigma(uf.truth_cluster_pt) * uf.truth_cluster_pt);

    pho_object maxpho(recoClusterPt, uf.cluster_e, uf.cluster_eta, uf.cluster_phi,
        uf.cluster_showershape[10], uf.cluster_showershape[11], uf.cluster_time,
        uf.cluster_bdt_scores[9],
        pho_object::get_showershape(uf.cluster_showershape, recoClusterPt));

    // No JES correction applied to Data here (on explicit request) - raw jet_pt_calib[ir],
    // not divided by ana::jesNominal[ir] as unfolder.cc's fill_matrix() does.
    float recoJetPt = isMCsample ? uf.jet_pt_smear_truth[ir] : uf.jet_pt_calib[ir];
    jet_object maxjet(recoJetPt, uf.jet_e[ir], uf.jet_eta[ir], uf.jet_phi[ir],
        uf.jet_emfrac[ir], 0, 0, uf.jet_time[ir]);

    if (!uf.check_pair(maxjet, ir, maxpho, true)) continue;
    // Every ana::ptBinsUsed-covered photon pT, not just Parts 1-2's single bin -
    // filterPtRange narrows this down for them; Part 3 uses the full range directly.
    if (maxpho.pt < ana::ptBinsUsed[0] || maxpho.pt >= ana::ptBinsUsed[ana::nPtBinsUsed]) continue;

    int iabcd = ana::findabcdBin(maxpho.iso4, maxpho.bdt, 0);
    if (iabcd != 0 && !(includeRegionC && iabcd == 2)) continue;

    events.push_back({maxjet.emfrac, (float)(maxjet.pt / maxpho.pt), mcWeight * (float)extraScale,
        maxpho.pt, iabcd});
  }
  return events;
}

// Combined Photon5+10+20 pythia sample, each cross-section-scaled via
// drawer::getScale(true, sample) - the same combination drawer::combineMC uses. Region A
// only (see collectEvents' includeRegionC comment).
vector<EmfracEvent> collectMCPhotonEvents() {
  drawer dScale; // only used for its scalemap accessor (getScale) - opens hists/*.root itself
  vector<EmfracEvent> all;
  for (int code : photonSampleCodes) {
    vector<EmfracEvent> ev = collectEvents(Form("Photon%d", code), "pythia", true, dScale.getScale(true, code));
    all.insert(all.end(), ev.begin(), ev.end());
  }
  return all;
}

vector<EmfracEvent> filterPtRange(const vector<EmfracEvent> & events, double ptlo, double pthi) {
  vector<EmfracEvent> out;
  for (const EmfracEvent & ev : events)
    if (ev.phopt >= ptlo && ev.phopt < pthi) out.push_back(ev);
  return out;
}

void splitByRegion(const vector<EmfracEvent> & events, vector<EmfracEvent> & regionA, vector<EmfracEvent> & regionC) {
  for (const EmfracEvent & ev : events) (ev.abcdRegion == 0 ? regionA : regionC).push_back(ev);
}

// Purity-corrects hA using its Region C sibling hC via the same two-purity method as
// temporary_study/draw_xj_data_purity.C (unfold_utility::purityCorrect, ana::getPurity/
// getPurityC at systag="nominal", this file's own ir). Falls back to raw Region A (hA
// itself has no statistics-cross-normalization problem, so the fallback just skips
// subtracting anything) if Region C has no statistics in [ptlo,pthi) -
// unfold_utility::purityCorrect already prints a WARNING in that case.
TH1D * purityCorrectHist(TH1D * hA, TH1D * hC, double ptlo, double pthi, const char * name) {
  float pA        = ana::getPurity(ptlo, pthi, "nominal", ir);
  float pAErrLow  = ana::getPurityErrorLow(ptlo, pthi, "nominal", ir);
  float pAErrHigh = ana::getPurityErrorHigh(ptlo, pthi, "nominal", ir);
  float pC        = ana::getPurityC(ptlo, pthi, "nominal", ir);
  float pCErrLow  = ana::getPurityCErrorLow(ptlo, pthi, "nominal", ir);
  float pCErrHigh = ana::getPurityCErrorHigh(ptlo, pthi, "nominal", ir);
  TH1D * hcorr = unfold_utility::purityCorrect(hA, hC, pA, pAErrLow, pAErrHigh, pC, pCErrLow, pCErrHigh, name);
  if (!hcorr) return (TH1D*)hA->Clone(name);
  return hcorr;
}

// --- Part 1: x_J split by the coarse ana::emfracBins, one page per sample ---

// MC (Region A only, raw - no background to subtract).
vector<TH1D*> binByEmfracCoarse(const vector<EmfracEvent> & events, const char * tag) {
  vector<TH1D*> h(3);
  for (int i = 0; i < 3; i++)
    h[i] = new TH1D(Form("hEmfracXj_%s_%d", tag, i), "", ana::nUnfoldXjBins, ana::unfoldXjBins);
  for (const EmfracEvent & ev : events) {
    int iemfrac = ana::findEmfracBin(ev.emfrac);
    if (iemfrac < 0) continue; // shouldn't happen (emfracBins spans [0,1] with no gap), but guard anyway
    h[iemfrac]->Fill(ev.xj, ev.weight);
  }
  return h;
}

// Data-only purity-corrected counterpart: bins Region A and Region C events separately by
// the same coarse emfrac slice, then purity-corrects each slice via purityCorrectHist -
// i.e. this assumes P_A/P_C (a function of photon p_T only, ana::getPurity/getPurityC)
// doesn't itself depend on jet emfrac, only that the emfrac split further partitions each
// region's own events.
vector<TH1D*> binByEmfracCoarsePurityCorrected(const vector<EmfracEvent> & dataEventsInBin, const char * tag) {
  vector<TH1D*> hA(3), hC(3), hCorr(3);
  for (int i = 0; i < 3; i++) {
    hA[i] = new TH1D(Form("hEmfracXjA_%s_%d", tag, i), "", ana::nUnfoldXjBins, ana::unfoldXjBins);
    hC[i] = new TH1D(Form("hEmfracXjC_%s_%d", tag, i), "", ana::nUnfoldXjBins, ana::unfoldXjBins);
  }
  for (const EmfracEvent & ev : dataEventsInBin) {
    int iemfrac = ana::findEmfracBin(ev.emfrac);
    if (iemfrac < 0) continue;
    (ev.abcdRegion == 0 ? hA : hC)[iemfrac]->Fill(ev.xj, ev.weight);
  }
  for (int i = 0; i < 3; i++) {
    hCorr[i] = purityCorrectHist(hA[i], hC[i], singleBinPtLo, singleBinPtHi, Form("hEmfracXjCorr_%s_%d", tag, i));
    delete hA[i];
    delete hC[i];
  }
  return hCorr;
}

// One shape-normalized, bin-width-divided x_J overlay page, one curve per coarse emfrac
// slice, for a given sample's three already-filled raw hIn[0..2] histograms.
void drawEmfracPage(TCanvas * c, const char * pdfPath, vector<TH1D*> & hIn,
    const vector<string> & samples, const vector<string> & extraFeatures) {
  const int colors[3] = {kBlue + 1, kGreen + 2, kRed + 1};
  const vector<string> labels = {
      Form("%.1f < emfrac < %.1f", ana::emfracBins[0], ana::emfracBins[1]),
      Form("%.1f < emfrac < %.1f", ana::emfracBins[1], ana::emfracBins[2]),
      Form("%.1f < emfrac < %.1f", ana::emfracBins[2], ana::emfracBins[3]),
  };

  vector<TH1D*> hDisp(3);
  double ymax = 0;
  for (int i = 0; i < 3; i++) {
    cout << "    " << labels[i] << ": " << hIn[i]->GetEntries() << " entries, "
         << hIn[i]->Integral() << " weighted" << endl;
    hDisp[i] = unfold_utility::densityForDisplay(hIn[i], Form("%s_disp", hIn[i]->GetName()));
    if (hDisp[i]->Integral() > 0) hDisp[i]->Scale(1. / hDisp[i]->Integral());
    // Sanitize NaN/Inf before drawing - poisons ROOT's axis auto-ranging otherwise.
    for (int b = 0; b <= hDisp[i]->GetNbinsX() + 1; b++) {
      double v = hDisp[i]->GetBinContent(b);
      if (std::isnan(v) || std::isinf(v)) hDisp[i]->SetBinContent(b, 0);
    }
    ymax = std::max(ymax, hDisp[i]->GetMaximum());
  }

  c->Clear();
  c->cd();
  gPad->SetTicks(1, 1);
  gPad->SetLeftMargin(.15);
  gPad->SetBottomMargin(.13);

  TLegend * leg = new TLegend(.5, .6, .85, .85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  for (int i = 0; i < 3; i++) {
    hDisp[i]->SetLineColor(colors[i]);
    hDisp[i]->SetLineWidth(2);
    hDisp[i]->GetXaxis()->SetTitle("x_{J#gamma}");
    hDisp[i]->GetYaxis()->SetTitle("(1/N) dN/dx_{J#gamma}");
    hDisp[i]->GetYaxis()->SetRangeUser(0, ymax * 1.4);
    hDisp[i]->Draw(i == 0 ? "hist" : "hist same");
    leg->AddEntry(hDisp[i], labels[i].c_str(), "l");
  }
  leg->Draw();
  vector<int> colorsVec(colors, colors + 3);
  drawMeans(hIn, colorsVec, "x_{J#gamma}", .5, .57);

  drawer d;
  d.drawAll(samples, extraFeatures, .18, .85, 16, gPad->GetWh() * 0.8);

  c->SaveAs(pdfPath);
}

// --- Part 2: derive a Data/MC emfrac reweighting and see what it does to MC's x_J ---

const int nFineEmfracBins = 10; // finer than ana::emfracBins - just for deriving w(emfrac)

TH1D * fillFineEmfracHist(const vector<EmfracEvent> & events, const char * name) {
  TH1D * h = new TH1D(name, ";Jet emfrac;Counts", nFineEmfracBins, 0, 1);
  for (const EmfracEvent & ev : events) h->Fill(ev.emfrac, ev.weight);
  return h;
}

// Page 3: Data's (purity-corrected) and MC's own shape-normalized emfrac distributions
// (top) and their per-bin ratio (bottom) - the ratio IS the reweighting factor w(emfrac)
// applied in Part 2 below, drawn here so it can be inspected on its own before trusting
// what it does to x_J. dataEventsInBin must already carry both Region A and Region C
// events (collectEvents' includeRegionC=true) in the single pT bin Parts 1-2 focus on.
TH1D * drawEmfracRatioPage(TCanvas * c, const char * pdfPath,
    const vector<EmfracEvent> & dataEventsInBin, const vector<EmfracEvent> & mcEventsInBin,
    const vector<string> & extraFeatures) {
  vector<EmfracEvent> dataA, dataC;
  splitByRegion(dataEventsInBin, dataA, dataC);
  TH1D * hDataRawA = fillFineEmfracHist(dataA, "hEmfracDataRegionA");
  TH1D * hDataRawC = fillFineEmfracHist(dataC, "hEmfracDataRegionC");
  TH1D * hDataRaw = purityCorrectHist(hDataRawA, hDataRawC, singleBinPtLo, singleBinPtHi, "hEmfracDataRaw");
  TH1D * hMCRaw = fillFineEmfracHist(mcEventsInBin, "hEmfracMCRaw");

  TH1D * hData = (TH1D*)hDataRaw->Clone("hEmfracDataNorm");
  TH1D * hMC = (TH1D*)hMCRaw->Clone("hEmfracMCNorm");
  if (hData->Integral() > 0) hData->Scale(1. / hData->Integral());
  if (hMC->Integral() > 0) hMC->Scale(1. / hMC->Integral());

  // w(emfrac) = normalized_Data / normalized_MC, per the coordinator's exact
  // specification - the per-event weight applied in Part 2 below. Where MC's normalized
  // bin is ~0 (no MC events at that emfrac to reweight in the first place), fall back to
  // w=1 (no reweight) rather than dividing by ~0 - flagged with a warning rather than
  // silently producing an enormous or infinite weight from noise.
  TH1D * hRatio = (TH1D*)hData->Clone("hEmfracRatio");
  hRatio->Reset("ICES");
  for (int b = 1; b <= nFineEmfracBins; b++) {
    double d = hData->GetBinContent(b);
    double m = hMC->GetBinContent(b);
    if (m > 1e-8) {
      hRatio->SetBinContent(b, d / m);
    } else {
      hRatio->SetBinContent(b, 1.0);
      if (d > 1e-8)
        cout << "    WARNING: emfrac bin " << b << " [" << hRatio->GetBinLowEdge(b) << ","
             << hRatio->GetBinLowEdge(b) + hRatio->GetBinWidth(b)
             << ") has Data content but ~0 MC content - falling back to w=1 (no reweight) there." << endl;
    }
  }

  cout << "  emfrac reweighting factors (Data/MC, normalized):" << endl;
  for (int b = 1; b <= nFineEmfracBins; b++) {
    cout << "    [" << hRatio->GetBinLowEdge(b) << ","
         << hRatio->GetBinLowEdge(b) + hRatio->GetBinWidth(b) << "): w=" << hRatio->GetBinContent(b)
         << "  (Data raw=" << hDataRaw->GetBinContent(b) << ", MC raw=" << hMCRaw->GetBinContent(b) << ")"
         << endl;
  }

  c->Clear();
  c->cd();
  TPad * p1 = new TPad("p1_ratio", "", 0, .35, 1, 1);
  TPad * p2 = new TPad("p2_ratio", "", 0, 0, 1, .35);
  p1->Draw();
  p2->Draw();

  p1->cd();
  p1->SetBottomMargin(0.02);
  p1->SetLeftMargin(.15);
  gPad->SetTicks(1, 1);
  double ymax = std::max(hData->GetMaximum(), hMC->GetMaximum());
  hData->SetLineColor(kBlack);
  hData->SetMarkerColor(kBlack);
  hData->SetMarkerStyle(20);
  hData->SetLineWidth(2);
  hData->GetXaxis()->SetLabelSize(0);
  hData->GetXaxis()->SetTitle("");
  hData->GetYaxis()->SetTitle("Normalized counts");
  hData->GetYaxis()->SetRangeUser(0, ymax * 1.4);
  hData->Draw("p e");
  hMC->SetLineColor(kRed + 1);
  hMC->SetMarkerColor(kRed + 1);
  hMC->SetMarkerStyle(21);
  hMC->SetLineWidth(2);
  hMC->Draw("p e same");
  TLegend * leg = new TLegend(.18, .65, .5, .85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hData, "p+p Run24 Data (purity-corr.)", "lep");
  leg->AddEntry(hMC, "Pythia8 #gamma+jet MC", "lep");
  leg->Draw();
  vector<TH1D*> hMeanSrc = {hDataRaw, hMCRaw};
  vector<int> meanColors = {kBlack, kRed + 1};
  drawMeans(hMeanSrc, meanColors, "emfrac", .18, .6);
  drawer d;
  d.drawAll({}, extraFeatures, .5, .85, 16, gPad->GetWh() * 0.8);

  p2->cd();
  p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.3);
  p2->SetLeftMargin(.15);
  gPad->SetTicks(1, 1);
  hRatio->SetLineColor(kBlack);
  hRatio->SetMarkerColor(kBlack);
  hRatio->SetMarkerStyle(20);
  hRatio->GetYaxis()->SetRangeUser(0, std::max(2.0, hRatio->GetMaximum() * 1.2));
  hRatio->GetYaxis()->SetTitle("Data / MC");
  hRatio->GetYaxis()->SetTitleSize(0.09);
  hRatio->GetYaxis()->SetTitleOffset(0.7);
  hRatio->GetYaxis()->SetLabelSize(0.08);
  hRatio->GetXaxis()->SetTitle("Jet emfrac");
  hRatio->GetXaxis()->SetTitleSize(0.09);
  hRatio->GetXaxis()->SetLabelSize(0.08);
  hRatio->Draw("p e");
  TLine * line = new TLine(0, 1, 1, 1);
  line->SetLineStyle(9);
  line->Draw("same");

  c->SaveAs(pdfPath);
  return hRatio;
}

// Page 4: MC's x_J shape before vs. after applying w(emfrac) as an extra per-event
// weight on top of its existing mcWeight*crossSectionScale, with purity-corrected Data
// alongside as the reference the reweighting is meant to approach.
void drawReweightedXjPage(TCanvas * c, const char * pdfPath, TH1D * hRatio,
    const vector<EmfracEvent> & dataEventsInBin, const vector<EmfracEvent> & mcEventsInBin,
    const vector<string> & extraFeatures) {
  vector<EmfracEvent> dataA, dataC;
  splitByRegion(dataEventsInBin, dataA, dataC);
  TH1D * hDataA = new TH1D("hXjDataA", "", ana::nUnfoldXjBins, ana::unfoldXjBins);
  TH1D * hDataC = new TH1D("hXjDataC", "", ana::nUnfoldXjBins, ana::unfoldXjBins);
  for (const EmfracEvent & ev : dataA) hDataA->Fill(ev.xj, ev.weight);
  for (const EmfracEvent & ev : dataC) hDataC->Fill(ev.xj, ev.weight);
  TH1D * hData = purityCorrectHist(hDataA, hDataC, singleBinPtLo, singleBinPtHi, "hXjData");
  delete hDataA;
  delete hDataC;

  TH1D * hMCBefore = new TH1D("hXjMCBefore", "", ana::nUnfoldXjBins, ana::unfoldXjBins);
  TH1D * hMCAfter = new TH1D("hXjMCAfter", "", ana::nUnfoldXjBins, ana::unfoldXjBins);
  for (const EmfracEvent & ev : mcEventsInBin) {
    hMCBefore->Fill(ev.xj, ev.weight);
    double w = hRatio->GetBinContent(hRatio->FindBin(ev.emfrac));
    hMCAfter->Fill(ev.xj, ev.weight * w);
  }

  vector<TH1D*> hRaw = {hData, hMCBefore, hMCAfter};
  vector<string> labels = {"p+p Run24 Data (purity-corr.)", "Pythia8 MC, before emfrac reweight", "Pythia8 MC, after emfrac reweight"};
  vector<int> colors = {kBlack, kRed + 1, kAzure + 2};
  vector<int> markers = {20, 21, 22};

  vector<TH1D*> hDisp(3);
  double ymax = 0;
  for (int i = 0; i < 3; i++) {
    cout << "    " << labels[i] << ": " << hRaw[i]->GetEntries() << " entries, "
         << hRaw[i]->Integral() << " weighted" << endl;
    hDisp[i] = unfold_utility::densityForDisplay(hRaw[i], Form("%s_disp", hRaw[i]->GetName()));
    if (hDisp[i]->Integral() > 0) hDisp[i]->Scale(1. / hDisp[i]->Integral());
    for (int b = 0; b <= hDisp[i]->GetNbinsX() + 1; b++) {
      double v = hDisp[i]->GetBinContent(b);
      if (std::isnan(v) || std::isinf(v)) hDisp[i]->SetBinContent(b, 0);
    }
    ymax = std::max(ymax, hDisp[i]->GetMaximum());
  }

  c->Clear();
  c->cd();
  gPad->SetTicks(1, 1);
  gPad->SetLeftMargin(.15);
  gPad->SetBottomMargin(.13);

  TLegend * leg = new TLegend(.45, .6, .85, .85);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  for (int i = 0; i < 3; i++) {
    hDisp[i]->SetLineColor(colors[i]);
    hDisp[i]->SetMarkerColor(colors[i]);
    hDisp[i]->SetMarkerStyle(markers[i]);
    hDisp[i]->SetLineWidth(2);
    hDisp[i]->GetXaxis()->SetTitle("x_{J#gamma}");
    hDisp[i]->GetYaxis()->SetTitle("(1/N) dN/dx_{J#gamma}");
    hDisp[i]->GetYaxis()->SetRangeUser(0, ymax * 1.4);
    hDisp[i]->Draw(i == 0 ? "p e" : "p e same");
    leg->AddEntry(hDisp[i], labels[i].c_str(), "lep");
  }
  leg->Draw();
  drawMeans(hRaw, colors, "x_{J#gamma}", .45, .57, 0.045, 0.40);

  drawer d;
  d.drawAll({}, extraFeatures, .18, .85, 16, gPad->GetWh() * 0.8);

  c->SaveAs(pdfPath);
}

// --- Part 3: <x_J> Data(purity-corr.)/MC ratio vs. photon p_T, by emfrac category ---

// Exact unbinned weighted mean/error accumulator - same as
// temporary_study/draw_xj_data_purity.C's identical Accum. Needed for the same reason:
// purity-correcting the sums directly (rather than a filled-then-corrected histogram's
// TH1::GetMean()) gives the exact, not a binned-approximation, mean.
struct Accum {
  double sumw = 0, sumw2 = 0, sumwx = 0, sumwx2 = 0;
  Long64_t n = 0;
  void fill(double x, double w) {
    sumw += w; sumw2 += w * w; sumwx += w * x; sumwx2 += w * x * x; n++;
  }
  double mean() const { return sumw > 0 ? sumwx / sumw : 0; }
  double meanErr() const {
    if (sumw <= 0 || n < 2) return 0;
    double var = sumwx2 / sumw - mean() * mean();
    if (var < 0) var = 0;
    double neff = sumw2 > 0 ? sumw * sumw / sumw2 : n;
    return std::sqrt(var / neff);
  }
};

const int nEmfracCategories = 4; // ana::emfracBins' 3 slices + "Total" (all emfrac), index 3

// Page 5: <x_J> Data(purity-corrected)/MC ratio vs. photon p_T, one curve per
// ana::emfracBins slice plus one for the inclusive "Total" sample, across every
// ana::ptBinsUsed bin - unlike Parts 1-2's single 15-20 GeV bin. Data's mean is the
// two-purity-corrected exact unbinned mean (purityCorrectCoeffs applied directly to the
// Accum sums, mirroring temporary_study/draw_xj_data_purity.C's accCorr construction
// exactly - not a TH1::GetMean() on a purity-corrected histogram); MC stays raw Region A
// signal, same as everywhere else in this file. Points with no valid purity-corrected
// Data mean (Region C empty and Region A also empty) or no MC statistics in that
// (p_T,emfrac) slice are simply omitted rather than plotted as a misleading 0.
void drawMeanXjVsPtPage(TCanvas * c, const char * pdfPath, const vector<EmfracEvent> & dataEvents,
    const vector<EmfracEvent> & mcEvents, const vector<string> & extraFeatures) {
  vector<vector<Accum>> accDataA(ana::nPtBinsUsed, vector<Accum>(nEmfracCategories));
  vector<vector<Accum>> accDataC(ana::nPtBinsUsed, vector<Accum>(nEmfracCategories));
  vector<vector<Accum>> accMC(ana::nPtBinsUsed, vector<Accum>(nEmfracCategories));

  auto fillAcc = [&](const vector<EmfracEvent> & events, vector<vector<Accum>> & acc, bool wantRegionC) {
    for (const EmfracEvent & ev : events) {
      if ((ev.abcdRegion == 2) != wantRegionC) continue;
      int ipt = ana::findPtBin(ev.phopt);
      int iused = ipt - ana::firstUsedPtBin;
      if (iused < 0 || iused >= ana::nPtBinsUsed) continue;
      int iemfrac = ana::findEmfracBin(ev.emfrac);
      if (iemfrac < 0) continue;
      acc[iused][iemfrac].fill(ev.xj, ev.weight);
      acc[iused][3].fill(ev.xj, ev.weight); // Total (all emfrac)
    }
  };
  fillAcc(dataEvents, accDataA, false);
  fillAcc(dataEvents, accDataC, true);
  fillAcc(mcEvents, accMC, false);

  const int colors[nEmfracCategories] = {kBlue + 1, kGreen + 2, kRed + 1, kBlack};
  const int markers[nEmfracCategories] = {20, 21, 22, 23};
  const vector<string> labels = {
      Form("%.1f < emfrac < %.1f", ana::emfracBins[0], ana::emfracBins[1]),
      Form("%.1f < emfrac < %.1f", ana::emfracBins[1], ana::emfracBins[2]),
      Form("%.1f < emfrac < %.1f", ana::emfracBins[2], ana::emfracBins[3]),
      "Total (all emfrac)",
  };

  vector<TGraphErrors*> g(nEmfracCategories);
  for (int cat = 0; cat < nEmfracCategories; cat++) {
    g[cat] = new TGraphErrors();
    g[cat]->SetName(Form("gMeanXjRatio_%d", cat));
    int ipoint = 0;
    for (int ip = 0; ip < ana::nPtBinsUsed; ip++) {
      double ptlo = ana::ptBinsUsed[ip], pthi = ana::ptBinsUsed[ip + 1];
      // Small per-category horizontal offset so the 4 series don't sit exactly on top of
      // each other at each pT bin's center.
      double ptx = 0.5 * (ptlo + pthi) + (cat - 1.5) * 0.06 * (pthi - ptlo);

      float pA = ana::getPurity(ptlo, pthi, "nominal", ir);
      float pC = ana::getPurityC(ptlo, pthi, "nominal", ir);
      double NA = accDataA[ip][cat].sumw, NC = accDataC[ip][cat].sumw;
      float coeffA, coeffC;
      unfold_utility::purityCorrectCoeffs(pA, pC, NA, NC, coeffA, coeffC);
      double Ncorr = coeffA * NA - coeffC * NC;
      if (Ncorr <= 0) continue;

      // Mirrors temporary_study/draw_xj_data_purity.C's accCorr construction exactly -
      // effective-N == Ncorr, not a separate sum-of-squared-weights.
      Accum accCorr;
      accCorr.sumw   = Ncorr;
      accCorr.sumw2  = Ncorr;
      accCorr.sumwx  = coeffA * accDataA[ip][cat].sumwx  - coeffC * accDataC[ip][cat].sumwx;
      accCorr.sumwx2 = coeffA * accDataA[ip][cat].sumwx2 - coeffC * accDataC[ip][cat].sumwx2;
      accCorr.n      = accDataA[ip][cat].n + accDataC[ip][cat].n;
      double meanDataCorr = accCorr.mean();
      double meanDataCorrErr = accCorr.meanErr();
      if (meanDataCorr == 0) continue;

      if (accMC[ip][cat].sumw <= 0) continue;
      double meanMC = accMC[ip][cat].mean();
      double meanMCErr = accMC[ip][cat].meanErr();
      if (meanMC == 0) continue;

      double ratio = meanDataCorr / meanMC;
      double relErrData = meanDataCorrErr / meanDataCorr;
      double relErrMC = meanMCErr / meanMC;
      double ratioErr = fabs(ratio) * std::sqrt(relErrData * relErrData + relErrMC * relErrMC);

      g[cat]->SetPoint(ipoint, ptx, ratio);
      g[cat]->SetPointError(ipoint, 0, ratioErr);
      ipoint++;
    }
  }

  c->Clear();
  c->cd();
  gPad->SetTicks(1, 1);
  gPad->SetLeftMargin(.15);
  gPad->SetBottomMargin(.13);

  TH1D * hFrame = new TH1D("hFrameMeanXjRatio",
      ";p_{T}^{#gamma} [GeV];#LTx_{J#gamma}#GT Data(purity-corr.) / MC",
      ana::nPtBinsUsed, ana::ptBinsUsed);
  hFrame->SetStats(0);
  hFrame->GetYaxis()->SetRangeUser(0.84, 1.3);
  hFrame->Draw("axis");

  TLine * line = new TLine(ana::ptBinsUsed[0], 1, ana::ptBinsUsed[ana::nPtBinsUsed], 1);
  line->SetLineStyle(9);
  line->Draw("same");

  TLegend * leg = new TLegend(.5, .65, .85, .88);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  for (int cat = 0; cat < nEmfracCategories; cat++) {
    g[cat]->SetLineColor(colors[cat]);
    g[cat]->SetMarkerColor(colors[cat]);
    g[cat]->SetMarkerStyle(markers[cat]);
    g[cat]->SetLineWidth(2);
    g[cat]->Draw("p same");
    leg->AddEntry(g[cat], labels[cat].c_str(), "lep");
  }
  leg->Draw();

  drawer d;
  d.drawAll({}, extraFeatures, .18, .85, 16, gPad->GetWh() * 0.8);

  c->SaveAs(pdfPath);
}

void draw_emfrac_xj() {
  gStyle->SetOptStat(0);
  TH1::AddDirectory(kFALSE);
  const char * pdfPath = ana::path("pdfs/draw_emfrac_xj.pdf");

  vector<string> commonFeatures = {
      Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", singleBinPtLo, singleBinPtHi),
      Form("Jet R=%.1f", ana::JetRs[ir]),
  };
  vector<string> dataFeatures = commonFeatures;
  dataFeatures.push_back("Data: purity-corrected (two-purity method)");
  vector<string> mcFeatures = commonFeatures;
  mcFeatures.push_back("MC: Region A only (signal, no background)");
  vector<string> ratioFeatures = {
      Form("Jet R=%.1f", ana::JetRs[ir]),
      "Data: purity-corrected (two-purity method)",
  };

  TCanvas * c = new TCanvas("c", "", 700, 700);
  c->SaveAs(Form("%s[", pdfPath));

  cout << "Collecting Data (Region A + Region C, full p_{T} range):" << endl;
  vector<EmfracEvent> dataEvents = collectEvents("Data", "pythia", false, 1.0, true);
  cout << "Collecting MC Photon (Photon5+10+20, pythia, Region A, full p_{T} range):" << endl;
  vector<EmfracEvent> mcEvents = collectMCPhotonEvents();

  vector<EmfracEvent> dataEventsInBin = filterPtRange(dataEvents, singleBinPtLo, singleBinPtHi);
  vector<EmfracEvent> mcEventsInBin = filterPtRange(mcEvents, singleBinPtLo, singleBinPtHi);

  // --- Part 1: x_J split by the coarse ana::emfracBins ---
  cout << "Page 1 (Data, coarse emfrac split, purity-corrected):" << endl;
  vector<TH1D*> hData = binByEmfracCoarsePurityCorrected(dataEventsInBin, "data");
  drawEmfracPage(c, pdfPath, hData, {"p+p Run24 Data (purity-corr.)"}, dataFeatures);

  cout << "Page 2 (MC, coarse emfrac split):" << endl;
  vector<TH1D*> hMC = binByEmfracCoarse(mcEventsInBin, "mc");
  drawEmfracPage(c, pdfPath, hMC, {"Pythia8 #gamma+jet MC"}, mcFeatures);

  // --- Part 2: derive w(emfrac) = normalized Data / normalized MC, reweight MC, compare ---
  cout << "Page 3 (Data vs MC emfrac distributions + ratio = reweighting factor):" << endl;
  TH1D * hRatio = drawEmfracRatioPage(c, pdfPath, dataEventsInBin, mcEventsInBin, commonFeatures);

  cout << "Page 4 (MC x_J before vs after emfrac reweighting):" << endl;
  drawReweightedXjPage(c, pdfPath, hRatio, dataEventsInBin, mcEventsInBin, commonFeatures);

  // --- Part 3: <x_J> Data/MC ratio vs pT, by emfrac category ---
  cout << "Page 5 (<x_J> Data(purity-corr.)/MC ratio vs p_T, by emfrac category):" << endl;
  drawMeanXjVsPtPage(c, pdfPath, dataEvents, mcEvents, ratioFeatures);

  c->SaveAs(Form("%s]", pdfPath));
  cout << "Wrote " << pdfPath << endl;
}
