#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
// The original gammajet project builds its OWN, differently-laid-out drawer/ana classes
// into /home/samson72/root/lib/libgammajet.so, sitting on the same library search path
// as this project's libgammajet_unfold.so. Without forcing which one loads first, ROOT's
// implicit symbol autoload can bind drawer/ana calls to the wrong (mismatched-layout)
// library and segfault - explicit load removes the ambiguity.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Systematic-uncertainty comparison: for each systematic variation, unfold Data's own
// purity-corrected xJ spectrum through THAT variation's own response + purity curve
// (exactly the draw_purity_corrected.C Data pipeline, just re-run once per systag via
// drawer's systag parameter - see src/unfolder.h's constructor comment for what each
// variation changes), and compare to the identically-built nominal result. "herwig" is
// the one exception: it holds Data/purity/dataSystag fixed at nominal and instead swaps
// which MC generator built the response matrix (drawer's sim parameter), to probe the
// unfolding's sensitivity to MC modeling rather than to a reco-level selection variant.
//
// Also includes the unfolding regularization (niter) choice as a systematic: niterLow/
// niterHigh reuse the exact same nominal response and purity-corrected Data as the
// nominal result (dataSystag="nominal" in systSources below) - the ONLY thing that
// changes is which niter unfoldOnce() is called with (niterate-1/niterate+1 instead of
// niterate). Unlike the other seven, this needs no separate production reprocessing.
//
// "priorSensitivity" is a different kind of source again - it doesn't reprocess Data
// through any different production sample or niter at all. It's read directly from
// hists/prior_sensitivity_nominal.root, written by drawing/draw_prior_sensitivity.C -
// run that macro (with its default "nominal" argument) before this one. See that file's
// header for the full method (data-informed prior reweighting, following the sPHENIX
// PPG08 dijet-xJ note and ATLAS's photon-jet xJ paper); only its "w"-variant fracDiff
// (hw_pt<N>, re-shape-normalized here the same way every other source is) is used, as
// discussed there - sqrt(w)/w^1.5 are diagnostic only in that file and aren't read here.
//
// Total uncertainty: JER/JES/emscale/EMR/niter are true two-point (high/low)
// systematics - each of their ten sources feeds ONLY the up or down total per bin,
// whichever matches its own sign that bin (see asymmetricSystematics below), rather
// than being symmetrized. The five ABCD sideband-boundary sources (narrowBDT/narrowISO/
// narrowBDTbkg/narrowISObkg/wideISObkg - see ana.h's isoBins/isoBinsHigh/bdtGoodLow/
// bdtBadLow comment) are first combined into their own symmetric purityUnc quadrature sum
// (purityMembers below - every one of the five is treated as an independent symmetrized
// source, including narrowISObkg/wideISObkg even though they're a genuine two-sided
// variation of the same boundary, NOT sign-split against each other) and enters the grand
// total as ONE more symmetric source, rather than five independent ones - this follows
// PPG12's treatment of the ABCD boundary variations as a single combined "Purity"
// systematic (sPHENIX isolated-photon analysis note, Sec. 5.3) rather than independently
// summing each cut shift. Every remaining source (threejet/herwig/priorSensitivity) is
// symmetrized the same way: its full magnitude feeds both the up and down total.
//
// Per (systematic, pT bin): a page with nominal vs variation overlaid (top) and their
// ratio (bottom) - both the ratio and the fractional difference (ratio-1) are also
// written to systematics.root for reuse as the actual systematic uncertainty numbers.
// This still runs (and is written to systematics.root) individually for each of the five
// purityMembers, even though they're combined into one Purity source below - useful for
// debugging any one boundary's own effect.
// Final page: fractional difference vs xJ, one panel per pT bin, all variations overlaid,
// so the relative size of each systematic is visible at a glance - grouped per
// displayGroups: each high/low pair shares one color/legend entry (both member curves
// drawn as-is), while each symmetric source is drawn twice (its curve and that curve's
// negation, same color) so the page visually matches the +/- treatment the total gives it.
// The five purityMembers are excluded from displayGroups entirely - the combined
// purityUnc curve gets its own single "Purity" entry on this page instead, drawn the same
// curve-plus-negation way as any other symmetric source.
//
// Both nominal and each variation are shape-normalized (unit area) before the ratio is
// taken, so a systematic that shifts the total accepted Data event count (narrowBDT/
// narrowISO/threejet all change which events pass selection) doesn't masquerade as a
// shape difference in xJ - only genuine shape effects survive into the ratio/fracdiff.

// Jet radius index - mutable (not const) so draw_systematics(int) can set it at the top
// of the function, before any of the code below (all written against this global) runs.
// Defaults to the nominal R=0.4 working point used throughout the note.
int ir = 2;
const int nPtBinsUsed = ana::nPtBinsUsed; // physics analysis only uses ana::ptBins[ana::firstUsedPtBin..]
const int niterate = 2; // matches draw_purity_corrected.C's best-iteration scan result
// Derived from ana::systags (src/ana.h) - the definitive systag reprocessing list, minus
// "nominal" (the baseline every source here is compared against, not a source itself) -
// plus the four sources that aren't systag reprocessings at all (herwig: different
// generator sample; niterLow/niterHigh: different unfolding iteration count;
// priorSensitivity: different unfolding prior - none of these have an insitu_tree
// equivalent, see ana.h's comment on ana::systags). A systag added to ana::systags
// propagates here automatically; it still needs its own systColors/systSources entry
// below (that .at() lookup throws loudly, rather than silently dropping it, if missing).
const vector<string> systematics = [] {
  vector<string> v;
  for (const string & s : ana::systags) if (s != "nominal") v.push_back(s);
  v.insert(v.end(), {"herwig", "niterLow", "niterHigh", "priorSensitivity"});
  return v;
}();
const map<string,int> systColors = {
  {"JERhigh",      kRed+1},
  {"JERlow",       kOrange+7},
  {"emscale_high", kAzure+1},
  {"emscale_low",  kTeal+2},
  {"EMRhigh",      kGray+2},
  {"EMRlow",       kGray+1},
  {"jes_high",     kPink+6},
  {"jes_low",      kYellow+2},
  {"threejet",     kGreen+2},
  {"narrowBDT",    kBlue+1},
  {"narrowISO",    kMagenta+1},
  {"narrowBDTbkg", kBlue+3},
  {"narrowISObkg", kMagenta+3},
  {"wideISObkg",   kMagenta-4},
  {"herwig",       kSpring+5},
  {"niterLow",     kCyan+2},
  {"niterHigh",    kViolet+1},
  {"priorSensitivity", kOrange+2},
};

// Per-source config: which production reprocessing (sim/dataSystag - see drawer's sim/
// systag params) and which unfoldOnce() iteration count each systematic uses. The first
// nine vary the reco-level production (their own dataSystag) against the pythia
// response, holding niter fixed at niterate; "herwig" instead holds dataSystag fixed at
// "nominal" and swaps the response-matrix generator (sim) to test the unfolding's
// sensitivity to MC modeling; the two niter sources hold the production fixed at
// "nominal"/pythia and vary niter instead. "priorSensitivity" has NO entry here - it
// doesn't reprocess Data through anything, so the main loop below special-cases it and
// never looks this map up for that key.
struct SystSource { string sim; string dataSystag; int niter; };
const map<string, SystSource> systSources = {
  {"JERhigh",      {"pythia", "JERhigh",      niterate}},
  {"JERlow",       {"pythia", "JERlow",       niterate}},
  {"emscale_high", {"pythia", "emscale_high", niterate}},
  {"emscale_low",  {"pythia", "emscale_low",  niterate}},
  {"EMRhigh",      {"pythia", "EMRhigh",      niterate}},
  {"EMRlow",       {"pythia", "EMRlow",       niterate}},
  {"jes_high",     {"pythia", "jes_high",     niterate}},
  {"jes_low",      {"pythia", "jes_low",      niterate}},
  {"threejet",     {"pythia", "threejet",     niterate}},
  {"narrowBDT",    {"pythia", "narrowBDT",    niterate}},
  {"narrowISO",    {"pythia", "narrowISO",    niterate}},
  {"narrowBDTbkg", {"pythia", "narrowBDTbkg", niterate}},
  {"narrowISObkg", {"pythia", "narrowISObkg", niterate}},
  {"wideISObkg",   {"pythia", "wideISObkg",   niterate}},
  {"herwig",       {"herwig", "nominal",      niterate}},
  {"niterLow",     {"pythia", "nominal",      niterate-1}},
  {"niterHigh",    {"pythia", "nominal",      niterate+1}},
};

// Two-point (high/low) systematic sources contributed asymmetrically to the total: per
// bin, each one's own signed fracDiff value feeds the "up" quadrature sum if positive or
// "down" if negative that bin (not paired/enveloped with its high/low counterpart - each
// of the eight contributes independently on its own sign). Every systag NOT listed here
// is symmetrized instead (see the total systematic uncertainty section below).
// Derived from ana::asymmetricSystagPairs (src/ana.h) flattened to individual names,
// plus niterHigh/niterLow (not a systags pair - see the systematics comment above) -
// a pair added to ana::asymmetricSystagPairs propagates here automatically.
const set<string> asymmetricSystematics = [] {
  set<string> s;
  for (const auto & pr : ana::asymmetricSystagPairs) { s.insert(pr.first); s.insert(pr.second); }
  s.insert("niterHigh"); s.insert("niterLow");
  return s;
}();

// Categorical palette for the 8 display groups below (plus the separately-drawn Purity
// curve - see purityMembers): a validated 9-hue, colorblind-safe ordering (fixed order,
// never cycled/reassigned) - each hex registered once as a ROOT color index via
// TColor::GetColor(). Chosen over plain kXXX constants because several of those
// (kAzure/kCyan/kTeal, kMagenta/kPink, kGreen/kSpring) sit too close in hue to reliably
// tell apart across several overlaid curves.
const int colorBlue    = TColor::GetColor("#2a78d6");
const int colorOrange  = TColor::GetColor("#eb6834");
const int colorAqua    = TColor::GetColor("#1baf7a");
const int colorYellow  = TColor::GetColor("#eda100");
const int colorMagenta = TColor::GetColor("#e87ba4");
const int colorPurity  = TColor::GetColor("#008300");
const int colorRed     = TColor::GetColor("#e34948");
const int colorBrown   = TColor::GetColor("#8c564b");
const int colorGrey    = TColor::GetColor("#767676");

// Display grouping for the final overlay page only (the per-(systag,pT) comparison pages
// above still use systColors, one distinct color per individual systag). Each
// asymmetricSystematics high/low pair shares one color and one legend entry, labeled by
// the group name rather than each half's systag - they're two views of one physical
// source, and both member curves are drawn as-is (their actual signed value, matching how
// the total treats them). Every symmetric source keeps its own color and legend entry,
// but is drawn TWICE - its real fracDiff curve and that curve's negation - so the page
// visually shows the same +/- treatment the total's quadrature sum already applies to it
// (a single symmetric source only measures one sign of deviation, but contributes
// symmetrically to the total).
//
// purityMembers (narrowBDT/narrowISO/narrowBDTbkg/narrowISObkg/wideISObkg) are
// deliberately NOT listed here - they're combined into one symmetric purityUnc curve
// below and drawn as their own single "Purity" entry further down (curve plus negation,
// same as any other symmetric source here), rather than getting individual DisplayGroup
// entries (which would show five sub-source curves instead of the one combined
// systematic they actually feed into the total as).
struct DisplayGroup { string label; int color; vector<string> members; bool symmetric; };
const vector<DisplayGroup> displayGroups = {
  {"JER",       colorBlue,    {"JERhigh", "JERlow"},          false},
  {"JES",       colorOrange,  {"jes_high", "jes_low"},        false},
  {"EM scale",  colorAqua,    {"emscale_high", "emscale_low"},false},
  {"EMR",       colorGrey,    {"EMRhigh", "EMRlow"},          false},
  {"niter",     colorYellow,  {"niterHigh", "niterLow"},      false},
  {"threejet",  colorMagenta, {"threejet"},                   true},
  {"herwig",    colorRed,     {"herwig"},                     true},
  {"prior sens.", colorBrown, {"priorSensitivity"},           true},
};
// The five ABCD sideband-boundary systematics combined into one "Purity" systematic
// before entering the grand total - see the header comment and the total-uncertainty
// section below. All five are treated as independent symmetrized sources here (none of
// them are in ana::asymmetricSystagPairs, so asymmetricSystematics below doesn't include
// them either) - narrowISObkg/wideISObkg are a genuine two-sided variation of the same
// isolation-gap boundary, but are deliberately NOT sign-split against each other like
// JER/emscale/jes/EMR are; each contributes its own full magnitude to both up and down.
const set<string> purityMembers = {
  "narrowBDT", "narrowISO", "narrowBDTbkg", "narrowISObkg", "wideISObkg"
};
// Set inside draw_systematics(int) from ir - the nominal R=0.4 default reproduces the
// unsuffixed filenames every other macro/main.tex reads; every other radius gets its own
// _<rname>-suffixed pair instead of clobbering the nominal file.
string pdfPathStr, rootPathStr;
const char * pdfPath;
const char * rootPath;

// densityForDisplay now lives in unfold_utility - see src/unfold_utility.h.
// buildFullyCorrected now lives in unfold_utility (purity-corrects all ana::nPtBins
// slices via unfold_utility::purityCorrect and reflattens for RooUnfold) - see
// src/unfold_utility.h.

// Data's purity-corrected xJ spectrum, unfolded through systag's own response - the same
// pipeline draw_purity_corrected.C uses for draw_one_sample(0, ...) ("data"), just
// factored out so it can be called once per systematic (including "nominal" itself).
// niter is explicit (rather than always the global niterate) so the niterLow/niterHigh
// regularization sources can reuse systag="nominal" and vary only the iteration count.
// sim selects which MC built the response matrix (pythia for every source except
// "herwig" - see systSources); purity always comes from dataSystag's own purity file
// (puritymaker.C is pythia-only), regardless of sim.
TH1D * getUnfoldedData(string sim, string systag, int niter, const char * name) {
  drawer d(sim, systag);
  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i",ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i",ir), 1);
  TH1D * flatA = d.get(Form("hrecoxj%i_0",ir), 0);
  TH1D * flatC = d.get(Form("hrecoxj%i_2",ir), 0);
  TH1D * flatCorrected = unfold_utility::buildFullyCorrected(flatA, flatC, systag.c_str(), systag);
  return unfold_utility::unfoldOnce(respRecoTemplate, respTruthTemplate, respMatrix2D, flatCorrected, niter, name);
}

void draw_systematics(int jetRadiusIndex = 2) {
  gStyle->SetOptStat(0);

  ir = jetRadiusIndex;
  pdfPathStr  = (ir == 2) ? "/home/samson72/sphnx/gammajet_unfold/pdfs/draw_systematics.pdf"
                          : Form("/home/samson72/sphnx/gammajet_unfold/pdfs/draw_systematics_%s.pdf", ana::rnames[ir]);
  rootPathStr = (ir == 2) ? "/home/samson72/sphnx/gammajet_unfold/hists/systematics.root"
                          : Form("/home/samson72/sphnx/gammajet_unfold/hists/systematics_%s.root", ana::rnames[ir]);
  pdfPath  = pdfPathStr.c_str();
  rootPath = rootPathStr.c_str();

  // drawText/drawAll don't touch any per-instance file data - one generic instance
  // (default sim/systag - irrelevant here) is reused purely for label drawing, since the
  // per-systag drawer built inside getUnfoldedData() goes out of scope with that call.
  drawer dLabel;

  TFile * fout = TFile::Open(rootPath, "RECREATE");
  TCanvas * c = new TCanvas("c","",700,900);
  c->SaveAs(Form("%s[", pdfPath));

  TH1D * flatNominal = getUnfoldedData("pythia", "nominal", niterate, "hUnfoldedNominal");
  vector<TH1D*> nominalDisp(ana::nPtBins);
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hNom = unfold_utility::unflattenXj(flatNominal, ipt, Form("hNominal_pt%d", ipt));
    nominalDisp[ipt] = unfold_utility::densityForDisplay(hNom, Form("hNominalDisp_pt%d", ipt));
    // Shape-normalize before comparing - a systematic that shifts the total accepted
    // event count (e.g. narrowBDT/narrowISO/threejet change which Data events pass
    // selection) shouldn't masquerade as a shape difference in xJ.
    nominalDisp[ipt]->Scale(1./nominalDisp[ipt]->Integral());
    nominalDisp[ipt]->GetYaxis()->SetTitle("Shape-normalized counts / bin width");
  }

  // fracDiff[systag][ipt]: (variation-nominal)/nominal vs xJ, kept for the final summary page.
  map<string, vector<TH1D*>> fracDiff;

  // priorSensitivity reads its per-pT-bin unfolded "w"-variant result from here rather
  // than reprocessing anything - see the header comment and drawing/draw_prior_sensitivity.C.
  TFile * fPriorSens = TFile::Open("/home/samson72/sphnx/gammajet_unfold/hists/prior_sensitivity_nominal.root");
  if (!fPriorSens || fPriorSens->IsZombie())
    cout << "WARNING: couldn't open prior_sensitivity_nominal.root - run "
         << "drawing/draw_prior_sensitivity.C(\"nominal\") first, or the priorSensitivity "
         << "source below will be empty." << endl;

  for (const string & systag : systematics) {
    bool isPriorSens = (systag == "priorSensitivity");
    const SystSource * src = isPriorSens ? nullptr : &systSources.at(systag);
    TH1D * flatVar = isPriorSens ? nullptr : getUnfoldedData(src->sim, src->dataSystag, src->niter, Form("hUnfolded_%s", systag.c_str()));
    fracDiff[systag] = vector<TH1D*>(ana::nPtBins);

    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      TH1D * hVar;
      if (isPriorSens) {
        TH1D * hRaw = fPriorSens ? (TH1D*)fPriorSens->Get(Form("hw_pt%d", ipt)) : nullptr;
        if (!hRaw) {
          cout << "WARNING: hw_pt" << ipt << " missing from prior_sensitivity_nominal.root - skipping pt" << ipt << "." << endl;
          continue;
        }
        hVar = (TH1D*)hRaw->Clone(Form("h%s_pt%d", systag.c_str(), ipt));
      } else {
        hVar = unfold_utility::unflattenXj(flatVar, ipt, Form("h%s_pt%d", systag.c_str(), ipt));
      }
      TH1D * hVarDisp = unfold_utility::densityForDisplay(hVar, Form("h%s_pt%d_disp", systag.c_str(), ipt));
      hVarDisp->Scale(1./hVarDisp->Integral());
      hVarDisp->GetYaxis()->SetTitle("Shape-normalized counts / bin width");

      TH1D * hratio = (TH1D*)hVarDisp->Clone(Form("hratio_%s_pt%d", systag.c_str(), ipt));
      hratio->Divide(nominalDisp[ipt]);
      // TH1::Divide already leaves a bin at 0 (not inf/NaN) whenever the denominator bin
      // is 0 - a bin where nominal or the variation is empty is physically meaningless to
      // compare, but a NaN-masking attempt here previously broke GetMaximum()/axis-range
      // calculations downstream and blanked the whole page. 0 is the simple, safe fallback.

      TH1D * hfrac = (TH1D*)hratio->Clone(Form("hfracdiff_%s_pt%d", systag.c_str(), ipt));
      for (int b = 1; b <= hfrac->GetNbinsX(); b++) {
        if (hratio->GetBinContent(b) != 0) hfrac->SetBinContent(b, hratio->GetBinContent(b) - 1);
      }
      hfrac->GetYaxis()->SetTitle("(Variation - Nominal) / Nominal");
      fracDiff[systag][ipt] = hfrac;

      c->Clear();
      c->cd();
      TPad * p1 = new TPad(Form("p1_%s_%d",systag.c_str(),ipt),"",0,.35,1,1);
      TPad * p2 = new TPad(Form("p2_%s_%d",systag.c_str(),ipt),"",0,0,1,.35);
      p1->Draw();
      p2->Draw();

      p1->cd();
      p1->SetBottomMargin(0.02);
      p1->SetLeftMargin(.15);
      gPad->SetTicks(1,1);
      nominalDisp[ipt]->SetLineColor(kBlack);
      nominalDisp[ipt]->SetMarkerColor(kBlack);
      nominalDisp[ipt]->SetMarkerStyle(20);
      nominalDisp[ipt]->SetLineWidth(2);
      nominalDisp[ipt]->GetXaxis()->SetLabelSize(0);
      nominalDisp[ipt]->GetXaxis()->SetTitle("");
      nominalDisp[ipt]->GetYaxis()->SetRangeUser(0, 0.25);
      nominalDisp[ipt]->Draw("p e");
      hVarDisp->SetLineColor(systColors.at(systag));
      hVarDisp->SetMarkerColor(systColors.at(systag));
      hVarDisp->SetMarkerStyle(21);
      hVarDisp->SetLineWidth(2);
      hVarDisp->Draw("p e same");
      TLegend * l = new TLegend(.55,.55,.85,.70);
      l->SetLineWidth(0);
      l->SetTextSize(0.035);
      l->AddEntry(nominalDisp[ipt], "Nominal");
      l->AddEntry(hVarDisp, systag.c_str());
      l->Draw();
      dLabel.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
          //isPriorSens ? Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV, nominal prior vs data-informed (w) prior, %d iter.", ana::JetRs[ir], ana::jet_calib_pt_cut[ir], niterate)
          //            : Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV, nominal %d iter. vs %s %d iter.", ana::JetRs[ir], ana::jet_calib_pt_cut[ir], niterate, systag.c_str(), src->niter)},
          },
          .5, .85, 16, gPad->GetWh()*0.8);

      p2->cd();
      p2->SetTopMargin(0.02);
      p2->SetBottomMargin(0.3);
      p2->SetLeftMargin(.15);
      gPad->SetTicks(1,1);
      hratio->SetLineColor(kBlack);
      hratio->SetMarkerColor(kBlack);
      hratio->SetMarkerStyle(20);
      hratio->GetYaxis()->SetRangeUser(0.5,1.5);
      hratio->GetYaxis()->SetTitle("Variation / Nominal");
      hratio->GetYaxis()->SetTitleSize(0.09);
      hratio->GetYaxis()->SetTitleOffset(0.7);
      hratio->GetYaxis()->SetLabelSize(0.08);
      hratio->GetXaxis()->SetTitle("x_{J#gamma}");
      hratio->GetXaxis()->SetTitleSize(0.09);
      hratio->GetXaxis()->SetLabelSize(0.08);
      hratio->Draw("p e");
      TLine * line = new TLine(ana::unfoldXjBins[0],1,ana::unfoldXjBins[ana::nUnfoldXjBins],1);
      line->SetLineStyle(9);
      line->Draw("same");
      c->SaveAs(pdfPath);

      fout->cd();
      hratio->Write();
      hfrac->Write();
    }
  }

  // Purity: combine the five ABCD sideband-boundary sources (purityMembers) into one
  // symmetric quadrature sum - every member here is treated as its own independent
  // symmetrized source (full magnitude feeds both up and down), including
  // narrowISObkg/wideISObkg even though they're a genuine two-sided variation of the same
  // isolation-gap boundary - they are NOT sign-split against each other the way
  // JER/emscale/jes/EMR are (see ana::asymmetricSystagPairs's comment). One combined
  // histogram is enough since there's no asymmetry left to carry: see the header comment
  // and PPG12 Sec. 5.3 for why these five are pre-combined into one systematic rather than
  // entering the grand total independently.
  vector<TH1D*> purityUnc(ana::nPtBins);
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hPurity = (TH1D*)fracDiff[systematics[0]][ipt]->Clone(Form("hpurity_pt%d", ipt));
    hPurity->Reset("ICES");
    for (int b = 1; b <= hPurity->GetNbinsX(); b++) {
      double sumsq = 0;
      for (const string & systag : purityMembers) {
        double v = fracDiff[systag][ipt]->GetBinContent(b);
        sumsq += v*v;
      }
      hPurity->SetBinContent(b, sqrt(sumsq));
    }
    hPurity->GetYaxis()->SetTitle("Purity systematic uncertainty");
    purityUnc[ipt] = hPurity;
    fout->cd();
    hPurity->Write();
  }

  // Total systematic uncertainty per pT bin: quadrature sum of the individual fractional
  // differences at each xJ bin, treating the sources as independent, PLUS the already-
  // combined Purity source above in place of its five individual purityMembers (not
  // double-counted - see the skip below; Purity itself is symmetric, so it feeds both
  // totals like any other symmetric source). Symmetric sources contribute their full
  // (signed-then-squared) magnitude to both totals; each asymmetricSystematics source
  // instead contributes only to whichever total matches its own sign that bin (0 to the
  // other) - up = sqrt(symmetric^2 + sum of positive asymmetric values squared), down =
  // sqrt(symmetric^2 + sum of negative asymmetric values squared). Both totals are always
  // >= 0 by construction, drawn as an asymmetric +/- band.
  vector<TH1D*> totalUncUp(ana::nPtBins), totalUncDown(ana::nPtBins);
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hTotalUp   = (TH1D*)fracDiff[systematics[0]][ipt]->Clone(Form("hquadsum_up_pt%d", ipt));
    TH1D * hTotalDown = (TH1D*)fracDiff[systematics[0]][ipt]->Clone(Form("hquadsum_down_pt%d", ipt));
    hTotalUp->Reset("ICES");
    hTotalDown->Reset("ICES");
    for (int b = 1; b <= hTotalUp->GetNbinsX(); b++) {
      double sumsqUp = 0, sumsqDown = 0;
      for (const string & systag : systematics) {
        if (purityMembers.count(systag)) continue; // folded into purityUnc instead
        double v = fracDiff[systag][ipt]->GetBinContent(b);
        if (asymmetricSystematics.count(systag)) {
          if (v > 0) sumsqUp   += v*v;
          else       sumsqDown += v*v;
        } else {
          sumsqUp   += v*v;
          sumsqDown += v*v;
        }
      }
      double p = purityUnc[ipt]->GetBinContent(b);
      sumsqUp   += p*p;
      sumsqDown += p*p;
      hTotalUp->SetBinContent(b, sqrt(sumsqUp));
      hTotalDown->SetBinContent(b, sqrt(sumsqDown));
    }
    hTotalUp->GetYaxis()->SetTitle("Total systematic uncertainty (up)");
    hTotalDown->GetYaxis()->SetTitle("Total systematic uncertainty (down)");
    totalUncUp[ipt] = hTotalUp;
    totalUncDown[ipt] = hTotalDown;
    fout->cd();
    hTotalUp->Write();
    hTotalDown->Write();
  }

  // Final page: fractional difference vs xJ, one panel per pT bin, every systematic
  // overlaid plus their quadrature sum - the relative size of each systematic source,
  // and the total, is visible at a glance. No error bars on this page - these are
  // per-bin ratios/differences of already-unfolded results, not independent
  // measurements, so per-curve error bars here would overstate what they represent.
  c->Clear();
  c->cd();
  vector<TPad*> pads(nPtBinsUsed);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    double y1 = 1.0 - (ipt+1)*(1.0/nPtBinsUsed);
    double y2 = 1.0 - ipt*(1.0/nPtBinsUsed);
    pads[ipt] = new TPad(Form("psum_%d",ipt),"",0,y1,1,y2);
    pads[ipt]->Draw();
  }
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    // idisplay: 0-based position among the nPtBinsUsed stacked pads (pads[]/"last panel"
    // checks need this), separate from ipt, the real ana::ptBins index (fracDiff[]/
    // totalUncUp[]/totalUncDown[]/ana::ptBins[] all still need the real index).
    int idisplay = ipt - ana::firstUsedPtBin;
    pads[idisplay]->cd();
    pads[idisplay]->SetLeftMargin(.15);
    pads[idisplay]->SetBottomMargin(idisplay == nPtBinsUsed-1 ? 0.2 : 0.02);
    pads[idisplay]->SetTopMargin(0.05);
    gPad->SetTicks(1,1);

    double ymax = 0.05; // headroom floor so a near-flat set of curves isn't over-zoomed
    for (const string & systag : systematics)
      ymax = std::max(ymax, fracDiff[systag][ipt]->GetMaximum());
    for (const string & systag : systematics)
      ymax = std::max(ymax, -fracDiff[systag][ipt]->GetMinimum());
    // purityUnc is a quadrature sum of its members, so it can exceed any single
    // purityMembers curve already covered by the two loops above.
    ymax = std::max(ymax, purityUnc[ipt]->GetMaximum());
    ymax = std::max(ymax, totalUncUp[ipt]->GetMaximum());
    ymax = std::max(ymax, totalUncDown[ipt]->GetMaximum());

    TH1D * frame = (TH1D*)fracDiff[systematics[0]][ipt]->Clone(Form("hsumframe_pt%d",ipt));
    frame->Reset("ICES");
    frame->SetLineColor(kWhite);
    frame->GetYaxis()->SetRangeUser(-1.3, 1.3);
    frame->GetYaxis()->SetTitle("(Var.-Nom.)/Nom.");
    frame->GetYaxis()->SetTitleSize(0.08);
    frame->GetYaxis()->SetTitleOffset(0.8);
    frame->GetYaxis()->SetLabelSize(0.07);
    frame->GetXaxis()->SetTitle(idisplay == nPtBinsUsed-1 ? "x_{J#gamma}" : "");
    frame->GetXaxis()->SetLabelSize(idisplay == nPtBinsUsed-1 ? 0.07 : 0);
    frame->GetXaxis()->SetTitleSize(0.08);
    frame->Draw("p");

    TLine * zero = new TLine(ana::unfoldXjBins[0],0,ana::unfoldXjBins[ana::nUnfoldXjBins],0);
    zero->SetLineStyle(9);
    zero->Draw("same");

    TLegend * ls = new TLegend(.4,.65,.68,.93);
    ls->SetLineWidth(0);
    ls->SetTextSize(0.05);
    for (const DisplayGroup & g : displayGroups) {
      bool firstMember = true;
      for (const string & systag : g.members) {
        // Zero-error display clone: fracDiff[][] itself (real Divide()-propagated errors)
        // still gets Write()'d to the ROOT file above with those errors intact - only the
        // plotted copy is stripped, since a per-bin ratio of two already-unfolded results
        // isn't an independent measurement and error bars here would overstate that.
        TH1D * hf = (TH1D*)fracDiff[systag][ipt]->Clone(Form("hfracdiff_%s_pt%d_disp", systag.c_str(), ipt));
        for (int b = 1; b <= hf->GetNbinsX(); b++) hf->SetBinError(b, 0);
        hf->SetLineColor(g.color);
        hf->SetLineWidth(1);
        hf->Draw("hist same");
        if (firstMember) { ls->AddEntry(hf, g.label.c_str(), "l"); firstMember = false; }

        if (g.symmetric) {
          TH1D * hfNeg = (TH1D*)hf->Clone(Form("hfracdiff_%s_pt%d_disp_neg", systag.c_str(), ipt));
          hfNeg->Scale(-1);
          hfNeg->Draw("hist same");
        }
      }
    }
    // Purity: a single combined symmetric magnitude (not a single systag's real
    // fracDiff), drawn the same curve-plus-negation way as any other symmetric
    // DisplayGroup member above, rather than through that generic loop directly.
    TH1D * hPurityDisp = (TH1D*)purityUnc[ipt]->Clone(Form("hpurity_pt%d_disp", ipt));
    for (int b = 1; b <= hPurityDisp->GetNbinsX(); b++) hPurityDisp->SetBinError(b, 0);
    hPurityDisp->SetLineColor(colorPurity);
    hPurityDisp->SetLineWidth(1);
    hPurityDisp->Draw("hist same");
    ls->AddEntry(hPurityDisp, "Purity", "l");

    TH1D * hPurityDispNeg = (TH1D*)hPurityDisp->Clone(Form("hpurity_pt%d_disp_neg", ipt));
    hPurityDispNeg->Scale(-1);
    hPurityDispNeg->Draw("hist same");

    // Asymmetric total: up and down are independent magnitudes now (not mirror images of
    // each other), so each is its own histogram - "down" negated only for display, to
    // draw on the same signed axis as the fracDiff curves above.
    TH1D * hTotalUpDisp = (TH1D*)totalUncUp[ipt]->Clone(Form("hquadsum_up_pt%d_disp", ipt));
    for (int b = 1; b <= hTotalUpDisp->GetNbinsX(); b++) hTotalUpDisp->SetBinError(b, 0);
    hTotalUpDisp->SetLineColor(kBlack);
    hTotalUpDisp->SetLineWidth(2);
    hTotalUpDisp->SetLineStyle(2);

    TH1D * hTotalDownDisp = (TH1D*)totalUncDown[ipt]->Clone(Form("hquadsum_down_pt%d_disp", ipt));
    for (int b = 1; b <= hTotalDownDisp->GetNbinsX(); b++) hTotalDownDisp->SetBinError(b, 0);
    hTotalDownDisp->Scale(-1);
    hTotalDownDisp->SetLineColor(kBlack);
    hTotalDownDisp->SetLineWidth(2);
    hTotalDownDisp->SetLineStyle(2);

    hTotalDownDisp->Draw("hist same");
    hTotalUpDisp->Draw("hist same");
    ls->AddEntry(hTotalUpDisp, "Total (asym. quad. sum)", "l");
    ls->Draw();

    TLatex * t = new TLatex(.18,.85,Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]));
    t->SetNDC();
    t->SetTextFont(43);
    t->SetTextSize(16);
    t->Draw();
  }
  c->SaveAs(pdfPath);

  c->SaveAs(Form("%s]", pdfPath));
  fout->Close();
  cout << "Done. Wrote " << pdfPath << " and " << rootPath << endl;
}
