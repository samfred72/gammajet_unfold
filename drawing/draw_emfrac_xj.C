#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/unfolder.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/unfold_utility.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/reweight_utility.h"
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
// Reuses unfolder::check_pair/check_keep_MC directly (by constructing a throwaway
// unfolder instance per sample and never calling fill_matrix/unfold/end on it) rather
// than reimplementing the pairing/eta/dphi/xJ-floor cuts by hand, so this can't silently
// drift from the production selection. Reconstructs the same nominal-systag reco-level
// quantities fill_matrix() would (unfolder.cc:149-232): MC gets the always-on 2% EMR
// smear on cluster pT and the truth-derived jet_pt_smear_truth[ir] jet pT; Data gets its
// per-radius JES correction (ana::jesNominal[ir], no jes_high/low offset).
//
// Deliberately simplified relative to the real analysis in two ways, both noted on every
// plot: (1) raw Region A only, no purity (ABCD background) correction - this is a
// detector/reconstruction-level question, not one where the background-subtraction
// machinery is expected to matter; (2) R=0.4 (ir=2) only, matching this project's
// standalone-diagnostic convention (draw_response_matrix.C, draw_refolding.C, etc. all
// hardcode the same radius).
//
// MC Photon5/10/20 combination uses the exact same per-sample cross-section scale
// factors as the real pipeline (drawer::getScale(true, sample) - reads the same
// `scalemap` drawer::combineMC uses), applied on top of the same per-event vz/cluster-pT
// reweighting (src/reweight_utility.h) unfolder.cc itself applies - not re-derived or
// approximated. The Part 2 emfrac reweighting is a further, additional per-event weight
// on top of both of those, not a replacement for either.

const int ir = 2; // R=0.4, nominal - see header comment
const set<int> photonSampleCodes = {5, 10, 20};

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

// One passing (paired, Region A, 13-15 GeV photon) event's emfrac, x_J, and the
// mcWeight*crossSectionScale weight it already carries (1 for Data) - collected once per
// sample and reused for both the coarse ana::emfracBins split (Part 1) and the fine-
// binned reweighting derivation (Part 2), instead of looping the tree twice.
struct EmfracEvent { float emfrac; float xj; float weight; };

// Loops one (trigger,sim) sample once and returns one EmfracEvent per event passing the
// same cuts drawEmfracPage's caller used to rely on fillEmfracXj for (13-15 GeV photon,
// Region A, R=0.4, check_pair) - see the file header comment for exactly what's
// reconstructed and why. extraScale is the cross-section combination factor for MC
// (drawer::getScale(true, sample)), 1 for Data.
vector<EmfracEvent> collectEvents(string trigger, string sim, bool isMCsample, double extraScale) {
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

    // Nominal-systag reco-level photon: MC carries the always-on 2% EMR smear
    // (unfolder.cc's systagEmrSigmaArr nominal value - see that file's header comment);
    // Data is used as-is. Not the same random sequence unfolder.cc itself draws (its own
    // TRandom member is private, and bit-for-bit reproducibility isn't needed for a
    // standalone diagnostic), just the same smearing recipe.
    float recoClusterPt = uf.cluster_pt;
    if (isMCsample) recoClusterPt = rnd.Gaus(recoClusterPt, recoClusterPt * 0.02);

    pho_object maxpho(recoClusterPt, uf.cluster_e, uf.cluster_eta, uf.cluster_phi,
        uf.cluster_showershape[10], uf.cluster_showershape[11], uf.cluster_time,
        uf.cluster_bdt_scores[9],
        pho_object::get_showershape(uf.cluster_showershape, recoClusterPt));

    float recoJetPt = isMCsample ? uf.jet_pt_smear_truth[ir]
                                  : uf.jet_pt_calib[ir] / ana::jesNominal[ir];
    jet_object maxjet(recoJetPt, uf.jet_e[ir], uf.jet_eta[ir], uf.jet_phi[ir],
        uf.jet_emfrac[ir], 0, 0, uf.jet_time[ir]);

    if (!uf.check_pair(maxjet, ir, maxpho, true)) continue;
    // ana::ptBins[0..1) = 13-15 GeV, the migration-only buffer bin below the lowest
    // reported physics bin (ana.h:90-95) - used here on explicit request.
    if (maxpho.pt < ana::ptBins[0] || maxpho.pt >= ana::ptBins[1]) continue;
    if (ana::findabcdBin(maxpho.iso4, maxpho.bdt, 0) != 0) continue; // Region A only - see header comment

    events.push_back({maxjet.emfrac, (float)(maxjet.pt / maxpho.pt), mcWeight * (float)extraScale});
  }
  return events;
}

// Combined Photon5+10+20 pythia sample, each cross-section-scaled via
// drawer::getScale(true, sample) - the same combination drawer::combineMC uses.
vector<EmfracEvent> collectMCPhotonEvents() {
  drawer dScale; // only used for its scalemap accessor (getScale) - opens hists/*.root itself
  vector<EmfracEvent> all;
  for (int code : photonSampleCodes) {
    vector<EmfracEvent> ev = collectEvents(Form("Photon%d", code), "pythia", true, dScale.getScale(true, code));
    all.insert(all.end(), ev.begin(), ev.end());
  }
  return all;
}

// --- Part 1: x_J split by the coarse ana::emfracBins, one page per sample ---

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

// Page 3: Data's and MC's own shape-normalized emfrac distributions (top) and their
// per-bin ratio (bottom) - the ratio IS the reweighting factor w(emfrac) applied in
// Part 2 below, drawn here so it can be inspected on its own before trusting what it does
// to x_J.
TH1D * drawEmfracRatioPage(TCanvas * c, const char * pdfPath,
    const vector<EmfracEvent> & dataEvents, const vector<EmfracEvent> & mcEvents,
    const vector<string> & extraFeatures) {
  TH1D * hDataRaw = fillFineEmfracHist(dataEvents, "hEmfracDataRaw");
  TH1D * hMCRaw = fillFineEmfracHist(mcEvents, "hEmfracMCRaw");

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
  leg->AddEntry(hData, "p+p Run24 Data", "lep");
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
// weight on top of its existing mcWeight*crossSectionScale, with Data alongside as the
// reference the reweighting is meant to approach.
void drawReweightedXjPage(TCanvas * c, const char * pdfPath, TH1D * hRatio,
    const vector<EmfracEvent> & dataEvents, const vector<EmfracEvent> & mcEvents,
    const vector<string> & extraFeatures) {
  TH1D * hData = new TH1D("hXjData", "", ana::nUnfoldXjBins, ana::unfoldXjBins);
  TH1D * hMCBefore = new TH1D("hXjMCBefore", "", ana::nUnfoldXjBins, ana::unfoldXjBins);
  TH1D * hMCAfter = new TH1D("hXjMCAfter", "", ana::nUnfoldXjBins, ana::unfoldXjBins);

  for (const EmfracEvent & ev : dataEvents) hData->Fill(ev.xj, ev.weight);
  for (const EmfracEvent & ev : mcEvents) {
    hMCBefore->Fill(ev.xj, ev.weight);
    double w = hRatio->GetBinContent(hRatio->FindBin(ev.emfrac));
    hMCAfter->Fill(ev.xj, ev.weight * w);
  }

  vector<TH1D*> hRaw = {hData, hMCBefore, hMCAfter};
  vector<string> labels = {"p+p Run24 Data", "Pythia8 MC, before emfrac reweight", "Pythia8 MC, after emfrac reweight"};
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

void draw_emfrac_xj() {
  gStyle->SetOptStat(0);
  TH1::AddDirectory(kFALSE);
  const char * pdfPath = "/home/samson72/sphnx/gammajet_unfold/pdfs/draw_emfrac_xj.pdf";

  vector<string> commonFeatures = {
      Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ana::ptBins[0], ana::ptBins[1]),
      Form("Jet R=%.1f", ana::JetRs[ir]),
      "Region A (no purity correction)",
  };

  TCanvas * c = new TCanvas("c", "", 700, 700);
  c->SaveAs(Form("%s[", pdfPath));

  cout << "Collecting Data:" << endl;
  vector<EmfracEvent> dataEvents = collectEvents("Data", "pythia", false, 1.0);
  cout << "Collecting MC Photon (Photon5+10+20, pythia):" << endl;
  vector<EmfracEvent> mcEvents = collectMCPhotonEvents();

  // --- Part 1: x_J split by the coarse ana::emfracBins ---
  cout << "Page 1 (Data, coarse emfrac split):" << endl;
  vector<TH1D*> hData = binByEmfracCoarse(dataEvents, "data");
  drawEmfracPage(c, pdfPath, hData, {"p+p Run24 Data"}, commonFeatures);

  cout << "Page 2 (MC, coarse emfrac split):" << endl;
  vector<TH1D*> hMC = binByEmfracCoarse(mcEvents, "mc");
  drawEmfracPage(c, pdfPath, hMC, {"Pythia8 #gamma+jet MC"}, commonFeatures);

  // --- Part 2: derive w(emfrac) = normalized Data / normalized MC, reweight MC, compare ---
  cout << "Page 3 (Data vs MC emfrac distributions + ratio = reweighting factor):" << endl;
  TH1D * hRatio = drawEmfracRatioPage(c, pdfPath, dataEvents, mcEvents, commonFeatures);

  cout << "Page 4 (MC x_J before vs after emfrac reweighting):" << endl;
  drawReweightedXjPage(c, pdfPath, hRatio, dataEvents, mcEvents, commonFeatures);

  c->SaveAs(Form("%s]", pdfPath));
  cout << "Wrote " << pdfPath << endl;
}
