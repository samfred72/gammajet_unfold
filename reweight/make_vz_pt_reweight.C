#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/treeuser.h"
#include "../src/pho_object.h"
#include "../src/jet_object.h"
#include "../src/unfold_utility.h"
// The original gammajet project builds its OWN, differently-laid-out drawer/ana classes
// into /home/samson72/root/lib/libgammajet.so, sitting on the same library search path
// as this project's libgammajet_unfold.so. Without forcing which one loads first, ROOT's
// implicit symbol autoload can bind drawer/ana calls to the wrong (mismatched-layout)
// library and segfault - explicit load removes the ambiguity.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Builds the Data/MC reweighting for v_z and photon-cluster p_T applied to every pythia
// MC event (src/reweight_utility.h, unfolder.cc's mcWeight). Both are ratios of shapes
// normalized to unit area.
//
// v_z: all Data events vs all (truth-pT-stitched, cross-section-weighted) Photon5/10/20
// MC events.
//
// Cluster p_T (rederived Sep 28 2026, PPG18 review issue 4): the weight must correct the
// MC SIGNAL photon spectrum, so it is the ratio of
//   Data: the PURITY-CORRECTED region-A cluster-pT spectrum, i.e. the same two-purity
//         subtraction the analysis applies (unfold_utility::purityCorrectCoeffs with the
//         committed ana::getPurity/getPurityC), done in 1 GeV bins with the coefficients
//         of the ana::ptBins purity bin each fine bin falls in
//   MC:   truth-matched (dR < 0.1) region-A cluster pT from Photon5/10/20, stitched and
//         cross-section weighted, with the v_z weight above applied
// both in the analysis selection: region A/C paired with an R = 0.4 jet exactly as
// unfolder::check_pair (jet pT = jet_pt_calib/jesNominal in Data, jet_pt_smear_truth in
// MC, > jet_calib_pt_cut, xJ floor, jet |eta| < 1.1-R, dphi > 7pi/8), photon |eta| <
// photonEtaMax. The previous version divided the shape of ALL Data clusters (any BDT or
// isolation, paired or not) by the leading MC cluster, so its slope followed the fake-
// photon fraction falling with pT rather than any signal mismodeling.
//
// Histogrammed and evaluated in the RAW cluster_pt (no EM-resolution smearing), since
// that is the argument unfolder.cc passes to Reweighter::GetWeight. The ratio is taken in
// all 1 GeV bins over the analysis range [13, 35) GeV and fit with "expo" over all of them
// (stays positive when extrapolated above 35 GeV). In the sparse high-pT bins a
// purity-corrected bin can have zero region-A counts; its Poisson variance would then be
// estimated from region C alone and come out near zero, letting a single empty bin drive
// the fit (30-31 GeV: A = 0, C = 3). Zero-count bins are therefore given a variance of one
// count (the standard floor for an empty Poisson bin) in both A and C. A fit to the same
// ratio in merged bins (1 GeV to 20 GeV, then the 20-25 and 25-35 GeV purity bins) is
// printed as a cross-check.
//
// hPtWeight's first bin is [10, 13) set to the fit value at 13 GeV - Reweighter holds the
// weight flat at hPtWeight's first-bin content below that bin's upper edge (13 GeV) and
// uses fPtWeight above it, so reweight_utility.h needs no change.
//
// Two stages: fillInputs() loops over the trees (~20 min) and writes the raw inputs to
// vz_pt_reweight_inputs.root; the derivation and plots then run from that file. Run
// make_vz_pt_reweight(true) to refill after the trees, purities or selection change.
//
// Circularity note: the purities (and their Photon-MC leakage templates) were themselves
// derived with the previous weight applied to MC. The leakage fractions are ratios within
// a purity bin, so this is a second-order effect; rerun puritymaker after the next
// unfolder pass and, if the purities move, rerun this macro once more.

// Photon acceptance for the weight's selection - the analysis photon acceptance.
const double photonEtaMax = ana::photonEtaCut;
const int reweightIr = 2; // R = 0.4 pairing, the nominal radius
const char * inputsPath = ana::path("reweight/vz_pt_reweight_inputs.root");
const char * outPath    = ana::path("reweight/vz_pt_reweight.root");

// unfolder::check_pair (without the photon-pT-bin lookup, done by the caller)
bool pairedR04(const pho_object & pho, const jet_object & jet, int ptbin) {
  const int ir = reweightIr;
  if (!(jet.pt > ana::jet_calib_pt_cut[ir])) return false;
  float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(ana::jet_calib_pt_cut[ir]/ana::ptBins[ptbin])+1];
  if (jet.pt/pho.pt < lowbin) return false;
  if (fabs(pho.eta) > photonEtaMax) return false;
  if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) return false;
  jet_object j = jet; pho_object p = pho;
  if (j.deltaPhi(p) < ana::oppcut) return false;
  return true;
}

// ---------------- stage 1: fill raw inputs from the trees ----------------
void fillInputs() {
  const double vzcut = ana::vzcut; // 60 cm
  const int nVzBins = 60; // 2 cm bins across the standard analysis vz window
  const int ir = reweightIr;
  const double ptLo = ana::ptBins[0], ptHi = ana::ptBins[ana::nPtBins-1]; // 13, 35
  const int nPtBins = (int)(ptHi - ptLo);

  TH1D * hVzData = new TH1D("hVzData", ";v_{z} [cm];Events (norm.)", nVzBins, -vzcut, vzcut);
  TH1D * hVzMC   = new TH1D("hVzMC",   ";v_{z} [cm];Events (norm.)", nVzBins, -vzcut, vzcut);
  TH1D * hPtA    = new TH1D("hPtA",    ";Cluster p_{T} [GeV];Clusters", nPtBins, ptLo, ptHi);
  TH1D * hPtC    = new TH1D("hPtC",    ";Cluster p_{T} [GeV];Clusters", nPtBins, ptLo, ptHi);
  // MC signal filled in (v_z, pT) so the v_z weight can be applied once it is known
  TH2D * hVzPtMC = new TH2D("hVzPtMC", ";v_{z} [cm];Cluster p_{T} [GeV]", nVzBins, -vzcut, vzcut, nPtBins, ptLo, ptHi);
  for (TH1 * h : std::initializer_list<TH1*>{hVzData, hVzMC, hPtA, hPtC, hVzPtMC}) h->Sumw2();

  {
    treeuser tu("Data");
    Long64_t nentries = tu.t->GetEntriesFast();
    for (Long64_t e = 0; e < nentries; e++) {
      tu.t->GetEntry(e);
      hVzData->Fill(tu.vz);
      if (fabs(tu.vz) > vzcut) continue;
      if (tu.cluster_pt < ptLo || tu.cluster_pt >= ptHi) continue;
      pho_object pho(tu.cluster_pt, tu.cluster_e, tu.cluster_eta, tu.cluster_phi, tu.cluster_showershape[10],
          tu.cluster_showershape[11], tu.cluster_time, tu.cluster_bdt_scores[9],
          pho_object::get_showershape(tu.cluster_showershape, tu.cluster_pt));
      int iabcd = ana::findabcdBin(pho.iso4, pho.bdt, 0);
      if (iabcd != 0 && iabcd != 2) continue;
      jet_object jet(tu.jet_pt_calib[ir]/ana::jesNominal[ir], tu.jet_e[ir], tu.jet_eta[ir], tu.jet_phi[ir], tu.jet_emfrac[ir], 0, 0, 0);
      if (!pairedR04(pho, jet, ana::findPtBin(pho.pt))) continue;
      (iabcd == 0 ? hPtA : hPtC)->Fill(tu.cluster_pt);
    }
    cout << "Data: " << nentries << " entries" << endl;
  }

  // Cross-section weights - same numbers as drawer.h's scalemap[isphoton=1][sample] and
  // insitu/grid_insitu.C's photon_scale, for sim="pythia". Photon cluster pT for the ABCD
  // shower-shape class uses the raw cluster pT, like the weight's own argument.
  const map<string,double> photonScale = {{"Photon5",146359.3},{"Photon10",6944.675},{"Photon20",130.4461}};
  const vector<string> mcSamples = {"Photon5", "Photon10", "Photon20"};
  for (const string & sample : mcSamples) {
    treeuser tu(sample, "pythia");
    double loThresh = tu.threshmap[-1][sample];
    double hiThresh = tu.threshmap_high[-1][sample];
    double scale = photonScale.at(sample);
    Long64_t nentries = tu.t->GetEntriesFast();
    Long64_t nkept = 0;
    for (Long64_t e = 0; e < nentries; e++) {
      tu.t->GetEntry(e);
      if (tu.truth_cluster_pt <= loThresh || tu.truth_cluster_pt >= hiThresh) continue;
      nkept++;
      hVzMC->Fill(tu.vz, scale);
      if (fabs(tu.vz) > vzcut) continue;
      if (tu.cluster_pt < ptLo || tu.cluster_pt >= ptHi) continue;
      pho_object pho(tu.cluster_pt, tu.cluster_e, tu.cluster_eta, tu.cluster_phi, tu.cluster_showershape[10],
          tu.cluster_showershape[11], tu.cluster_time, tu.cluster_bdt_scores[9],
          pho_object::get_showershape(tu.cluster_showershape, tu.cluster_pt));
      if (ana::findabcdBin(pho.iso4, pho.bdt, 0) != 0) continue;
      pho_object tru(tu.truth_cluster_pt, tu.truth_cluster_e, tu.truth_cluster_eta, tu.truth_cluster_phi,
          tu.truth_cluster_iso3, tu.truth_cluster_iso4, 0, 0.99, 2);
      if (pho.deltaR(tru) >= 0.1) continue; // truth-matched signal (unfolder::check_match)
      jet_object jet(tu.jet_pt_smear_truth[ir], tu.jet_e[ir], tu.jet_eta[ir], tu.jet_phi[ir], tu.jet_emfrac[ir], 0, 0, 0);
      if (!pairedR04(pho, jet, ana::findPtBin(pho.pt))) continue;
      hVzPtMC->Fill(tu.vz, tu.cluster_pt, scale);
    }
    cout << sample << ": " << nkept << "/" << nentries
         << " entries kept (truth_cluster_pt in (" << loThresh << "," << hiThresh << ")), scale=" << scale << endl;
  }

  // MC signal pT spectrum with the v_z weight (Data/MC shape ratio) applied
  TH1D * hVzW = (TH1D*)hVzData->Clone("hVzW_tmp");
  hVzW->Scale(1.0/hVzData->Integral());
  TH1D * hVzMCn = (TH1D*)hVzMC->Clone("hVzMCn_tmp");
  hVzMCn->Scale(1.0/hVzMC->Integral());
  hVzW->Divide(hVzMCn);
  TH1D * hPtMC = new TH1D("hPtMC", ";Cluster p_{T} [GeV];Clusters", nPtBins, ptLo, ptHi);
  hPtMC->Sumw2();
  for (int bv = 1; bv <= nVzBins; bv++) {
    double wv = hVzW->GetBinContent(bv);
    if (!std::isfinite(wv)) wv = 1;
    for (int bp = 1; bp <= nPtBins; bp++) {
      double c = hVzPtMC->GetBinContent(bv, bp), e = hVzPtMC->GetBinError(bv, bp);
      hPtMC->SetBinContent(bp, hPtMC->GetBinContent(bp) + wv*c);
      hPtMC->SetBinError(bp, hypot(hPtMC->GetBinError(bp), wv*e));
    }
  }

  TFile * f = new TFile(inputsPath, "recreate");
  for (TH1 * h : std::initializer_list<TH1*>{hVzData, hVzMC, hPtA, hPtC, hPtMC}) h->Write();
  f->Close();
  cout << "Saved " << inputsPath << endl;
}

namespace {
  void sanitize(TH1 * h, double fill = 1) {
    for (int b = 0; b <= h->GetNbinsX() + 1; b++)
      if (!std::isfinite(h->GetBinContent(b)) || !std::isfinite(h->GetBinError(b))) { h->SetBinContent(b, fill); h->SetBinError(b, 0); }
  }
  // top (distributions) / bottom (ratio) pads with explicit margins
  void twoPads(TCanvas * c, TPad *& top, TPad *& bot) {
    c->cd();
    top = new TPad(Form("%s_top", c->GetName()), "", 0, 0.32, 1, 1);
    bot = new TPad(Form("%s_bot", c->GetName()), "", 0, 0, 1, 0.32);
    top->SetLeftMargin(0.15); top->SetRightMargin(0.04); top->SetTopMargin(0.05); top->SetBottomMargin(0.02);
    bot->SetLeftMargin(0.15); bot->SetRightMargin(0.04); bot->SetTopMargin(0.03); bot->SetBottomMargin(0.33);
    top->SetTicks(1, 1); bot->SetTicks(1, 1);
    top->Draw(); bot->Draw();
  }
  void styleBottomAxes(TH1 * h, const char * ytitle, double ylo, double yhi) {
    h->GetYaxis()->SetTitle(ytitle); h->GetYaxis()->SetRangeUser(ylo, yhi); h->GetYaxis()->SetNdivisions(505);
    h->GetYaxis()->SetTitleSize(0.1); h->GetYaxis()->SetLabelSize(0.09); h->GetYaxis()->SetTitleOffset(0.65);
    h->GetXaxis()->SetTitleSize(0.12); h->GetXaxis()->SetLabelSize(0.1); h->GetXaxis()->SetTitleOffset(1.1);
  }
  void styleTopAxes(TH1 * h) {
    h->GetXaxis()->SetLabelSize(0); h->GetXaxis()->SetTitleSize(0);
    h->GetYaxis()->SetTitleSize(0.055); h->GetYaxis()->SetLabelSize(0.045); h->GetYaxis()->SetTitleOffset(1.25);
  }
}

// ---------------- stage 2: derive the weights and draw ----------------
void make_vz_pt_reweight(bool refill = false) {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  const int ir = reweightIr;
  const double ptLo = ana::ptBins[0], ptHi = ana::ptBins[ana::nPtBins-1]; // 13, 35

  if (refill || gSystem->AccessPathName(inputsPath)) fillInputs();
  TFile * fin = TFile::Open(inputsPath, "read");
  TH1D * hVzData = (TH1D*)fin->Get("hVzData")->Clone("hVzData");
  TH1D * hVzMC   = (TH1D*)fin->Get("hVzMC")->Clone("hVzMC");
  TH1D * hPtA    = (TH1D*)fin->Get("hPtA")->Clone("hPtA");
  TH1D * hPtC    = (TH1D*)fin->Get("hPtC")->Clone("hPtC");
  TH1D * hPtMC   = (TH1D*)fin->Get("hPtMC")->Clone("hPtMC");
  for (TH1 * h : std::initializer_list<TH1*>{hVzData, hVzMC, hPtA, hPtC, hPtMC}) h->SetDirectory(0);
  fin->Close();
  const int nPtBins = hPtA->GetNbinsX();

  // previous (pre-Sep-28, all-cluster) weight, kept for the comparison curve
  TF1 * fOld = nullptr;
  {
    TFile * fprev = TFile::Open(outPath, "read");
    if (fprev && !fprev->IsZombie()) {
      TObject * o = fprev->Get("fPtWeight_previous");
      if (o) fOld = (TF1*)o->Clone("fPtWeight_previous");
      fprev->Close();
    }
  }

  // ---- v_z weight ----
  hVzData->Scale(1.0 / hVzData->Integral());
  hVzMC->Scale(1.0 / hVzMC->Integral());
  TH1D * hVzWeight = (TH1D*)hVzData->Clone("hVzWeight");
  hVzWeight->SetTitle(";v_{z} [cm];Data / MC weight");
  hVzWeight->Divide(hVzMC);
  sanitize(hVzWeight); // no data to reweight against - leave MC unweighted there

  // ---- purity-corrected Data spectrum, all 1 GeV bins ----
  TH1D * hPtData = (TH1D*)hPtA->Clone("hPtData");
  hPtData->Reset("ICES");
  hPtData->SetTitle(";Cluster p_{T} [GeV];Clusters (norm.)");
  printf("\nPurity-corrected Data cluster-pT spectrum (R = %.1f pairing):\n", ana::JetRs[ir]);
  for (int ib = 0; ib < ana::nPtBins - 1; ib++) { // 13-15, 15-20, 20-25, 25-35
    double lo = ana::ptBins[ib], hi = ana::ptBins[ib+1];
    int b0 = hPtA->FindBin(lo + 1e-3), b1 = hPtA->FindBin(hi - 1e-3);
    double NA = hPtA->Integral(b0, b1), NC = hPtC->Integral(b0, b1);
    float pA = ana::getPurity(lo, hi, "nominal", ir), pC = ana::getPurityC(lo, hi, "nominal", ir);
    float cA, cC;
    bool exact = unfold_utility::purityCorrectCoeffs(pA, pC, NA, NC, cA, cC);
    for (int b = b0; b <= b1; b++) {
      double a = hPtA->GetBinContent(b), c = hPtC->GetBinContent(b);
      hPtData->SetBinContent(b, cA*a - cC*c);
      // counting statistics only; empty bins get a one-count variance (see header)
      hPtData->SetBinError(b, sqrt(cA*cA*std::max(a, 1.0) + cC*cC*std::max(c, 1.0)));
    }
    printf("  %4.0f-%-4.0f GeV: N_A %6.0f N_C %6.0f  P_A %.3f P_C %.3f%s -> signal %.0f\n", lo, hi, NA, NC, pA, pC,
           exact ? "" : " (single-purity fallback)", cA*NA - cC*NC);
  }
  double nD = hPtData->Integral(), nM = hPtMC->Integral();
  hPtData->Scale(1.0/nD);
  hPtMC->Scale(1.0/nM);
  hPtMC->SetTitle(";Cluster p_{T} [GeV];Clusters (norm.)");

  TH1D * hPtRatio = (TH1D*)hPtData->Clone("hPtRatio");
  hPtRatio->SetTitle(";Cluster p_{T} [GeV];Data / MC weight");
  hPtRatio->Divide(hPtMC);
  sanitize(hPtRatio);

  TF1 * fPtWeight = new TF1("fPtWeight", "expo", ptLo, ptHi);
  hPtRatio->Fit(fPtWeight, "RQ0");
  printf("\nexpo fit, all %d 1-GeV bins over %.0f-%.0f GeV: slope %.4f +- %.4f /GeV, chi2/ndf %.1f/%d, w(15)=%.3f w(20)=%.3f w(35)=%.3f\n",
         nPtBins, ptLo, ptHi, fPtWeight->GetParameter(1), fPtWeight->GetParError(1), fPtWeight->GetChisquare(), fPtWeight->GetNDF(),
         fPtWeight->Eval(15), fPtWeight->Eval(20), fPtWeight->Eval(35));
  { // cross-check: merged bins (1 GeV to 20 GeV, then the 20-25 and 25-35 GeV purity bins)
    const vector<double> e = {13, 14, 15, 16, 17, 18, 19, 20, 25, 35};
    TH1D * d = (TH1D*)hPtData->Rebin(e.size()-1, "hPtDataMerged", e.data());
    TH1D * m = (TH1D*)hPtMC->Rebin(e.size()-1, "hPtMCMerged", e.data());
    d->Divide(m);
    TF1 * fm = new TF1("fPtWeightMerged", "expo", ptLo, ptHi);
    d->Fit(fm, "RQ0");
    printf("cross-check, merged bins: slope %.4f +- %.4f /GeV, chi2/ndf %.1f/%d, w(15)=%.3f w(35)=%.3f\n",
           fm->GetParameter(1), fm->GetParError(1), fm->GetChisquare(), fm->GetNDF(), fm->Eval(15), fm->Eval(35));
  }
  if (fOld) printf("previous weight (all clusters): w(15)=%.3f w(20)=%.3f w(35)=%.3f\n", fOld->Eval(15), fOld->Eval(20), fOld->Eval(35));

  // hPtWeight in the layout Reweighter expects: bin 1 = [10, 13) at the fit value at 13 GeV
  // (held flat below 13 GeV), then the measured 1 GeV ratio bins (for display; above 13 GeV
  // Reweighter evaluates fPtWeight).
  vector<double> edges = {ana::cluster_pt_cut};
  for (int b = 0; b <= nPtBins; b++) edges.push_back(ptLo + b);
  TH1D * hPtWeight = new TH1D("hPtWeight", ";Cluster p_{T} [GeV];Data / MC weight", edges.size()-1, edges.data());
  hPtWeight->SetBinContent(1, fPtWeight->Eval(ptLo));
  for (int b = 1; b <= nPtBins; b++) { hPtWeight->SetBinContent(b+1, hPtRatio->GetBinContent(b)); hPtWeight->SetBinError(b+1, hPtRatio->GetBinError(b)); }

  TFile * fout = new TFile(outPath, "recreate");
  for (TH1 * h : std::initializer_list<TH1*>{hVzData, hVzMC, hVzWeight, hPtA, hPtC, hPtData, hPtMC, hPtRatio, hPtWeight}) h->Write();
  fPtWeight->Write("fPtWeight");
  if (fOld) fOld->Write("fPtWeight_previous");
  fout->Close();
  cout << "Saved " << outPath << endl;

  drawer d("pythia", "nominal"); // ctor args only pick files for other helpers; labels come from drawAll

  // ---- plot 1: v_z ----
  {
    TCanvas * c = new TCanvas("cvz", "", 700, 700);
    TPad * top, * bot; twoPads(c, top, bot);
    top->cd();
    hVzData->SetLineColor(kBlack); hVzData->SetMarkerColor(kBlack); hVzData->SetMarkerStyle(20); hVzData->SetMarkerSize(0.8);
    hVzMC->SetLineColor(kRed+1); hVzMC->SetLineWidth(2);
    styleTopAxes(hVzData);
    hVzData->GetYaxis()->SetTitle("Events (norm.)");
    hVzData->GetYaxis()->SetRangeUser(0, std::max(hVzData->GetMaximum(), hVzMC->GetMaximum()) * 1.6);
    hVzData->Draw("e");
    hVzMC->Draw("hist same");
    hVzData->Draw("e same");
    TLegend * l = new TLegend(0.6, 0.72, 0.95, 0.88); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.045);
    l->AddEntry(hVzData, "p+p Run24 Data", "lp");
    l->AddEntry(hVzMC, "Pythia8 #gamma+jet MC", "l");
    l->Draw();
    d.drawAll({}, {"all events, |v_{z}| < 60 cm"}, .19, .88, 16, gPad->GetWh() * 0.8);
    bot->cd();
    TH1D * r = (TH1D*)hVzWeight->Clone("hVzWeightDraw");
    r->SetLineColor(kBlack); r->SetMarkerColor(kBlack); r->SetMarkerStyle(20); r->SetMarkerSize(0.8);
    styleBottomAxes(r, "Data / MC", 0, 1.15*r->GetMaximum());
    r->Draw("e");
    TLine * one = new TLine(r->GetXaxis()->GetXmin(), 1, r->GetXaxis()->GetXmax(), 1); one->SetLineStyle(2); one->Draw();
    c->SaveAs(ana::path("pdfs/reweight_vz.pdf"));
  }

  // ---- plot 2: cluster p_T ----
  {
    const double xLo = ana::cluster_pt_cut, xHi = ptHi; // show the flat extension below 13 GeV too
    TCanvas * c = new TCanvas("cpt", "", 700, 700);
    TPad * top, * bot; twoPads(c, top, bot);
    top->cd(); top->SetLogy();
    TH1D * fr = new TH1D("frPtTop", ";Cluster p_{T} [GeV];Clusters (norm.)", 1, xLo, xHi);
    styleTopAxes(fr);
    double ymin = 1e30;
    for (int b = 1; b <= nPtBins; b++) { if (hPtData->GetBinContent(b) > 0) ymin = std::min(ymin, hPtData->GetBinContent(b)); if (hPtMC->GetBinContent(b) > 0) ymin = std::min(ymin, hPtMC->GetBinContent(b)); }
    fr->SetMinimum(0.3*ymin); fr->SetMaximum(20*std::max(hPtData->GetMaximum(), hPtMC->GetMaximum()));
    fr->Draw("axis");
    hPtMC->SetLineColor(kRed+1); hPtMC->SetLineWidth(2);
    hPtData->SetLineColor(kBlack); hPtData->SetMarkerColor(kBlack); hPtData->SetMarkerStyle(20); hPtData->SetMarkerSize(0.8);
    hPtMC->Draw("hist same");
    hPtData->Draw("e same");
    TLegend * l = new TLegend(0.42, 0.64, 0.95, 0.8); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.038);
    l->AddEntry(hPtData, "Data, purity-corrected region A", "lp");
    l->AddEntry(hPtMC, "Pythia8 #gamma+jet, truth-matched region A", "l");
    l->Draw();
    d.drawAll({}, {Form("paired, Jet R=%.1f, |#eta^{#gamma}| < %.1f", ana::JetRs[ir], photonEtaMax),
                   "MC v_{z}-weighted"}, .19, .88, 16, gPad->GetWh() * 0.8);
    bot->cd();
    TH1D * frb = new TH1D("frPtBot", ";Cluster p_{T} [GeV];Data / MC", 1, xLo, xHi);
    double rmin = 0, rmax = 2.2;
    for (int b = 1; b <= nPtBins; b++) { rmin = std::min(rmin, hPtRatio->GetBinContent(b) - hPtRatio->GetBinError(b)); }
    styleBottomAxes(frb, "Data / MC", std::max(rmin*1.1, -1.0), rmax);
    frb->Draw("axis");
    hPtRatio->SetLineColor(kBlack); hPtRatio->SetMarkerColor(kBlack); hPtRatio->SetMarkerStyle(20); hPtRatio->SetMarkerSize(0.8);
    hPtRatio->Draw("e same");
    TF1 * fDraw = new TF1("fDraw", "expo", ptLo, xHi); fDraw->SetParameters(fPtWeight->GetParameters()); fDraw->SetLineColor(kBlue); fDraw->SetLineWidth(2);
    TF1 * fFlat = new TF1("fFlat", Form("%f", fPtWeight->Eval(ptLo)), xLo, ptLo); fFlat->SetLineColor(kBlue); fFlat->SetLineStyle(2); fFlat->SetLineWidth(2);
    fDraw->Draw("same"); fFlat->Draw("same");
    if (fOld) { fOld->SetLineColor(kGray+1); fOld->SetLineStyle(7); fOld->SetLineWidth(2); fOld->SetRange(11, xHi); fOld->Draw("same"); }
    TLine * one = new TLine(xLo, 1, xHi, 1); one->SetLineStyle(3); one->Draw();
    TLegend * lb = new TLegend(0.17, 0.36, 0.6, 0.56); lb->SetBorderSize(0); lb->SetFillStyle(0); lb->SetTextSize(0.07);
    lb->AddEntry(fDraw, "weight: expo fit (flat below 13 GeV)", "l");
    if (fOld) lb->AddEntry(fOld, "previous weight", "l");
    lb->Draw();
    c->SaveAs(ana::path("pdfs/reweight_pt.pdf"));
  }
}
