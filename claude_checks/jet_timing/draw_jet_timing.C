#include "../../src/ana.h"
#include "../../src/drawer.h"
#include "../../src/pho_object.h"
#include "../../src/jet_object.h"
// Explicit load - see drawing/draw_final_result.C: the sibling gammajet project's
// libgammajet.so has same-named classes, and implicit autoload can bind to it. Run
// interpreted (root -b -q draw_jet_timing.C), never with ACLiC "+".
R__LOAD_LIBRARY(libgammajet_unfold.so);

// PPG18 review round 2, issue 1: jet timing window 0 < t_MBD - t_jet < 4 ns.
//
// The window is applied at treemaking (gammajet_treemaking CaloAna.cc:433, data only): a
// jet failing it is skipped and the next-hardest in-window jet is stored instead. So the
// trees here only ever contain jets INSIDE the window - what was removed cannot be seen
// directly. This macro characterizes the surviving distribution to answer: are the
// entries piled against the 0 ns edge in-time recoil jets (i.e. the window is cutting
// into signal), or out-of-time background (the window doing its job)?
//
//   p1  Delta t for photon cluster and recoil jet, region A paired, window drawn
//   p2  Delta t_jet in slices of corrected jet pT
//   p3  Delta t_jet in slices of jet EM fraction (EMCal vs HCal timing offset)
//   p4  mean Delta t_jet vs jet pT and vs EM fraction (profiles)
//   p5  Delta t_jet vs Delta t_cluster (common MBD-t0 jitter / correlation)
//   p6  photon-jet Delta phi for edge (0-0.5 ns) vs core (1-3 ns) jets, no Delta phi cut
//   p7  xJ for edge vs core jets, paired
//   p8  estimated fraction of in-time jets lost, vs jet pT and vs EM fraction
//       (truncated Gaussian fit inside the window, integrated outside it)
//   p9  size of the implied reco-level xJ shape change per photon-pT bin: region-A paired
//       xJ reweighted by 1/(1 - loss(jet pT)) over nominal. Treats a lost jet as a lost
//       event (in reality the next in-window jet may be stored instead), and uses the
//       Gaussian model of p8 - an order-of-magnitude estimate, not a correction.
//
// Jet pT is the analysis-level corrected pT, jet_pt_calib[ir]/ana::jesNominal[ir], as in
// unfolder.cc. "Paired" follows unfolder::check_pair (photon |eta| < 1.1, jet |eta| <
// 1.1-R, Delta phi > 7pi/8, same xJ floor). Photon regions from ana::findabcdBin(iso4,
// bdt, 0) exactly as unfolder.cc builds pho_object.

namespace {
  const double tLow = ana::tlowcut, tHigh = ana::thighcut; // 0, 4 ns (the treemaking window)

  // NaN/Inf poisons ROOT's auto-ranging (CLAUDE.md) - zero such bins before drawing.
  void sanitize(TH1 * h) {
    for (int b = 0; b <= h->GetNcells(); b++) {
      if (!std::isfinite(h->GetBinContent(b)) || !std::isfinite(h->GetBinError(b))) {
        h->SetBinContent(b, 0); h->SetBinError(b, 0);
      }
    }
  }
  void unitNorm(TH1 * h) { if (h->Integral() > 0) h->Scale(1.0/h->Integral(), "width"); }
  void windowLines(double ymin, double ymax) {
    for (double x : {tLow, tHigh}) {
      TLine * l = new TLine(x, ymin, x, ymax);
      l->SetLineStyle(2); l->SetLineColor(kGray+2); l->Draw();
    }
  }

  // Fraction of a Gaussian lying outside [tLow, tHigh], fitted only to the in-window part
  // of h. Only meaningful under the assumption that the in-time population is Gaussian -
  // an upper-bound-style model estimate, not a measurement (see README.md).
  struct LossEst { double frac, err, mean, sigma; };
  LossEst gausLoss(TH1 * h, const char * name) {
    LossEst r{0, 0, 0, 0};
    if (h->GetEntries() < 50) return r;
    TF1 * g = new TF1(name, "gaus", tLow, tHigh);
    g->SetParameters(h->GetMaximum(), h->GetMean(), h->GetRMS());
    h->Fit(g, "QRN0L");
    r.mean = g->GetParameter(1); r.sigma = fabs(g->GetParameter(2));
    auto inside = [&](double m, double s) {
      return 0.5*(TMath::Erf((tHigh-m)/(sqrt(2)*s)) - TMath::Erf((tLow-m)/(sqrt(2)*s)));
    };
    r.frac = 1 - inside(r.mean, r.sigma);
    // error: propagate mean/sigma fit errors numerically (uncorrelated approx.)
    double dm = g->GetParError(1), ds = g->GetParError(2);
    double fm = 1 - inside(r.mean+dm, r.sigma), fs = 1 - inside(r.mean, r.sigma+ds);
    r.err = sqrt(pow(fm-r.frac,2) + pow(fs-r.frac,2));
    return r;
  }
}

void draw_jet_timing(int ir = 2)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();

  TFile * fin = TFile::Open(ana::path("trees/gammajet_Data.root"), "read");
  TTree * t = (TTree*)fin->Get("towerntup");

  float vz, mbd_time, cluster_pt, cluster_e, cluster_eta, cluster_phi, cluster_time;
  float cluster_showershape[12], cluster_bdt_scores[11];
  float jet_pt_calib[7], jet_e[7], jet_eta[7], jet_phi[7], jet_emfrac[7], jet_time[7];
  t->SetBranchStatus("*", 0);
  for (const char * b : {"vz","mbd_time","cluster_pt","cluster_e","cluster_eta","cluster_phi",
        "cluster_time","cluster_showershape","cluster_bdt_scores","jet_pt_calib","jet_e",
        "jet_eta","jet_phi","jet_emfrac","jet_time"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("vz", &vz);
  t->SetBranchAddress("mbd_time", &mbd_time);
  t->SetBranchAddress("cluster_pt", &cluster_pt);
  t->SetBranchAddress("cluster_e", &cluster_e);
  t->SetBranchAddress("cluster_eta", &cluster_eta);
  t->SetBranchAddress("cluster_phi", &cluster_phi);
  t->SetBranchAddress("cluster_time", &cluster_time);
  t->SetBranchAddress("cluster_showershape", cluster_showershape);
  t->SetBranchAddress("cluster_bdt_scores", cluster_bdt_scores);
  t->SetBranchAddress("jet_pt_calib", jet_pt_calib);
  t->SetBranchAddress("jet_e", jet_e);
  t->SetBranchAddress("jet_eta", jet_eta);
  t->SetBranchAddress("jet_phi", jet_phi);
  t->SetBranchAddress("jet_emfrac", jet_emfrac);
  t->SetBranchAddress("jet_time", jet_time);

  const double phoLow = 15, phoHigh = 35; // reported photon-pT range
  const int nJpt = 6;
  const double jptEdges[nJpt+1] = {5, 8, 11, 14, 18, 24, 35};
  const int nEmf = 5;
  const double emfEdges[nEmf+1] = {0, 0.3, 0.5, 0.7, 0.85, 1.0};
  const int nDt = 60; const double dtLo = -1, dtHi = 5;

  // region A paired (headline) and all-ABCD paired (for slices - jet timing is a
  // property of the jet, not of the photon ID region, so this just adds statistics)
  TH1D * hDtJetA = new TH1D("hDtJetA", ";t_{MBD} - t_{object} [ns];normalized", nDt, dtLo, dtHi);
  TH1D * hDtClA  = new TH1D("hDtClA",  ";t_{MBD} - t_{object} [ns];normalized", nDt, dtLo, dtHi);
  TH1D * hDtJetPt[nJpt], * hDtJetEmf[nEmf];
  for (int i = 0; i < nJpt; i++) hDtJetPt[i] = new TH1D(Form("hDtJetPt%d",i), ";t_{MBD} - t_{jet} [ns];normalized", nDt, dtLo, dtHi);
  for (int i = 0; i < nEmf; i++) hDtJetEmf[i] = new TH1D(Form("hDtJetEmf%d",i), ";t_{MBD} - t_{jet} [ns];normalized", nDt, dtLo, dtHi);
  TProfile * pDtVsPt  = new TProfile("pDtVsPt",  ";corrected jet p_{T} [GeV];#LT t_{MBD} - t_{jet} #GT [ns]", nJpt, jptEdges);
  TProfile * pDtVsEmf = new TProfile("pDtVsEmf", ";jet EM fraction;#LT t_{MBD} - t_{jet} #GT [ns]", nEmf, emfEdges);
  TH2D * hJetVsCl = new TH2D("hJetVsCl", ";t_{MBD} - t_{cluster} [ns];t_{MBD} - t_{jet} [ns]", 50, -1, 5, 50, -1, 5);
  TH1D * hDphiEdge = new TH1D("hDphiEdge", ";#Delta#phi(#gamma, jet);normalized", 32, 0, M_PI);
  TH1D * hDphiCore = new TH1D("hDphiCore", ";#Delta#phi(#gamma, jet);normalized", 32, 0, M_PI);
  TH1D * hDphiHigh = new TH1D("hDphiHigh", ";#Delta#phi(#gamma, jet);normalized", 32, 0, M_PI);
  TH1D * hXjEdge = new TH1D("hXjEdge", ";x_{J#gamma};normalized", ana::nUnfoldXjBins, ana::unfoldXjBins);
  TH1D * hXjCore = new TH1D("hXjCore", ";x_{J#gamma};normalized", ana::nUnfoldXjBins, ana::unfoldXjBins);
  long nA = 0, nAedge = 0, nClEdge = 0;
  struct PairA { int ipt; float xj, jetPt; };
  std::vector<PairA> pairsA; // region-A paired, for p9

  Long64_t n = t->GetEntries();
  for (Long64_t e = 0; e < n; e++) {
    t->GetEntry(e);
    if (fabs(vz) > ana::vzcut) continue;
    if (cluster_pt < phoLow || cluster_pt >= phoHigh) continue;
    if (jet_pt_calib[ir] <= 0) continue; // no in-window jet stored for this radius

    pho_object pho(cluster_pt, cluster_e, cluster_eta, cluster_phi,
        cluster_showershape[10], cluster_showershape[11], cluster_time,
        cluster_bdt_scores[9], pho_object::get_showershape(cluster_showershape, cluster_pt));
    int iabcd = ana::findabcdBin(pho.iso4, pho.bdt, 0);
    if (iabcd < 0) continue;
    if (fabs(pho.eta) > ana::etacut) continue;

    float jetPt = jet_pt_calib[ir] / ana::jesNominal[ir];
    jet_object jet(jetPt, jet_e[ir], jet_eta[ir], jet_phi[ir], jet_emfrac[ir], 0, 0, jet_time[ir]);
    if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) continue;

    float dtJet = mbd_time - jet_time[ir];
    float dtCl  = mbd_time - cluster_time;
    float dphi  = jet.deltaPhi(pho);
    float xj    = jetPt / pho.pt;

    // same xJ floor as unfolder::check_pair
    int ptbin = ana::findPtBin(pho.pt);
    float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(ana::jet_calib_pt_cut[ir]/ana::ptBins[ptbin])+1];
    bool aboveFloor = xj >= lowbin;

    bool edge = dtJet < 0.5, core = dtJet > 1 && dtJet < 3, high = dtJet > 3.5;
    if (iabcd == 0 && aboveFloor) {
      if (edge) hDphiEdge->Fill(dphi);
      if (core) hDphiCore->Fill(dphi);
      if (high) hDphiHigh->Fill(dphi);
    }
    if (dphi < ana::oppcut || !aboveFloor) continue;

    // --- paired from here on ---
    int ijp = -1; for (int i = 0; i < nJpt; i++) if (jetPt >= jptEdges[i] && jetPt < jptEdges[i+1]) ijp = i;
    int iem = -1; for (int i = 0; i < nEmf; i++) if (jet_emfrac[ir] >= emfEdges[i] && jet_emfrac[ir] <= emfEdges[i+1]) { iem = i; break; }
    if (ijp >= 0) hDtJetPt[ijp]->Fill(dtJet);
    if (iem >= 0) hDtJetEmf[iem]->Fill(dtJet);
    if (jetPt >= jptEdges[0] && jetPt < jptEdges[nJpt]) pDtVsPt->Fill(jetPt, dtJet);
    pDtVsEmf->Fill(jet_emfrac[ir], dtJet);

    if (iabcd != 0) continue;
    nA++;
    if (dtJet < 0.25) nAedge++;
    if (dtCl < 0.25) nClEdge++;
    hDtJetA->Fill(dtJet);
    hDtClA->Fill(dtCl);
    hJetVsCl->Fill(dtCl, dtJet);
    { int ipt = -1; for (int i = 0; i < ana::nPtBinsUsed; i++) if (pho.pt >= ana::ptBinsUsed[i] && pho.pt < ana::ptBinsUsed[i+1]) ipt = i;
      if (ipt >= 0) pairsA.push_back({ipt, xj, jetPt}); }
    if (edge) hXjEdge->Fill(xj);
    if (core) hXjCore->Fill(xj);
  }

  // ---------------- loss estimates ----------------
  TGraphErrors * gLossPt = new TGraphErrors();
  TGraphErrors * gLossEmf = new TGraphErrors();
  cout << "\nGaussian-extrapolated fraction outside [0,4] ns (in-time hypothesis):" << endl;
  for (int i = 0; i < nJpt; i++) {
    LossEst r = gausLoss(hDtJetPt[i], Form("gPt%d",i));
    double x = 0.5*(jptEdges[i]+jptEdges[i+1]);
    gLossPt->SetPoint(i, x, r.frac); gLossPt->SetPointError(i, 0.5*(jptEdges[i+1]-jptEdges[i]), r.err);
    printf("  jet pT %4.0f-%4.0f GeV: N=%6.0f  mean=%.2f sigma=%.2f  lost=%.3f +- %.3f\n",
           jptEdges[i], jptEdges[i+1], hDtJetPt[i]->GetEntries(), r.mean, r.sigma, r.frac, r.err);
  }
  for (int i = 0; i < nEmf; i++) {
    LossEst r = gausLoss(hDtJetEmf[i], Form("gEmf%d",i));
    double x = 0.5*(emfEdges[i]+emfEdges[i+1]);
    gLossEmf->SetPoint(i, x, r.frac); gLossEmf->SetPointError(i, 0.5*(emfEdges[i+1]-emfEdges[i]), r.err);
    printf("  EM frac %.2f-%.2f: N=%6.0f  mean=%.2f sigma=%.2f  lost=%.3f +- %.3f\n",
           emfEdges[i], emfEdges[i+1], hDtJetEmf[i]->GetEntries(), r.mean, r.sigma, r.frac, r.err);
  }
  printf("\nregion A paired, %.0f-%.0f GeV: %ld jets; first 0.25 ns of window: %ld jets vs %ld clusters\n",
         phoLow, phoHigh, nA, nAedge, nClEdge);
  double fEdge = hDphiEdge->Integral(hDphiEdge->FindBin(ana::oppcut+1e-3), hDphiEdge->GetNbinsX()) / std::max(1.0, hDphiEdge->Integral());
  double fCore = hDphiCore->Integral(hDphiCore->FindBin(ana::oppcut+1e-3), hDphiCore->GetNbinsX()) / std::max(1.0, hDphiCore->Integral());
  double fHigh = hDphiHigh->Integral(hDphiHigh->FindBin(ana::oppcut+1e-3), hDphiHigh->GetNbinsX()) / std::max(1.0, hDphiHigh->Integral());
  printf("fraction with dphi > 7pi/8 (region A, no dphi cut): edge(<0.5ns) %.3f (N=%.0f), core(1-3ns) %.3f (N=%.0f), high(>3.5ns) %.3f (N=%.0f)\n",
         fEdge, hDphiEdge->Integral(), fCore, hDphiCore->Integral(), fHigh, hDphiHigh->Integral());
  printf("mean xJ: edge %.3f +- %.3f (N=%.0f), core %.3f +- %.3f (N=%.0f)\n",
         hXjEdge->GetMean(), hXjEdge->GetMeanError(), hXjEdge->GetEntries(),
         hXjCore->GetMean(), hXjCore->GetMeanError(), hXjCore->GetEntries());
  printf("KS(xJ edge vs core) = %.3f;  KS(dphi edge vs core) = %.3f\n",
         hXjEdge->KolmogorovTest(hXjCore), hDphiEdge->KolmogorovTest(hDphiCore));

  // p9 inputs: reweight region-A paired xJ by 1/(1 - loss(jet pT)), loss from the p8 fit in
  // the jet's pT slice (lowest/highest slice used below/above the sliced range).
  TH1D * hXjNom[ana::nPtBinsUsed], * hXjRw[ana::nPtBinsUsed];
  for (int i = 0; i < ana::nPtBinsUsed; i++) {
    hXjNom[i] = new TH1D(Form("hXjNom%d",i), ";x_{J#gamma};reweighted / nominal (shape)", ana::nUnfoldXjBins, ana::unfoldXjBins);
    hXjRw[i]  = new TH1D(Form("hXjRw%d",i),  ";x_{J#gamma};reweighted / nominal (shape)", ana::nUnfoldXjBins, ana::unfoldXjBins);
  }
  for (auto & p : pairsA) {
    int ijp = 0; for (int i = 0; i < nJpt; i++) if (p.jetPt >= jptEdges[i]) ijp = i;
    double loss = gLossPt->GetY()[ijp];
    hXjNom[p.ipt]->Fill(p.xj);
    hXjRw[p.ipt]->Fill(p.xj, loss < 1 ? 1.0/(1.0-loss) : 1.0);
  }
  cout << "\np9: reco-level region-A xJ shape, reweighted/nominal (Gaussian loss model):" << endl;
  for (int i = 0; i < ana::nPtBinsUsed; i++) {
    unitNorm(hXjNom[i]); unitNorm(hXjRw[i]);
    hXjRw[i]->Divide(hXjNom[i]);
    printf("  %.0f-%.0f GeV:", ana::ptBinsUsed[i], ana::ptBinsUsed[i+1]);
    for (int b = 1; b <= hXjRw[i]->GetNbinsX(); b++)
      if (hXjNom[i]->GetBinContent(b) > 0) printf(" [%.1f-%.1f] %.3f", hXjRw[i]->GetBinLowEdge(b), hXjRw[i]->GetXaxis()->GetBinUpEdge(b), hXjRw[i]->GetBinContent(b));
    printf("\n");
  }

  // ---------------- drawing ----------------
  drawer d("pythia", "nominal"); // ctor args only pick files for other helpers; labels come from drawAll
  string rLabel = Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir]);
  string phoLabel = Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", phoLow, phoHigh);
  string outdir = ana::path("claude_checks/jet_timing/pdfs");
  gSystem->mkdir(outdir.c_str(), true);
  string pdf = outdir + Form("/draw_jet_timing_r%02d.pdf", (int)(ana::JetRs[ir]*10+0.5));
  TCanvas * c = new TCanvas("c", "", 700, 600);
  auto pad = [&]() { c->Clear(); c->SetLeftMargin(0.15); c->SetBottomMargin(0.13); c->SetRightMargin(0.05); c->SetTopMargin(0.05); c->SetLogy(0); c->SetLogz(0); c->SetRightMargin(0.05); };
  auto legend = [](double x1, double y1, double x2, double y2) { TLegend * l = new TLegend(x1,y1,x2,y2); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.03); return l; };
  int cols[8] = {kBlack, kBlue+1, kAzure+7, kGreen+2, kOrange+1, kRed+1, kMagenta+1, kGray+2};
  c->SaveAs((pdf+"[").c_str());

  // p1
  pad();
  unitNorm(hDtJetA); unitNorm(hDtClA); sanitize(hDtJetA); sanitize(hDtClA);
  hDtClA->SetLineColor(kBlue+1); hDtJetA->SetLineColor(kRed+1);
  hDtClA->SetLineWidth(2); hDtJetA->SetLineWidth(2);
  hDtClA->SetMaximum(1.5*std::max(hDtClA->GetMaximum(), hDtJetA->GetMaximum())); hDtClA->SetMinimum(0);
  hDtClA->Draw("hist"); hDtJetA->Draw("hist same");
  windowLines(0, hDtClA->GetMaximum());
  { TLegend * l = legend(0.58, 0.62, 0.93, 0.74);
    l->AddEntry(hDtClA, "photon cluster", "l"); l->AddEntry(hDtJetA, "recoil jet", "l"); l->Draw(); }
  d.drawAll({"p+p Run24 Data"}, {phoLabel, rLabel, "Region A, paired", "dashed: treemaking window (data only)"}, .18, .9, 14, gPad->GetWh()*0.8);
  c->SaveAs(pdf.c_str());

  // p2, p3: slices
  auto drawSlices = [&](TH1D ** hs, int nh, const double * edges, const char * fmt, const string & extra) {
    pad();
    double ymax = 0;
    for (int i = 0; i < nh; i++) { unitNorm(hs[i]); sanitize(hs[i]); ymax = std::max(ymax, hs[i]->GetMaximum()); }
    TLegend * l = legend(0.62, 0.55, 0.93, 0.9);
    for (int i = 0; i < nh; i++) {
      hs[i]->SetLineColor(cols[i+1]); hs[i]->SetLineWidth(2);
      hs[i]->SetMaximum(1.5*ymax); hs[i]->SetMinimum(0);
      hs[i]->Draw(i ? "hist same" : "hist");
      l->AddEntry(hs[i], Form(fmt, edges[i], edges[i+1], hs[i]->GetEntries()), "l");
    }
    windowLines(0, 1.5*ymax); l->Draw();
    d.drawAll({"p+p Run24 Data"}, {phoLabel, rLabel, "Regions A+B+C+D, paired", extra}, .18, .9, 14, gPad->GetWh()*0.8);
    c->SaveAs(pdf.c_str());
  };
  drawSlices(hDtJetPt, nJpt, jptEdges, "p_{T}^{jet} %.0f-%.0f GeV (%.0f)", "slices of corrected jet p_{T}");
  drawSlices(hDtJetEmf, nEmf, emfEdges, "EM frac. %.2f-%.2f (%.0f)", "slices of jet EM fraction");

  // p4: profiles
  for (TProfile * p : {pDtVsPt, pDtVsEmf}) {
    pad();
    p->SetMarkerStyle(20); p->SetMarkerColor(kBlack); p->SetLineColor(kBlack);
    p->SetMinimum(0); p->SetMaximum(5);
    p->Draw("e");
    TLine * l = new TLine(p->GetXaxis()->GetXmin(), 2, p->GetXaxis()->GetXmax(), 2); l->SetLineStyle(2); l->SetLineColor(kGray+2); l->Draw();
    d.drawAll({"p+p Run24 Data"}, {phoLabel, rLabel, "Regions A+B+C+D, paired", "dashed: window center (2 ns)", "mean of in-window entries only"}, .18, .9, 14, gPad->GetWh()*0.8);
    c->SaveAs(pdf.c_str());
  }

  // p5: jet vs cluster
  pad(); c->SetRightMargin(0.13);
  sanitize(hJetVsCl);
  hJetVsCl->Draw("colz");
  windowLines(-1, 5);
  for (double y : {tLow, tHigh}) { TLine * l = new TLine(-1, y, 5, y); l->SetLineStyle(2); l->SetLineColor(kGray+2); l->Draw(); }
  d.drawAll({"p+p Run24 Data"}, {phoLabel, rLabel, "Region A, paired", Form("corr. = %.2f", hJetVsCl->GetCorrelationFactor())}, .18, .9, 14, gPad->GetWh()*0.8);
  c->SaveAs(pdf.c_str());

  // p6: dphi edge vs core
  pad();
  for (TH1D * h : {hDphiEdge, hDphiCore, hDphiHigh}) { unitNorm(h); sanitize(h); h->SetLineWidth(2); }
  hDphiCore->SetLineColor(kBlack); hDphiEdge->SetLineColor(kRed+1); hDphiHigh->SetLineColor(kBlue+1);
  c->SetLogy(1);
  hDphiCore->SetMinimum(0.005); hDphiCore->SetMaximum(50*std::max(hDphiCore->GetMaximum(), hDphiEdge->GetMaximum()));
  hDphiCore->Draw("hist"); hDphiEdge->Draw("hist same"); hDphiHigh->Draw("hist same");
  { TLine * l = new TLine(ana::oppcut, 0.005, ana::oppcut, hDphiCore->GetMaximum()); l->SetLineStyle(2); l->SetLineColor(kGray+2); l->Draw();
    TLegend * lg = legend(0.18, 0.45, 0.6, 0.6);
    lg->AddEntry(hDphiEdge, Form("0 < #Deltat_{jet} < 0.5 ns (%.0f)", hDphiEdge->GetEntries()), "l");
    lg->AddEntry(hDphiCore, Form("1 < #Deltat_{jet} < 3 ns (%.0f)", hDphiCore->GetEntries()), "l");
    lg->AddEntry(hDphiHigh, Form("3.5 < #Deltat_{jet} < 4 ns (%.0f)", hDphiHigh->GetEntries()), "l");
    lg->Draw(); }
  d.drawAll({"p+p Run24 Data"}, {phoLabel, rLabel, "Region A, no #Delta#phi cut", "dashed: 7#pi/8"}, .18, .9, 14, gPad->GetWh()*0.8);
  c->SaveAs(pdf.c_str());

  // p7: xJ edge vs core
  pad();
  for (TH1D * h : {hXjEdge, hXjCore}) { unitNorm(h); sanitize(h); h->SetLineWidth(2); h->SetMarkerStyle(20); }
  hXjCore->SetLineColor(kBlack); hXjCore->SetMarkerColor(kBlack);
  hXjEdge->SetLineColor(kRed+1); hXjEdge->SetMarkerColor(kRed+1); hXjEdge->SetMarkerStyle(24);
  hXjCore->GetXaxis()->SetRangeUser(0, 2);
  hXjCore->SetMinimum(0); hXjCore->SetMaximum(1.6*std::max(hXjCore->GetMaximum(), hXjEdge->GetMaximum()));
  hXjCore->Draw("e"); hXjEdge->Draw("e same");
  { TLegend * lg = legend(0.6, 0.6, 0.93, 0.72);
    lg->AddEntry(hXjEdge, Form("0 < #Deltat_{jet} < 0.5 ns (%.0f)", hXjEdge->GetEntries()), "lp");
    lg->AddEntry(hXjCore, Form("1 < #Deltat_{jet} < 3 ns (%.0f)", hXjCore->GetEntries()), "lp");
    lg->Draw(); }
  d.drawAll({"p+p Run24 Data"}, {phoLabel, rLabel, "Region A, paired", "reco level, not purity corrected", Form("KS prob. = %.2f", hXjEdge->KolmogorovTest(hXjCore))}, .18, .9, 14, gPad->GetWh()*0.8);
  c->SaveAs(pdf.c_str());

  // p8: loss estimates
  for (int k = 0; k < 2; k++) {
    TGraphErrors * g = k ? gLossEmf : gLossPt;
    pad();
    TH1D * fr = new TH1D(Form("frLoss%d",k), k ? ";jet EM fraction;fraction outside [0,4] ns (Gaussian model)" :
                                           ";corrected jet p_{T} [GeV];fraction outside [0,4] ns (Gaussian model)",
                         1, k ? 0 : jptEdges[0], k ? 1 : jptEdges[nJpt]);
    fr->SetMinimum(0); fr->SetMaximum(0.5); fr->Draw();
    g->SetMarkerStyle(20); g->SetMarkerColor(kRed+1); g->SetLineColor(kRed+1); g->Draw("p same");
    d.drawAll({"p+p Run24 Data"}, {phoLabel, rLabel, "Regions A+B+C+D, paired",
              "Gaussian fit inside window, integrated outside", "valid only if edge entries are in-time jets"}, .18, .9, 14, gPad->GetWh()*0.8);
    c->SaveAs(pdf.c_str());
  }
  // p9
  pad();
  { TH1D * fr = new TH1D("frRw", ";x_{J#gamma};reweighted / nominal (shape)", 1, 0.3, 1.7);
    fr->SetMinimum(0.7); fr->SetMaximum(1.6); fr->Draw();
    TLine * one = new TLine(0.3, 1, 1.7, 1); one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();
    TLegend * lg = legend(0.55, 0.2, 0.93, 0.35);
    int mk[3] = {20, 21, 22};
    for (int i = 0; i < ana::nPtBinsUsed; i++) {
      sanitize(hXjRw[i]);
      // ratio of two shapes from the same events - its statistical error is not the
      // naive one; drawn without errors, it only indicates the size of the effect
      for (int b = 0; b <= hXjRw[i]->GetNbinsX()+1; b++) { hXjRw[i]->SetBinError(b, 0); if (hXjNom[i]->GetBinContent(b) <= 0) hXjRw[i]->SetBinContent(b, -1); }
      hXjRw[i]->SetMarkerStyle(mk[i]); hXjRw[i]->SetMarkerColor(cols[i+1]); hXjRw[i]->SetLineColor(cols[i+1]);
      hXjRw[i]->Draw("p same");
      lg->AddEntry(hXjRw[i], Form("%.0f < p_{T}^{#gamma} < %.0f GeV", ana::ptBinsUsed[i], ana::ptBinsUsed[i+1]), "p");
    }
    lg->Draw(); }
  d.drawAll({"p+p Run24 Data"}, {rLabel, "Region A, paired, reco level",
            "weight 1/(1 - loss(p_{T}^{jet})), Gaussian model (p8)", "lost jet treated as lost event"}, .18, .9, 14, gPad->GetWh()*0.8);
  c->SaveAs(pdf.c_str());
  c->SaveAs((pdf+"]").c_str());

  TFile * fout = TFile::Open((outdir + Form("/draw_jet_timing_r%02d.root", (int)(ana::JetRs[ir]*10+0.5))).c_str(), "recreate");
  for (TObject * o : std::vector<TObject*>{hDtJetA, hDtClA, pDtVsPt, pDtVsEmf, hJetVsCl, hDphiEdge, hDphiCore, hDphiHigh, hXjEdge, hXjCore, gLossPt, gLossEmf}) o->Write();
  for (int i = 0; i < nJpt; i++) hDtJetPt[i]->Write();
  for (int i = 0; i < nEmf; i++) hDtJetEmf[i]->Write();
  gLossPt->Write("gLossPt"); gLossEmf->Write("gLossEmf");
  fout->Close();
}
