#include "../../src/ana.h"
#include "../../src/pho_object.h"
#include "../../src/jet_object.h"
#include "../../src/insitu_utility.h"
// Explicit load - see drawing/draw_final_result.C. Run interpreted, never with ACLiC "+".
R__LOAD_LIBRARY(libgammajet_unfold.so);

// How the Data three-jet veto drives the R = 0.2 threejet in-situ result to the scan edge.
// Only the Data side depends on the assumed p_a (the threejet systag's own table value): the recoil jet is
// jet_pt_calib/p_a and the veto is thirdjet_pt/p_a > cut, i.e. a calibrated third-jet threshold of cut x p_a.
// The trees store Data jets fully only above calibrated 4.5 GeV (CaloAna storage rule, calib/0.90 >= 5 GeV).
// Here the Data selection is redone at each assumed p_a (pairing as unfolder::check_pair, x_J floor, region
// from ana::findabcdBin), R = 0.2.
//   p1  Data third-jet calibrated pT (paired at p_a = 0.90), incomplete band and veto thresholds
//   p2  fraction of paired events vetoed vs assumed p_a, per ABCD region (5 GeV veto)
//   p3  simple ABCD purity 1 - N_B N_C / (N_A N_D), 25-35 GeV, vs assumed p_a, per veto variant. No MC leakage
//       correction (puritymaker has it), so a trend, not the pipeline's value; the pipeline's P_A at each
//       run's start p_a is overlaid
//   p4  the pipeline runs: purity-corrected threejet R = 0.2 result vs the p_a it started from

namespace {
  struct Ev { int k; float phoPt, phoEta, phoPhi, jetCalib, jetEta, jetPhi, t3Calib, t3dr; };
  const int ir = 0; // R = 0.2
  const int nVar = 4;
  const char * varName[nVar] = {"no veto", "5 GeV veto", "5 GeV veto, #DeltaR(#gamma,j_{3}) #geq 0.4", "7 GeV veto"};
  const int varCol[nVar] = {kGray+2, kRed+1, kOrange+1, kBlue+1};
  const char * reg[4] = {"A (iso, tight)", "B (non-iso, tight)", "C (iso, non-tight)", "D (non-iso, non-tight)"};
  const int regCol[4] = {kBlack, kRed+1, kBlue+1, kMagenta+1};

  bool paired(const Ev & e, float pa) {
    float jetPt = e.jetCalib / pa;
    if (jetPt <= ana::jet_calib_pt_cut[ir]) return false;
    int ptbin = ana::findPtBin(e.phoPt);
    float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(ana::jet_calib_pt_cut[ir]/ana::ptBins[ptbin])+1];
    return jetPt / e.phoPt >= lowbin;
  }
  bool vetoed(const Ev & e, float pa, int var) {
    if (var == 0) return false;
    float cut = (var == 3) ? 7.0 : 5.0;
    if (var == 2 && e.t3dr < 0.4) return false;
    return e.t3Calib / pa > cut;
  }
  void label(vector<string> f, TPad * p, double x = .18, double y = .9) {
    insitu_utility::drawSPhenixLabel({"p+p Run24 Data"}, f, x, y, 14, p->GetWh());
  }
}

void threejet_runaway() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();

  TFile * fin = TFile::Open(ana::path("trees/gammajet_Data.root"), "read");
  TTree * t = (TTree*)fin->Get("towerntup");
  float vz, cluster_pt, cluster_e, cluster_eta, cluster_phi, cluster_time;
  float cluster_showershape[12], cluster_bdt_scores[11];
  float jet_pt_calib[7], jet_e[7], jet_eta[7], jet_phi[7], jet_emfrac[7], jet_time[7];
  float thirdjet_pt[7], thirdjet_eta[7], thirdjet_phi[7];
  t->SetBranchStatus("*", 0);
  for (const char * b : {"vz","cluster_pt","cluster_e","cluster_eta","cluster_phi","cluster_time","cluster_showershape",
        "cluster_bdt_scores","jet_pt_calib","jet_e","jet_eta","jet_phi","jet_emfrac","jet_time",
        "thirdjet_pt","thirdjet_eta","thirdjet_phi"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("vz", &vz);
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
  t->SetBranchAddress("thirdjet_pt", thirdjet_pt);
  t->SetBranchAddress("thirdjet_eta", thirdjet_eta);
  t->SetBranchAddress("thirdjet_phi", thirdjet_phi);

  // Cache the p_a-independent part of the selection.
  vector<Ev> evs;
  for (Long64_t e = 0; e < t->GetEntries(); e++) {
    t->GetEntry(e);
    if (fabs(vz) > ana::vzcut) continue;
    if (cluster_pt < 15 || cluster_pt >= 35) continue;
    if (jet_pt_calib[ir] <= 0) continue;
    pho_object pho(cluster_pt, cluster_e, cluster_eta, cluster_phi,
        cluster_showershape[10], cluster_showershape[11], cluster_time,
        cluster_bdt_scores[9], pho_object::get_showershape(cluster_showershape, cluster_pt));
    int k = ana::findabcdBin(pho.iso4, pho.bdt, 0);
    if (k < 0 || fabs(pho.eta) > ana::photonEtaCut) continue;
    jet_object jet(jet_pt_calib[ir], jet_e[ir], jet_eta[ir], jet_phi[ir], jet_emfrac[ir], 0, 0, jet_time[ir]);
    if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) continue;
    if (jet.deltaPhi(pho) < ana::oppcut) continue;
    float deta = thirdjet_eta[ir] - pho.eta, dphi = fabs(thirdjet_phi[ir] - pho.phi);
    if (dphi > M_PI) dphi = 2*M_PI - dphi;
    evs.push_back({k, pho.pt, pho.eta, pho.phi, jet_pt_calib[ir], jet.eta, jet.phi, thirdjet_pt[ir], (float)sqrt(deta*deta + dphi*dphi)});
  }
  printf("cached %zu photon-jet candidates (R=0.2, 15-35 GeV, regions A-D)\n", evs.size());

  // p_a scan
  const int nPa = 33; const double paLo = 0.80, paStep = 0.005;
  TGraph * gVeto[2][4];   // [all 15-35, 25-35][region], 5 GeV veto
  TGraphErrors * gPur[nVar];
  for (int s = 0; s < 2; s++) for (int k = 0; k < 4; k++) gVeto[s][k] = new TGraph();
  for (int v = 0; v < nVar; v++) gPur[v] = new TGraphErrors();
  printf("p_a   N(A,B,C,D) 25-35 GeV after 5 GeV veto   simple P_A: none / 5 / 5+dR / 7\n");
  for (int ip = 0; ip < nPa; ip++) {
    float pa = paLo + ip*paStep;
    double nP[2][4] = {}, nV[2][4] = {};
    double n[nVar][4] = {};
    for (auto & e : evs) {
      if (!paired(e, pa)) continue;
      bool hi = e.phoPt >= 25;
      for (int s = 0; s < 2; s++) {
        if (s == 1 && !hi) continue;
        nP[s][e.k]++;
        if (vetoed(e, pa, 1)) nV[s][e.k]++;
      }
      if (!hi) continue;
      for (int v = 0; v < nVar; v++) if (!vetoed(e, pa, v)) n[v][e.k]++;
    }
    for (int s = 0; s < 2; s++) for (int k = 0; k < 4; k++) gVeto[s][k]->SetPoint(ip, pa, nP[s][k] > 0 ? nV[s][k]/nP[s][k] : 0);
    double pur[nVar];
    for (int v = 0; v < nVar; v++) {
      double A = n[v][0], B = n[v][1], C = n[v][2], D = n[v][3];
      pur[v] = (A > 0 && D > 0) ? 1 - B*C/(A*D) : 0;
      // Poisson error on B C / (A D)
      double r = (A > 0 && D > 0) ? B*C/(A*D) : 0;
      double rel = (A > 0 && B > 0 && C > 0 && D > 0) ? sqrt(1/A + 1/B + 1/C + 1/D) : 0;
      gPur[v]->SetPoint(ip, pa, pur[v]);
      gPur[v]->SetPointError(ip, 0, r*rel);
    }
    if (ip % 4 == 0)
      printf("%.3f  %4.0f %3.0f %4.0f %3.0f      %.3f / %.3f / %.3f / %.3f\n", pa, n[1][0], n[1][1], n[1][2], n[1][3], pur[0], pur[1], pur[2], pur[3]);
  }

  string outdir = ana::path("claude_checks/insitu_variations/pdfs");
  gSystem->mkdir(outdir.c_str(), true);
  string pdf = outdir + "/threejet_runaway.pdf";
  TCanvas * c = new TCanvas("c", "", 900, 700);
  auto pad = [&]() { c->Clear(); c->SetLeftMargin(0.13); c->SetBottomMargin(0.12); c->SetRightMargin(0.04); c->SetTopMargin(0.05); c->SetLogy(0); };
  auto legend = [](double x1, double y1, double x2, double y2) { TLegend * l = new TLegend(x1, y1, x2, y2); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.032); return l; };

  // ---- p1: third-jet calibrated pT ----
  {
    pad();
    TH1D * h = new TH1D("h3", ";third-jet calibrated p_{T} [GeV];paired events / 0.25 GeV", 40, 2, 12);
    for (auto & e : evs) if (paired(e, 0.90) && e.t3Calib > 0) h->Fill(e.t3Calib);
    h->SetLineColor(kBlack); h->SetLineWidth(2);
    h->SetMaximum(1.6*h->GetMaximum()); h->SetMinimum(0);
    h->Draw("hist");
    TBox * band = new TBox(4.0, 0, 4.5, h->GetMaximum());
    band->SetFillColor(kGray); band->SetFillStyle(1001); band->Draw();
    h->Draw("hist same");
    struct L { double x; int col; const char * txt; };
    vector<L> lines = {{5*0.9013, kRed+1, "5 GeV #times 0.901"}, {5*0.868, kOrange+1, "5 GeV #times 0.868"},
                       {5*0.80, kMagenta+1, "5 GeV #times 0.80"}, {7*0.80, kBlue+1, "7 GeV #times 0.80"}};
    TLegend * l = legend(0.5, 0.6, 0.95, 0.9); l->SetTextSize(0.027);
    l->AddEntry(band, "Data storage incomplete (4.0-4.5 GeV)", "f");
    for (auto & x : lines) {
      TLine * ln = new TLine(x.x, 0, x.x, h->GetMaximum()); ln->SetLineColor(x.col); ln->SetLineWidth(3); ln->SetLineStyle(2); ln->Draw();
      l->AddEntry(ln, Form("veto threshold %s", x.txt), "l");
    }
    l->Draw();
    label({"15 GeV < p_{T}^{#gamma} < 35 GeV", "Jet R=0.2, paired at p_{a} = 0.90, regions A-D",
           "veto: calibrated p_{T,3} > cut #times p_{a}"}, (TPad*)gPad, .62, .55);
    c->SaveAs((pdf + "(").c_str());
  }

  // ---- p2: veto fraction vs p_a ----
  for (int s = 0; s < 2; s++) {
    pad();
    TH1D * fr = new TH1D(Form("fr2_%d", s), ";assumed p_{a} (threejet);fraction of paired events vetoed (5 GeV)", 1, 0.79, 0.97);
    fr->SetMinimum(0); fr->SetMaximum(1.1); fr->Draw("axis");
    TBox * band = new TBox(0.79, 0, 0.90, 1.1); band->SetFillColor(kGray); band->SetFillStyle(1001); band->Draw();
    fr->Draw("axis same");
    TLegend * l = legend(0.55, 0.62, 0.95, 0.85);
    for (int k = 0; k < 4; k++) {
      gVeto[s][k]->SetLineColor(regCol[k]); gVeto[s][k]->SetMarkerColor(regCol[k]); gVeto[s][k]->SetMarkerStyle(20); gVeto[s][k]->SetMarkerSize(0.7); gVeto[s][k]->SetLineWidth(2);
      gVeto[s][k]->Draw("lp same");
      l->AddEntry(gVeto[s][k], reg[k], "lp");
    }
    l->Draw();
    label({s == 0 ? "15 GeV < p_{T}^{#gamma} < 35 GeV" : "25 GeV < p_{T}^{#gamma} < 35 GeV", "Jet R=0.2, paired at each p_{a}",
           "grey: veto threshold below 4.5 GeV (storage incomplete)"}, (TPad*)gPad);
    c->SaveAs(pdf.c_str());
  }

  // ---- p3: simple ABCD purity vs p_a ----
  {
    pad();
    TH1D * fr = new TH1D("fr3", ";assumed p_{a} (threejet);simple ABCD purity, 25-35 GeV", 1, 0.79, 0.97);
    fr->SetMinimum(0.3); fr->SetMaximum(1.5); fr->Draw("axis");
    TBox * band = new TBox(0.79, 0.3, 0.90, 1.5); band->SetFillColor(kGray); band->SetFillStyle(1001); band->Draw();
    fr->Draw("axis same");
    TLegend * l = legend(0.5, 0.68, 0.95, 0.94); l->SetTextSize(0.026); l->SetNColumns(2);
    for (int v = 0; v < nVar; v++) {
      gPur[v]->SetLineColor(varCol[v]); gPur[v]->SetMarkerColor(varCol[v]); gPur[v]->SetMarkerStyle(20); gPur[v]->SetMarkerSize(0.6); gPur[v]->SetLineWidth(2);
      gPur[v]->Draw("lp same");
      l->AddEntry(gPur[v], varName[v], "lp");
    }
    // pipeline (puritymaker, with MC leakage) P_A in 25-35 GeV at each run's start p_a
    TGraph * gPipe5 = new TGraph(); gPipe5->SetPoint(0, 0.9013, 0.711); gPipe5->SetPoint(1, 0.8073, 0.664); gPipe5->SetPoint(2, 0.800, 0.664);
    TGraph * gPipeDr = new TGraph(); gPipeDr->SetPoint(0, 0.9013, 0.737); gPipeDr->SetPoint(1, 0.868, 0.673);
    TGraph * gPipe7 = new TGraph(); gPipe7->SetPoint(0, 0.800, 0.798);
    TGraph * gNom = new TGraph(); gNom->SetPoint(0, 0.9013, 0.833);
    struct P { TGraph * g; int col; int sty; const char * txt; };
    for (auto & p : vector<P>{{gNom, kGray+2, 29, "pipeline P_{A}: no veto (nominal)"}, {gPipe5, kRed+1, 29, "pipeline P_{A}: 5 GeV veto"},
                              {gPipeDr, kOrange+1, 29, "pipeline P_{A}: 5 GeV, #DeltaR #geq 0.4"}, {gPipe7, kBlue+1, 29, "pipeline P_{A}: 7 GeV veto"}}) {
      p.g->SetMarkerStyle(p.sty); p.g->SetMarkerSize(2.2); p.g->SetMarkerColor(p.col); p.g->Draw("p same");
      l->AddEntry(p.g, p.txt, "p");
    }
    l->Draw();
    label({"25 GeV < p_{T}^{#gamma} < 35 GeV, Jet R=0.2", "lines: 1 - N_{B}N_{C}/(N_{A}N_{D}), no MC leakage",
           "stars: puritymaker P_{A} at the run's start p_{a}"}, (TPad*)gPad, .16, .9);
    c->SaveAs(pdf.c_str());
  }

  // ---- p4: pipeline iterations ----
  {
    pad();
    TH1D * fr = new TH1D("fr4", ";p_{a} the run started from (threejet R=0.2);purity-corrected result", 1, 0.79, 0.93);
    fr->SetMinimum(0.79); fr->SetMaximum(0.95); fr->Draw("axis");
    TLine * diag = new TLine(0.79, 0.79, 0.93, 0.93); diag->SetLineStyle(2); diag->Draw();
    TLine * edge = new TLine(0.79, 0.80, 0.93, 0.80); edge->SetLineStyle(3); edge->SetLineColor(kGray+2); edge->Draw();
    TLine * comp = new TLine(0.90, 0.79, 0.90, 0.95); comp->SetLineStyle(3); comp->SetLineColor(kGray+2); comp->Draw();
    struct Run { vector<pair<double,double>> pts; int col; const char * txt; };
    vector<Run> runs = {
      {{{0.9013, 0.8073}, {0.8073, 0.800}}, kRed+1, "5 GeV veto (start 0.901)"},
      {{{0.800, 0.800}}, kRed+3, "5 GeV veto, full pipeline (3 scans)"},
      {{{0.9013, 0.868}, {0.868, 0.800}}, kOrange+1, "5 GeV veto, #DeltaR #geq 0.4"},
      {{{0.800, 0.869}}, kGreen+2, "5 GeV veto + timing cut"},
      {{{0.800, 0.9055}}, kBlue+1, "7 GeV veto"},
      {{{0.800, 0.8982}}, kAzure+7, "9 GeV veto"},
    };
    TLegend * l = legend(0.45, 0.62, 0.95, 0.9);
    for (auto & r : runs) {
      TGraph * g = new TGraph();
      for (size_t i = 0; i < r.pts.size(); i++) g->SetPoint(i, r.pts[i].first, r.pts[i].second);
      g->SetMarkerStyle(20); g->SetMarkerSize(1.3); g->SetMarkerColor(r.col); g->SetLineColor(r.col); g->SetLineWidth(2);
      g->Draw(r.pts.size() > 1 ? "lp same" : "p same");
      for (size_t i = 0; i + 1 < r.pts.size(); i++) { // result -> next start
        TArrow * a = new TArrow(r.pts[i].second, r.pts[i].second, r.pts[i+1].first, r.pts[i+1].second, 0.015, "|>");
        a->SetLineColor(r.col); a->SetFillColor(r.col); a->SetLineStyle(2); a->Draw();
      }
      l->AddEntry(g, r.txt, "p");
    }
    l->Draw();
    label({"Jet R=0.2, threejet systag, one scan per point", "dashed: result = start; dotted: scan edge 0.80,",
           "and p_{a} = 0.90 (5 GeV veto reaches 4.5 GeV)"}, (TPad*)gPad, .16, .9);
    c->SaveAs((pdf + ")").c_str());
  }
  printf("Wrote %s\n", pdf.c_str());
}
