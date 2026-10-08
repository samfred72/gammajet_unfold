#include "../../src/ana.h"
#include "../../src/pho_object.h"
#include "../../src/jet_object.h"
#include "../../src/insitu_utility.h"
// Explicit load - see drawing/draw_final_result.C. Run interpreted, never with ACLiC "+".
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Why the three-jet variation fails at R = 0.2: photon - third-jet Delta R by ABCD region, in data.
// The treemaker drops jets within Delta R < R of the photon cluster (CaloAna.cc jet loop), while the
// photon isolation (iso_topo_04) uses a 0.4 cone. At R = 0.2 a jet at Delta R = 0.2-0.4 can be the
// third jet. Regions (ana::findabcdBin): A iso+tight, B non-iso+tight, C iso+non-tight, D non-iso+non-tight.
// Pairing as unfolder::check_pair: |v_z| cut, photon |eta| < photonEtaCut, jet |eta| < etacut - R,
// jet pT = jet_pt_calib/jesNominal > jet_calib_pt_cut, Delta phi > oppcut, x_J floor. Veto as the threejet
// systag on Data: thirdjet_pt / p_a > thirdJetPtCut, with p_a = jesNominal (not the threejet table value).

void threejet_dr() {
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

  const int nR = 2; const int irs[nR] = {0, 2};          // R = 0.2, 0.4
  const int nPho = 2; const double phoLo[nPho] = {15, 25}, phoHi[nPho] = {35, 35};
  const char * reg[4] = {"A (iso, tight)", "B (non-iso, tight)", "C (iso, non-tight)", "D (non-iso, non-tight)"};
  const int cols[4] = {kBlack, kRed+1, kBlue+1, kMagenta+1};
  TH1D * hdr[nR][nPho][4];
  double nPair[nR][nPho][4] = {}, nVeto[nR][nPho][4] = {}, nVetoNear[nR][nPho][4] = {};
  for (int r = 0; r < nR; r++) for (int p = 0; p < nPho; p++) for (int k = 0; k < 4; k++)
    hdr[r][p][k] = new TH1D(Form("hdr_%d_%d_%d", r, p, k), ";#DeltaR(#gamma, third jet);vetoing third jets / paired event", 30, 0, 3);

  Long64_t n = t->GetEntries();
  for (Long64_t e = 0; e < n; e++) {
    t->GetEntry(e);
    if (fabs(vz) > ana::vzcut) continue;
    if (cluster_pt < 15 || cluster_pt >= 35) continue;
    pho_object pho(cluster_pt, cluster_e, cluster_eta, cluster_phi,
        cluster_showershape[10], cluster_showershape[11], cluster_time,
        cluster_bdt_scores[9], pho_object::get_showershape(cluster_showershape, cluster_pt));
    int k = ana::findabcdBin(pho.iso4, pho.bdt, 0);
    if (k < 0) continue;
    if (fabs(pho.eta) > ana::photonEtaCut) continue;
    for (int r = 0; r < nR; r++) {
      int ir = irs[r];
      if (jet_pt_calib[ir] <= 0) continue;
      float jetPt = jet_pt_calib[ir] / ana::jesNominal[ir];
      if (jetPt <= ana::jet_calib_pt_cut[ir]) continue;
      jet_object jet(jetPt, jet_e[ir], jet_eta[ir], jet_phi[ir], jet_emfrac[ir], 0, 0, jet_time[ir]);
      if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) continue;
      if (jet.deltaPhi(pho) < ana::oppcut) continue;
      int ptbin = ana::findPtBin(pho.pt);
      float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(ana::jet_calib_pt_cut[ir]/ana::ptBins[ptbin])+1];
      if (jetPt / pho.pt < lowbin) continue;
      bool veto = thirdjet_pt[ir] / ana::jesNominal[ir] > ana::thirdJetPtCut;
      float deta = thirdjet_eta[ir] - pho.eta, dphi = fabs(thirdjet_phi[ir] - pho.phi);
      if (dphi > M_PI) dphi = 2*M_PI - dphi;
      float dr = sqrt(deta*deta + dphi*dphi);
      for (int p = 0; p < nPho; p++) {
        if (pho.pt < phoLo[p] || pho.pt >= phoHi[p]) continue;
        nPair[r][p][k]++;
        if (!veto) continue;
        nVeto[r][p][k]++;
        if (dr < 0.4) nVetoNear[r][p][k]++;
        hdr[r][p][k]->Fill(dr);
      }
    }
  }

  printf("Paired events vetoed by the three-jet cut (third jet / p_a > %.0f GeV), and the part with DeltaR(photon, third jet) < 0.4\n", ana::thirdJetPtCut);
  for (int r = 0; r < nR; r++) for (int p = 0; p < nPho; p++) {
    printf("R=%.1f, %.0f-%.0f GeV:\n", ana::JetRs[irs[r]], phoLo[p], phoHi[p]);
    for (int k = 0; k < 4; k++)
      printf("  %-24s pairs %6.0f  vetoed %5.1f%%  vetoed with dR<0.4 %5.1f%%\n", reg[k], nPair[r][p][k],
             100*nVeto[r][p][k]/std::max(nPair[r][p][k], 1.), 100*nVetoNear[r][p][k]/std::max(nPair[r][p][k], 1.));
  }

  string outdir = ana::path("claude_checks/insitu_variations/pdfs");
  gSystem->mkdir(outdir.c_str(), true);
  string pdf = outdir + "/threejet_dr.pdf";
  TCanvas * c = new TCanvas("c", "", 1200, 1000);
  c->Divide(2, 2);
  for (int r = 0; r < nR; r++) for (int p = 0; p < nPho; p++) {
    TPad * pad = (TPad*)c->cd(1 + 2*r + p);
    pad->SetLeftMargin(0.15); pad->SetBottomMargin(0.13); pad->SetRightMargin(0.04); pad->SetTopMargin(0.05);
    double ymax = 0;
    for (int k = 0; k < 4; k++) { if (nPair[r][p][k] > 0) hdr[r][p][k]->Scale(1.0/nPair[r][p][k]); ymax = std::max(ymax, hdr[r][p][k]->GetMaximum()); }
    TLegend * l = new TLegend(0.5, 0.55, 0.95, 0.75); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.035);
    for (int k = 0; k < 4; k++) {
      TH1D * h = hdr[r][p][k];
      h->SetLineColor(cols[k]); h->SetMarkerColor(cols[k]); h->SetMarkerStyle(20); h->SetMarkerSize(0.6); h->SetLineWidth(2);
      h->SetMaximum(1.7*ymax); h->SetMinimum(0);
      h->Draw(k == 0 ? "e" : "e same");
      l->AddEntry(h, reg[k], "lp");
    }
    l->Draw();
    TLine * lr = new TLine(ana::JetRs[irs[r]], 0, ana::JetRs[irs[r]], 1.7*ymax); lr->SetLineStyle(2); lr->Draw();
    TLine * l4 = new TLine(0.4, 0, 0.4, 1.7*ymax); l4->SetLineStyle(3); l4->SetLineColor(kGray+2); l4->Draw();
    insitu_utility::drawSPhenixLabel({"p+p Run24 Data"},
        {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", phoLo[p], phoHi[p]), Form("Jet R=%.1f, paired", ana::JetRs[irs[r]]),
         Form("third jet / p_{a} > %.0f GeV", ana::thirdJetPtCut), "dashed: #DeltaR = R, dotted: 0.4 (iso cone)"}, .19, .9, 14, pad->GetWh());
  }
  c->SaveAs(pdf.c_str());
  printf("Wrote %s\n", pdf.c_str());
}
