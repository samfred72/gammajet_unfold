#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/pho_object.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/jet_object.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/reweight_utility.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/unfold_utility.h"
// Explicit load; run interpreted, never with ACLiC "+" (sibling libgammajet.so collision).
R__LOAD_LIBRARY(libgammajet_unfold.so);

// PPG18 review round 2, issue 3 (in-situ scans stuck at the 0.900 scan edge at R = 0.2/0.3).
//
// The purity-corrected in-situ mean (insitu_utility::computeCorrectedMeans) and the
// unfolding input (unfold_utility::purityCorrect) both assume
//     A = P_A S + (1-P_A) B,   C = P_C S + (1-P_C) B
// with ONE background xJ shape B shared by regions A and C, and one signal shape S. This
// macro checks that assumption in the QCD Pythia8 "Jet12_long" sample, where every reco
// cluster can be classified with truth:
//   signal     = truth-tagged reco photon: reco cluster within dR < 0.1 of the tree's
//                truth photon (unfolder::check_match) - the purity method's own signal
//                definition (unfolder.cc hclusterpt_abcd_truthmatched)
//   background = every other cluster (decay photons, merged pi0s, hadrons)
// Selection/kinematics follow unfolder.cc's nominal MC branch: reco photon pT =
// cluster_pt + N(0, emResolutionSigma(truth pT) * truth pT), jet pT = jet_pt_smear_truth[ir],
// check_pair cuts incl. the xJ floor, vz/cluster-pT mcWeight, Jet12_long truth-jet window.
//
// Output per R (0.2, 0.3, 0.4 for reference) and reported photon-pT bin:
//   page 1..3  normalized xJ: S_A, B_A, B_C, (and C total), with B_C/B_A and S_A/B_A ratios
//   page 4     <xJ> of each component vs photon pT, per R
//   log        <xJ> table and the bias the B_A != B_C mismatch alone induces in the
//              purity-corrected mean, using the DATA purities (hists/purity_nominal.root)
//              for that R: build A_mix = P_A S_A + (1-P_A) B_A and C_mix = P_C S_A +
//              (1-P_C) B_C from the MC shapes, apply the analysis coefficients
//              (purityCorrectCoeffs), compare to <S_A>. Since the in-situ fit scales jet
//              pT by 1/p_a, a mean bias delta translates to delta p_a / p_a ~ delta/<xJ>.

namespace {
  void sanitize(TH1 * h) {
    for (int b = 0; b <= h->GetNcells(); b++)
      if (!std::isfinite(h->GetBinContent(b)) || !std::isfinite(h->GetBinError(b))) { h->SetBinContent(b, 0); h->SetBinError(b, 0); }
  }
  TH1D * shape(TH1D * h, const char * name) {
    TH1D * s = (TH1D*)h->Clone(name);
    if (s->Integral() > 0) s->Scale(1.0/s->Integral(), "width");
    sanitize(s);
    return s;
  }
  // weighted mean and its error from a filled histogram's own moments (same bins for all)
  struct Mean { double m, e, neff; };
  Mean meanOf(TH1D * h) {
    double sw = 0, sw2 = 0;
    for (int b = 1; b <= h->GetNbinsX(); b++) { sw += h->GetBinContent(b); sw2 += pow(h->GetBinError(b), 2); }
    double neff = sw2 > 0 ? sw*sw/sw2 : 0;
    return {h->GetMean(), neff > 1 ? h->GetStdDev()/sqrt(neff) : 0, neff};
  }
}

void draw_jet12_sig_bkg_xj()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  // Unbinned means inside TH1 (GetMean) use the stats accumulated at Fill time, which is
  // what we want - keep them even after Scale().
  TH1::StatOverflows(false);

  const int irs[3] = {0, 1, 2}; // R = 0.2, 0.3, 0.4
  const int nR = 3, nPt = ana::nPtBinsUsed;
  enum { kSA, kBA, kSC, kBC, kCall, kN };
  const char * cname[kN] = {"SA", "BA", "SC", "BC", "Call"};
  TH1D * h[nR][nPt][kN];
  for (int r = 0; r < nR; r++) for (int p = 0; p < nPt; p++) for (int k = 0; k < kN; k++)
    h[r][p][k] = new TH1D(Form("h_%s_r%d_pt%d", cname[k], irs[r], p), ";x_{J#gamma};(1/N) dN/dx_{J#gamma}", ana::nUnfoldXjBins, ana::unfoldXjBins);

  TFile * fin = TFile::Open("/home/samson72/sphnx/gammajet_unfold/trees/gammajet_pythia_Jet12_long.root", "read");
  TTree * t = (TTree*)fin->Get("towerntup");
  float vz, cluster_pt, cluster_e, cluster_eta, cluster_phi, cluster_time;
  float cluster_showershape[12], cluster_bdt_scores[11];
  float truth_cluster_pt, truth_cluster_e, truth_cluster_eta, truth_cluster_phi, truth_cluster_iso3, truth_cluster_iso4;
  float jet_e[7], jet_eta[7], jet_phi[7], jet_emfrac[7], jet_pt_smear_truth[7], truth_jet_pt[7];
  t->SetBranchStatus("*", 0);
  for (const char * b : {"vz","cluster_pt","cluster_e","cluster_eta","cluster_phi","cluster_time","cluster_showershape",
        "cluster_bdt_scores","truth_cluster_pt","truth_cluster_e","truth_cluster_eta","truth_cluster_phi",
        "truth_cluster_iso3","truth_cluster_iso4","jet_e","jet_eta","jet_phi","jet_emfrac","jet_pt_smear_truth","truth_jet_pt"})
    t->SetBranchStatus(b, 1);
  t->SetBranchAddress("vz", &vz);
  t->SetBranchAddress("cluster_pt", &cluster_pt);
  t->SetBranchAddress("cluster_e", &cluster_e);
  t->SetBranchAddress("cluster_eta", &cluster_eta);
  t->SetBranchAddress("cluster_phi", &cluster_phi);
  t->SetBranchAddress("cluster_time", &cluster_time);
  t->SetBranchAddress("cluster_showershape", cluster_showershape);
  t->SetBranchAddress("cluster_bdt_scores", cluster_bdt_scores);
  t->SetBranchAddress("truth_cluster_pt", &truth_cluster_pt);
  t->SetBranchAddress("truth_cluster_e", &truth_cluster_e);
  t->SetBranchAddress("truth_cluster_eta", &truth_cluster_eta);
  t->SetBranchAddress("truth_cluster_phi", &truth_cluster_phi);
  t->SetBranchAddress("truth_cluster_iso3", &truth_cluster_iso3);
  t->SetBranchAddress("truth_cluster_iso4", &truth_cluster_iso4);
  t->SetBranchAddress("jet_e", jet_e);
  t->SetBranchAddress("jet_eta", jet_eta);
  t->SetBranchAddress("jet_phi", jet_phi);
  t->SetBranchAddress("jet_emfrac", jet_emfrac);
  t->SetBranchAddress("jet_pt_smear_truth", jet_pt_smear_truth);
  t->SetBranchAddress("truth_jet_pt", truth_jet_pt);

  Reweighter rw;
  TRandom3 rand(18);
  Long64_t n = t->GetEntries();
  for (Long64_t e = 0; e < n; e++) {
    t->GetEntry(e);
    if (fabs(vz) > ana::vzcut) continue;
    float w = rw.GetWeight(vz, cluster_pt);
    float recoPt = cluster_pt + rand.Gaus(0, ana::emResolutionSigma(truth_cluster_pt, ana::emrNominal)*truth_cluster_pt);
    pho_object pho(recoPt, cluster_e, cluster_eta, cluster_phi, cluster_showershape[10], cluster_showershape[11],
        cluster_time, cluster_bdt_scores[9], pho_object::get_showershape(cluster_showershape, recoPt));
    int iabcd = ana::findabcdBin(pho.iso4, pho.bdt, 0);
    if (iabcd != 0 && iabcd != 2) continue; // only regions A and C enter the subtraction
    if (fabs(pho.eta) > ana::etacut) continue;
    int ipt = -1; for (int i = 0; i < nPt; i++) if (pho.pt >= ana::ptBinsUsed[i] && pho.pt < ana::ptBinsUsed[i+1]) ipt = i;
    if (ipt < 0) continue;

    bool isSignal = false;
    if (truth_cluster_pt > 0) {
      pho_object tru(truth_cluster_pt, truth_cluster_e, truth_cluster_eta, truth_cluster_phi,
          truth_cluster_iso3, truth_cluster_iso4, 0, 0.99, 2);
      isSignal = pho.deltaR(tru) < 0.1;
    }

    for (int r = 0; r < nR; r++) {
      int ir = irs[r];
      if (!(truth_jet_pt[ir] > 0 && truth_jet_pt[ir] < 100)) continue; // Jet12_long keep window (treeuser.h)
      if (jet_pt_smear_truth[ir] <= 0) continue;
      jet_object jet(jet_pt_smear_truth[ir], jet_e[ir], jet_eta[ir], jet_phi[ir], jet_emfrac[ir], 0, 0, 0);
      // unfolder::check_pair
      int ptbin = ana::findPtBin(pho.pt);
      float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(ana::jet_calib_pt_cut[ir]/ana::ptBins[ptbin])+1];
      float xj = jet.pt/pho.pt;
      if (xj < lowbin) continue;
      if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) continue;
      if (jet.deltaPhi(pho) < ana::oppcut) continue;

      if (iabcd == 0) h[r][ipt][isSignal ? kSA : kBA]->Fill(xj, w);
      else { h[r][ipt][isSignal ? kSC : kBC]->Fill(xj, w); h[r][ipt][kCall]->Fill(xj, w); }
    }
  }

  // ---------------- numbers ----------------
  FILE * out = stdout;
  TGraphErrors * gMean[nR][kN];
  for (int r = 0; r < nR; r++) for (int k = 0; k < kN; k++) gMean[r][k] = new TGraphErrors();
  for (int r = 0; r < nR; r++) {
    int ir = irs[r];
    fprintf(out, "\n===== R = %.1f =====\n", ana::JetRs[ir]);
    fprintf(out, "%-9s %8s %8s %8s %8s | %7s %7s %7s %7s %7s | %6s %6s | %9s %9s\n", "pT bin", "Neff SA", "Neff BA", "Neff SC", "Neff BC",
            "<S_A>", "<B_A>", "<S_C>", "<B_C>", "<C>", "P_A", "P_C", "bias", "dpa/pa");
    for (int p = 0; p < nPt; p++) {
      Mean m[kN]; for (int k = 0; k < kN; k++) m[k] = meanOf(h[r][p][k]);
      for (int k = 0; k < kN; k++) { gMean[r][k]->SetPoint(p, 0.5*(ana::ptBinsUsed[p]+ana::ptBinsUsed[p+1]) + (k-2)*0.3, m[k].m);
                                     gMean[r][k]->SetPointError(p, 0, m[k].e); }
      float lo = ana::ptBinsUsed[p], hi = ana::ptBinsUsed[p+1];
      float pA = ana::getPurity(lo, hi, "nominal", ir), pC = ana::getPurityC(lo, hi, "nominal", ir);
      // Mixture with the data purities; A and C normalized to unit count, so K = NA/NC = 1.
      // The signal shape in C is set equal to S_A so that ONLY the background-shape
      // mismatch B_A vs B_C enters the bias (S_A vs S_C is shown separately in the table).
      float cA, cC; unfold_utility::purityCorrectCoeffs(pA, pC, 1, 1, cA, cC);
      double mA = pA*m[kSA].m + (1-pA)*m[kBA].m;
      double mC = pC*m[kSA].m + (1-pC)*m[kBC].m;
      double corr = (cA*mA - cC*mC)/(cA - cC);
      double bias = corr - m[kSA].m;
      // error: B_A and B_C mean errors propagated through the linear combination
      double dcorr_dBA = cA*(1-pA)/(cA-cC), dcorr_dBC = -cC*(1-pC)/(cA-cC);
      double biasErr = hypot(dcorr_dBA*m[kBA].e, dcorr_dBC*m[kBC].e);
      fprintf(out, "%3.0f-%-5.0f %8.0f %8.0f %8.0f %8.0f | %7.3f %7.3f %7.3f %7.3f %7.3f | %6.3f %6.3f | %+6.3f+-%.3f %+6.1f%%\n",
              lo, hi, m[kSA].neff, m[kBA].neff, m[kSC].neff, m[kBC].neff,
              m[kSA].m, m[kBA].m, m[kSC].m, m[kBC].m, m[kCall].m, pA, pC, bias, biasErr, 100*bias/m[kSA].m);
      fprintf(out, "          mean errors: S_A %.3f  B_A %.3f  S_C %.3f  B_C %.3f\n", m[kSA].e, m[kBA].e, m[kSC].e, m[kBC].e);
      fprintf(out, "          KS(B_A,B_C)=%.3f  KS(S_A,B_A)=%.3f   MC signal fraction in A (Jet12 only, not physical): %.3f\n",
              h[r][p][kBA]->KolmogorovTest(h[r][p][kBC]), h[r][p][kSA]->KolmogorovTest(h[r][p][kBA]),
              h[r][p][kSA]->Integral()/std::max(1e-9, h[r][p][kSA]->Integral()+h[r][p][kBA]->Integral()));
    }
  }

  // ---------------- drawing ----------------
  drawer d("pythia", "nominal");
  string outdir = "/home/samson72/sphnx/gammajet_unfold/claude_checks/insitu_scan_edge/pdfs";
  gSystem->mkdir(outdir.c_str(), true);
  string pdf = outdir + "/draw_jet12_sig_bkg_xj.pdf";
  TCanvas * c = new TCanvas("c", "", 1500, 800);
  c->SaveAs((pdf+"[").c_str());
  int col[kN] = {kBlack, kRed+1, kGray+1, kBlue+1, kAzure+7};
  int mk[kN]  = {20, 21, 24, 22, 26};
  const char * lab[kN] = {"truth-tagged, region A (S_{A})", "background, region A (B_{A})", "truth-tagged, region C (S_{C})",
                          "background, region C (B_{C})", "region C, all"};
  const char * shortLab[kN] = {"S_{A} (tagged, A)", "B_{A} (bkg., A)", "S_{C} (tagged, C)", "B_{C} (bkg., C)", "C (all)"};
  for (int r = 0; r < nR; r++) {
    int ir = irs[r];
    c->Clear();
    for (int p = 0; p < nPt; p++) {
      c->cd();
      TPad * top = new TPad(Form("top%d_%d",r,p), "", p/3.0, 0.35, (p+1)/3.0, 1.0);
      TPad * bot = new TPad(Form("bot%d_%d",r,p), "", p/3.0, 0.0, (p+1)/3.0, 0.35);
      top->SetLeftMargin(0.16); top->SetRightMargin(0.03); top->SetBottomMargin(0.02); top->SetTopMargin(0.05);
      bot->SetLeftMargin(0.16); bot->SetRightMargin(0.03); bot->SetTopMargin(0.02); bot->SetBottomMargin(0.3);
      top->Draw(); bot->Draw();
      top->cd();
      TH1D * s[kN]; double ymax = 0;
      for (int k = 0; k < kN; k++) { s[k] = shape(h[r][p][k], Form("s_%s_r%d_pt%d", cname[k], ir, p)); ymax = std::max(ymax, s[k]->GetMaximum()); }
      TLegend * l = new TLegend(0.6, 0.5, 0.97, 0.78); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.045);
      bool first = true;
      for (int k : {kSA, kBA, kBC, kCall}) {
        s[k]->SetLineColor(col[k]); s[k]->SetMarkerColor(col[k]); s[k]->SetMarkerStyle(mk[k]); s[k]->SetMarkerSize(0.8);
        s[k]->GetXaxis()->SetRangeUser(0, 2); s[k]->GetXaxis()->SetLabelSize(0);
        s[k]->GetYaxis()->SetTitleSize(0.055); s[k]->GetYaxis()->SetLabelSize(0.045); s[k]->GetYaxis()->SetTitleOffset(1.3);
        s[k]->SetMinimum(0); s[k]->SetMaximum(1.7*ymax);
        s[k]->Draw(first ? "e" : "e same"); first = false;
        l->AddEntry(s[k], shortLab[k], "lp"); // means: summary page and log
      }
      l->Draw();
      d.drawAll({"Pythia8 QCD (Jet12_long) MC"}, {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ana::ptBinsUsed[p], ana::ptBinsUsed[p+1]),
                Form("Jet R=%.1f", ana::JetRs[ir])}, .2, .88, 13, gPad->GetWh()*0.8);
      bot->cd();
      TH1D * rBC = (TH1D*)s[kBC]->Clone(Form("rBC_%d_%d", r, p)); rBC->Divide(s[kBA]);
      TH1D * rSA = (TH1D*)s[kSA]->Clone(Form("rSA_%d_%d", r, p)); rSA->Divide(s[kBA]);
      sanitize(rBC); sanitize(rSA);
      rBC->SetMinimum(0); rBC->SetMaximum(2.5);
      rBC->GetYaxis()->SetTitle("ratio to B_{A}"); rBC->GetYaxis()->SetNdivisions(505);
      rBC->GetYaxis()->SetTitleSize(0.1); rBC->GetYaxis()->SetLabelSize(0.08); rBC->GetYaxis()->SetTitleOffset(0.7);
      rBC->GetXaxis()->SetTitleSize(0.11); rBC->GetXaxis()->SetLabelSize(0.09); rBC->GetXaxis()->SetLabelSize(0.09);
      rBC->GetXaxis()->SetRangeUser(0, 2);
      rBC->Draw("e"); rSA->Draw("e same");
      TLine * one = new TLine(0, 1, 2, 1); one->SetLineStyle(2); one->Draw();
    }
    c->SaveAs(pdf.c_str());
  }
  // mean summary
  c->Clear(); c->Divide(3, 1);
  for (int r = 0; r < nR; r++) {
    c->cd(r+1); gPad->SetLeftMargin(0.16); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05);
    TH1D * fr = new TH1D(Form("frm%d", r), ";p_{T}^{#gamma} [GeV];#LTx_{J#gamma}#GT", 1, ana::ptBinsUsed[0], ana::ptBinsUsed[nPt]);
    fr->SetMinimum(0.3); fr->SetMaximum(1.3);
    fr->GetYaxis()->SetTitleSize(0.05); fr->GetXaxis()->SetTitleSize(0.05); fr->GetYaxis()->SetTitleOffset(1.4); fr->Draw();
    TLegend * l = new TLegend(0.2, 0.15, 0.7, 0.35); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.04);
    for (int k : {kSA, kBA, kSC, kBC}) {
      gMean[r][k]->SetMarkerStyle(mk[k]); gMean[r][k]->SetMarkerColor(col[k]); gMean[r][k]->SetLineColor(col[k]);
      gMean[r][k]->Draw("p same"); l->AddEntry(gMean[r][k], lab[k], "p");
    }
    l->Draw();
    d.drawAll({"Pythia8 QCD (Jet12_long) MC"}, {Form("Jet R=%.1f", ana::JetRs[irs[r]])}, .2, .88, 13, gPad->GetWh()*0.8);
  }
  c->SaveAs(pdf.c_str());
  c->SaveAs((pdf+"]").c_str());

  TFile * fout = TFile::Open((outdir + "/draw_jet12_sig_bkg_xj.root").c_str(), "recreate");
  for (int r = 0; r < nR; r++) for (int p = 0; p < nPt; p++) for (int k = 0; k < kN; k++) h[r][p][k]->Write();
  fout->Close();
}
