#include "../../src/ana.h"
#include "../../src/drawer.h"
#include "../../src/pho_object.h"
#include "../../src/jet_object.h"
#include "../../src/reweight_utility.h"
// Reuse the production purity solver itself (leakage-corrected ABCD quadratic + bootstrap).
#include "../../macros/puritymaker.C"
// Explicit load; run interpreted, never with ACLiC "+" (sibling libgammajet.so collision).
R__LOAD_LIBRARY(libgammajet_unfold.so);

// PPG18 review round 2, issue 10: photon purity as a function of |eta^gamma|.
//
// The production purity (macros/puritymaker.C) is one number per photon-pT bin and R,
// integrated over |eta| < 1.1. The reviewer notes region C is twice as concentrated at
// |eta| > 0.9 as region A (the isolation cone leaks out of the EMCal for |eta| > 0.7), so
// the background composition - and the ABCD independence assumption - may depend on eta.
//
// Here the SAME chain is run in |eta| bins:
//   data ABCD counts   paired clusters per photon-pT bin (= hclusterpt_abcd%i_%i), from
//                      trees/gammajet_Data.root with unfolder.cc's data branch:
//                      jet pT = jet_pt_calib/ana::jesNominal, ispaired incl. jet pT >
//                      jet_calib_pt_cut and check_pair (xJ floor, |eta|, dphi)
//   leakage templates  truth-matched (dR < 0.1) paired clusters from Photon5+10+20
//                      (= hclusterpt_abcd_truthmatched%i_%i), unfolder.cc's MC branch:
//                      EMR-smeared photon pT, jet_pt_smear_truth, treeuser.h truth-photon
//                      windows (0-12 / 12-24 / 24-100 GeV), vz/cluster-pT mcWeight,
//                      stitched with drawer's scalemap weights
//   purity             puritymaker.C::combine_hists (bootstrap median, 16/84% errors)
// The all-|eta| set is validated against the committed hists/*_nominal_unfolding.root
// hclusterpt_abcd counts and hists/purity_nominal.root.
//
// Key number: if purity varies with eta, the eta-integrated ABCD purity need not equal the
// region-A-weighted average of the per-eta purities (the solve is non-linear in the counts
// and the leakage fractions differ with eta). That difference is the bias of using one
// eta-integrated purity.

namespace {
  const int nEta = 3;
  const double etaEdges[nEta+1] = {0, 0.35, 0.7, 1.1};
  const int nSets = nEta + 1;                  // 0..nEta-1 = |eta| bins, nEta = all |eta|
  const int irs[2] = {0, 2};                   // R = 0.2, 0.4
  const int nR = 2;
  int etaBin(float aeta) { for (int i = 0; i < nEta; i++) if (aeta >= etaEdges[i] && aeta < etaEdges[i+1]) return i; return -1; }
  const char * setName(int s) { return s == nEta ? "all" : Form("eta%d", s); }
}

void purity_vs_eta()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();

  // [r][set][abcd]
  TH1D * hD[nR][nSets][4], * hM[nR][nSets][4];
  for (int r = 0; r < nR; r++) for (int s = 0; s < nSets; s++) for (int k = 0; k < 4; k++) {
    hD[r][s][k] = new TH1D(Form("data_abcd_r%d_%s_%d", irs[r], setName(s), k), ";p_{T}^{lead cluster};Counts", ana::nPtBins, ana::ptBins);
    hM[r][s][k] = new TH1D(Form("mc_abcd_tm_r%d_%s_%d", irs[r], setName(s), k), ";p_{T}^{lead cluster};Counts", ana::nPtBins, ana::ptBins);
  }

  auto process = [&](const char * fname, bool isMC, double scale, float truthLo, float truthHi) {
    TFile * f = TFile::Open(fname, "read");
    TTree * t = (TTree*)f->Get("towerntup");
    float vz, cluster_pt, cluster_e, cluster_eta, cluster_phi, cluster_time;
    float cluster_showershape[12], cluster_bdt_scores[11];
    float truth_cluster_pt = 0, truth_cluster_e = 0, truth_cluster_eta = 0, truth_cluster_phi = 0, truth_cluster_iso3 = 0, truth_cluster_iso4 = 0;
    float jet_pt_calib[7], jet_pt_smear_truth[7], jet_e[7], jet_eta[7], jet_phi[7], jet_emfrac[7];
    t->SetBranchStatus("*", 0);
    vector<const char*> br = {"vz","cluster_pt","cluster_e","cluster_eta","cluster_phi","cluster_time","cluster_showershape",
                              "cluster_bdt_scores","jet_e","jet_eta","jet_phi","jet_emfrac"};
    if (isMC) for (const char * b : {"truth_cluster_pt","truth_cluster_e","truth_cluster_eta","truth_cluster_phi",
                                     "truth_cluster_iso3","truth_cluster_iso4","jet_pt_smear_truth"}) br.push_back(b);
    else br.push_back("jet_pt_calib");
    for (auto b : br) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("vz", &vz);
    t->SetBranchAddress("cluster_pt", &cluster_pt);
    t->SetBranchAddress("cluster_e", &cluster_e);
    t->SetBranchAddress("cluster_eta", &cluster_eta);
    t->SetBranchAddress("cluster_phi", &cluster_phi);
    t->SetBranchAddress("cluster_time", &cluster_time);
    t->SetBranchAddress("cluster_showershape", cluster_showershape);
    t->SetBranchAddress("cluster_bdt_scores", cluster_bdt_scores);
    t->SetBranchAddress("jet_e", jet_e);
    t->SetBranchAddress("jet_eta", jet_eta);
    t->SetBranchAddress("jet_phi", jet_phi);
    t->SetBranchAddress("jet_emfrac", jet_emfrac);
    if (isMC) {
      t->SetBranchAddress("truth_cluster_pt", &truth_cluster_pt);
      t->SetBranchAddress("truth_cluster_e", &truth_cluster_e);
      t->SetBranchAddress("truth_cluster_eta", &truth_cluster_eta);
      t->SetBranchAddress("truth_cluster_phi", &truth_cluster_phi);
      t->SetBranchAddress("truth_cluster_iso3", &truth_cluster_iso3);
      t->SetBranchAddress("truth_cluster_iso4", &truth_cluster_iso4);
      t->SetBranchAddress("jet_pt_smear_truth", jet_pt_smear_truth);
    } else t->SetBranchAddress("jet_pt_calib", jet_pt_calib);

    Reweighter rw;
    TRandom3 rand(18);
    Long64_t n = t->GetEntries();
    for (Long64_t e = 0; e < n; e++) {
      t->GetEntry(e);
      if (fabs(vz) > ana::vzcut) continue;
      if (isMC && !(truth_cluster_pt > truthLo && truth_cluster_pt < truthHi)) continue; // treeuser.h photon window
      float w = isMC ? scale*rw.GetWeight(vz, cluster_pt) : 1.0;
      float recoPt = isMC ? cluster_pt + rand.Gaus(0, ana::emResolutionSigma(truth_cluster_pt, ana::emrNominal)*truth_cluster_pt) : cluster_pt;
      pho_object pho(recoPt, cluster_e, cluster_eta, cluster_phi, cluster_showershape[10], cluster_showershape[11],
          cluster_time, cluster_bdt_scores[9], pho_object::get_showershape(cluster_showershape, recoPt));
      if (!(pho.pt >= ana::ptBins[0] && pho.pt < ana::ptBins[ana::nPtBins])) continue;
      int iabcd = ana::findabcdBin(pho.iso4, pho.bdt, 0);
      if (iabcd < 0) continue;
      int ie = etaBin(fabs(pho.eta));
      bool matched = false;
      if (isMC) {
        pho_object tru(truth_cluster_pt, truth_cluster_e, truth_cluster_eta, truth_cluster_phi,
            truth_cluster_iso3, truth_cluster_iso4, 0, 0.99, 2);
        matched = pho.deltaR(tru) < 0.1; // unfolder::check_match
        if (!matched) continue;          // MC is only used for the truth-matched templates
      }
      int ptbin = ana::findPtBin(pho.pt);
      for (int r = 0; r < nR; r++) {
        int ir = irs[r];
        float jpt = isMC ? jet_pt_smear_truth[ir] : jet_pt_calib[ir]/ana::jesNominal[ir];
        if (!(jpt > ana::jet_calib_pt_cut[ir])) continue;
        jet_object jet(jpt, jet_e[ir], jet_eta[ir], jet_phi[ir], jet_emfrac[ir], 0, 0, 0);
        // unfolder::check_pair
        float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(ana::jet_calib_pt_cut[ir]/ana::ptBins[ptbin])+1];
        if (jet.pt/pho.pt < lowbin) continue;
        if (fabs(pho.eta) > ana::etacut) continue;
        if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) continue;
        if (jet.deltaPhi(pho) < ana::oppcut) continue;
        TH1D * (*h)[nSets][4] = isMC ? hM : hD;
        h[r][nEta][iabcd]->Fill(pho.pt, w);
        if (ie >= 0) h[r][ie][iabcd]->Fill(pho.pt, w);
      }
    }
    f->Close();
  };

  const string tdir = ana::path("trees/");
  process((tdir + "gammajet_Data.root").c_str(), false, 1, 0, 0);
  cout << "data done" << endl;
  // drawer.h scalemap[isphoton=1] for sim="pythia" (same numbers as grid_insitu.C), and
  // treeuser.h threshmap/threshmap_high[-1] truth-photon windows
  process((tdir + "gammajet_pythia_Photon5.root").c_str(),  true, 146359.3,  0, 12);
  process((tdir + "gammajet_pythia_Photon10.root").c_str(), true, 6944.675, 12, 24);
  process((tdir + "gammajet_pythia_Photon20.root").c_str(), true, 130.4461, 24, 100);
  cout << "MC done" << endl;

  // ---------------- validation of the all-|eta| set ----------------
  TFile * fData = TFile::Open(ana::path("hists/Data_nominal_unfolding.root"), "read");
  drawer d("pythia", "nominal");
  printf("\nValidation (all |eta|): counts here vs committed hclusterpt_abcd; leakage c here vs committed\n");
  for (int r = 0; r < nR; r++) {
    int ir = irs[r];
    for (int ib = ana::firstUsedPtBin; ib < ana::firstUsedPtBin + ana::nPtBinsUsed; ib++) {
      printf("R=%.1f %2.0f-%-2.0f data:", ana::JetRs[ir], ana::ptBins[ib], ana::ptBins[ib+1]);
      for (int k = 0; k < 4; k++) {
        TH1D * hc = (TH1D*)fData->Get(Form("hclusterpt_abcd%d_%d", ir, k));
        printf(" %c %5.0f(%5.0f)", 'A'+k, hD[r][nEta][k]->GetBinContent(ib+1), hc ? hc->GetBinContent(ib+1) : -1.);
      }
      TH1D * c2 = d.get(Form("hclusterpt_abcd_truthmatched%d_2", ir), 1), * c0 = d.get(Form("hclusterpt_abcd_truthmatched%d_0", ir), 1);
      printf(" | c %.3f(%.3f)\n", hM[r][nEta][2]->GetBinContent(ib+1)/hM[r][nEta][0]->GetBinContent(ib+1), c2->GetBinContent(ib+1)/c0->GetBinContent(ib+1));
    }
  }

  // ---------------- purity per set ----------------
  string outdir = ana::path("claude_checks/purity/pdfs");
  gSystem->mkdir(outdir.c_str(), true);
  TFile * fout = TFile::Open((outdir + "/purity_vs_eta.root").c_str(), "recreate");
  const int nPt = ana::nPtBinsUsed;
  double P[nR][nSets][nPt], Plo[nR][nSets][nPt], Phi[nR][nSets][nPt], PC[nR][nSets][nPt];
  double fb[nR][nSets][nPt], fc[nR][nSets][nPt], fd[nR][nSets][nPt];
  for (int r = 0; r < nR; r++) for (int s = 0; s < nSets; s++) {
    TH1D * fp[4];
    for (int k = 0; k < 4; k++) { fp[k] = (TH1D*)hM[r][s][k]->Clone(Form("fp_r%d_%s_%d", irs[r], setName(s), k)); fp[k]->Divide(hM[r][s][k], hM[r][s][0]); }
    fout->mkdir(Form("r%d_%s", irs[r], setName(s)))->cd();
    for (int k = 0; k < 4; k++) { hD[r][s][k]->Write(); hM[r][s][k]->Write(); }
    TGraphAsymmErrors * odC = nullptr;
    TGraphAsymmErrors * od = combine_hists(hD[r][s], fp, &odC);
    fout->cd();
    for (int p = 0; p < nPt; p++) {
      int ib = ana::firstUsedPtBin + p;
      double x;
      od->GetPoint(ib, x, P[r][s][p]); Plo[r][s][p] = od->GetErrorYlow(ib); Phi[r][s][p] = od->GetErrorYhigh(ib);
      if (odC) odC->GetPoint(ib, x, PC[r][s][p]); else PC[r][s][p] = 0;
      fb[r][s][p] = fp[1]->GetBinContent(ib+1); fc[r][s][p] = fp[2]->GetBinContent(ib+1); fd[r][s][p] = fp[3]->GetBinContent(ib+1);
    }
  }

  // ---------------- numbers ----------------
  printf("\nPurity vs |eta| (bootstrap median -lo/+hi); leakage b/c/d from Photon MC in the same |eta| bin\n");
  for (int r = 0; r < nR; r++) {
    int ir = irs[r];
    printf("\n===== R = %.1f =====  (committed purity_nominal.root, all |eta|: %.3f %.3f %.3f)\n", ana::JetRs[ir],
           ana::getPurity(15, 20, "nominal", ir), ana::getPurity(20, 25, "nominal", ir), ana::getPurity(25, 35, "nominal", ir));
    for (int p = 0; p < nPt; p++) {
      int ib = ana::firstUsedPtBin + p;
      printf("-- %2.0f-%2.0f GeV --\n", ana::ptBinsUsed[p], ana::ptBinsUsed[p+1]);
      double sumA = 0, sumPA = 0;
      for (int s = 0; s < nSets; s++) {
        double A = hD[r][s][0]->GetBinContent(ib+1), B = hD[r][s][1]->GetBinContent(ib+1), C = hD[r][s][2]->GetBinContent(ib+1), D = hD[r][s][3]->GetBinContent(ib+1);
        printf("  %-15s A %5.0f B %4.0f C %5.0f D %5.0f | C/A %.2f B/D %.2f | b %.3f c %.3f d %.3f | P_A %.3f -%.3f +%.3f  P_C %.3f\n",
               s == nEta ? "all |eta|" : Form("%.2f-%.2f", etaEdges[s], etaEdges[s+1]), A, B, C, D, A > 0 ? C/A : 0, D > 0 ? B/D : 0,
               fb[r][s][p], fc[r][s][p], fd[r][s][p], P[r][s][p], Plo[r][s][p], Phi[r][s][p], PC[r][s][p]);
        if (s < nEta) { sumA += A; sumPA += A*P[r][s][p]; }
      }
      printf("  A-weighted average of per-|eta| P_A = %.3f   vs  |eta|-integrated P_A = %.3f   (difference %+.3f)\n",
             sumA > 0 ? sumPA/sumA : 0, P[r][nEta][p], (sumA > 0 ? sumPA/sumA : 0) - P[r][nEta][p]);
    }
  }

  // ---------------- drawing ----------------
  string pdf = outdir + "/purity_vs_eta.pdf";
  TCanvas * c = new TCanvas("c", "", 1400, 600);
  c->SaveAs((pdf+"[").c_str());
  int cols[nSets] = {kBlue+1, kGreen+2, kRed+1, kBlack};
  int mks[nSets] = {21, 22, 23, 20};
  // page 1: P_A vs pT, one curve per |eta| bin, one pad per R
  c->Clear(); c->Divide(2, 1);
  for (int r = 0; r < nR; r++) {
    c->cd(r+1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05);
    TH1D * fr = new TH1D(Form("frPe%d", r), ";p_{T}^{#gamma} [GeV];purity P_{A}", 1, ana::ptBinsUsed[0], ana::ptBinsUsed[nPt]);
    fr->SetMinimum(0); fr->SetMaximum(1.5); fr->GetYaxis()->SetTitleSize(0.05); fr->GetXaxis()->SetTitleSize(0.05); fr->Draw();
    TLegend * l = new TLegend(0.55, 0.62, 0.97, 0.88); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.035);
    for (int s = 0; s < nSets; s++) {
      TGraphAsymmErrors * g = new TGraphAsymmErrors();
      for (int p = 0; p < nPt; p++) {
        g->SetPoint(p, 0.5*(ana::ptBinsUsed[p]+ana::ptBinsUsed[p+1]) + (s-1.5)*0.5, P[r][s][p]);
        g->SetPointError(p, 0, 0, Plo[r][s][p], Phi[r][s][p]);
      }
      g->SetMarkerStyle(mks[s]); g->SetMarkerColor(cols[s]); g->SetLineColor(cols[s]); g->SetMarkerSize(1.2);
      g->Draw("p same");
      l->AddEntry(g, s == nEta ? "all |#eta^{#gamma}| < 1.1" : Form("%.2f < |#eta^{#gamma}| < %.2f", etaEdges[s], etaEdges[s+1]), "p");
    }
    l->Draw();
    d.drawAll({"p+p Run24 Data"}, {Form("Jet R=%.1f, paired", ana::JetRs[irs[r]]), "leakage-corrected ABCD (puritymaker.C)", "Photon MC leakage per |#eta| bin"}, .17, .9, 13, gPad->GetWh()*0.8);
  }
  c->SaveAs(pdf.c_str());
  // page 2: C/A and leakage c vs |eta| (15-20 GeV), R = 0.4
  c->Clear(); c->Divide(2, 1);
  {
    int r = 1, p = 0, ib = ana::firstUsedPtBin;
    TGraph * gCA = new TGraph(), * gBD = new TGraph(), * gc = new TGraph(), * gb = new TGraph();
    for (int s = 0; s < nEta; s++) {
      double x = 0.5*(etaEdges[s]+etaEdges[s+1]);
      double A = hD[r][s][0]->GetBinContent(ib+1), B = hD[r][s][1]->GetBinContent(ib+1), C = hD[r][s][2]->GetBinContent(ib+1), D = hD[r][s][3]->GetBinContent(ib+1);
      gCA->SetPoint(s, x, C/A); gBD->SetPoint(s, x, B/D); gc->SetPoint(s, x, fc[r][s][p]); gb->SetPoint(s, x, 10*fb[r][s][p]);
    }
    c->cd(1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.05); gPad->SetRightMargin(0.03);
    TH1D * fr1 = new TH1D("frCA", ";|#eta^{#gamma}|;ratio", 1, 0, 1.1); fr1->SetMinimum(0); fr1->SetMaximum(2.0); fr1->Draw();
    gCA->SetMarkerStyle(20); gBD->SetMarkerStyle(24); gCA->Draw("pl same"); gBD->Draw("pl same");
    TLegend * l1 = new TLegend(0.55, 0.65, 0.95, 0.8); l1->SetBorderSize(0); l1->SetFillStyle(0); l1->SetTextSize(0.04);
    l1->AddEntry(gCA, "data C/A", "p"); l1->AddEntry(gBD, "data B/D", "p"); l1->Draw();
    d.drawAll({"p+p Run24 Data"}, {"15 GeV < p_{T}^{#gamma} < 20 GeV", "Jet R=0.4, paired"}, .17, .88, 13, gPad->GetWh()*0.8);
    c->cd(2); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetTopMargin(0.05); gPad->SetRightMargin(0.03);
    TH1D * fr2 = new TH1D("frc", ";|#eta^{#gamma}|;leakage fraction", 1, 0, 1.1); fr2->SetMinimum(0); fr2->SetMaximum(0.6); fr2->Draw();
    gc->SetMarkerStyle(21); gb->SetMarkerStyle(25); gc->Draw("pl same"); gb->Draw("pl same");
    TLegend * l2 = new TLegend(0.55, 0.65, 0.95, 0.8); l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextSize(0.04);
    l2->AddEntry(gc, "c = N_{sig}^{C}/N_{sig}^{A}", "p"); l2->AddEntry(gb, "b #times 10", "p"); l2->Draw();
    d.drawAll({"Pythia8 #gamma+jet MC"}, {"15 GeV < p_{T}^{#gamma} < 20 GeV", "Jet R=0.4, truth-matched"}, .17, .88, 13, gPad->GetWh()*0.8);
  }
  c->SaveAs(pdf.c_str());
  c->SaveAs((pdf+"]").c_str());
  fout->Close();
}
