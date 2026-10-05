#include "../../src/ana.h"
#include "../../src/drawer.h"
#include "../../src/pho_object.h"
#include "../../src/jet_object.h"
#include "../../src/unfold_utility.h"
// Reuse the production purity solver itself (leakage-corrected ABCD quadratic + bootstrap).
#include "../../macros/puritymaker.C"
// Explicit load; run interpreted, never with ACLiC "+" (sibling libgammajet.so collision).
R__LOAD_LIBRARY(libgammajet_unfold.so);

// PPG18 review round 2, issue 10: R = 0.4 data xJ with the nominal photon acceptance
// (|eta^gamma| < 1.1) vs a PPG12-style |eta^gamma| < 0.7 cut.
//
// Reco level (the unfolding input), not unfolded - unfolding with the cut needs the MC
// response rebuilt with the same cut (pipeline rerun).
//   data        trees/gammajet_Data.root, unfolder.cc data branch: jet pT =
//               jet_pt_calib/ana::jesNominal, check_pair (xJ floor, jet |eta| < 1.1-R,
//               dphi > 7pi/8), regions A and C from ana::findabcdBin(iso4, bdt, 0)
//   purity      |eta| < 1.1: committed hists/purity_nominal.root (ana::getPurity/C)
//               |eta| < 0.7: puritymaker.C::combine_hists on the data ABCD counts and
//               Photon5+10+20 truth-matched leakage templates restricted to |eta| < 0.7,
//               summed from purity_vs_eta.C's 0-0.35 and 0.35-0.7 bins (run that first)
//   correction  unfold_utility::purityCorrect(A, C, P_A, P_C) - the production call
// Stat and purity uncertainties are drawn as the corrected histograms' own bin errors
// (as purityCorrect stores them); no systematic band is drawn.

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
  double binMean(TH1D * h, double * err = nullptr) {
    double s = 0, sx = 0, ve = 0;
    for (int b = 1; b <= h->GetNbinsX(); b++) { s += h->GetBinContent(b); sx += h->GetBinContent(b)*h->GetBinCenter(b); }
    if (s <= 0) { if (err) *err = 0; return 0; }
    double m = sx/s;
    for (int b = 1; b <= h->GetNbinsX(); b++) ve += pow((h->GetBinCenter(b)-m)*h->GetBinError(b)/s, 2);
    if (err) *err = sqrt(ve);
    return m;
  }
}

void draw_xj_eta07()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  const int ir = 2;           // R = 0.4
  const int nPt = ana::nPtBinsUsed;
  const int nSel = 2;         // 0: |eta| < 1.1 (nominal), 1: |eta| < 0.7
  const double etaMax[nSel] = {ana::etacut, 0.7};
  const char * selLab[nSel] = {"|#eta^{#gamma}| < 1.1", "|#eta^{#gamma}| < 0.7"};

  TH1D * hA[nSel][nPt], * hC[nSel][nPt];
  for (int s = 0; s < nSel; s++) for (int p = 0; p < nPt; p++) {
    hA[s][p] = new TH1D(Form("hA_s%d_pt%d", s, p), ";x_{J#gamma};(1/N) dN/dx_{J#gamma}", ana::nUnfoldXjBins, ana::unfoldXjBins);
    hC[s][p] = new TH1D(Form("hC_s%d_pt%d", s, p), ";x_{J#gamma};(1/N) dN/dx_{J#gamma}", ana::nUnfoldXjBins, ana::unfoldXjBins);
  }

  TFile * fin = TFile::Open(ana::path("trees/gammajet_Data.root"), "read");
  TTree * t = (TTree*)fin->Get("towerntup");
  float vz, cluster_pt, cluster_e, cluster_eta, cluster_phi, cluster_time;
  float cluster_showershape[12], cluster_bdt_scores[11];
  float jet_pt_calib[7], jet_e[7], jet_eta[7], jet_phi[7], jet_emfrac[7];
  t->SetBranchStatus("*", 0);
  for (const char * b : {"vz","cluster_pt","cluster_e","cluster_eta","cluster_phi","cluster_time",
        "cluster_showershape","cluster_bdt_scores","jet_pt_calib","jet_e","jet_eta","jet_phi","jet_emfrac"})
    t->SetBranchStatus(b, 1);
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
  for (Long64_t e = 0; e < t->GetEntries(); e++) {
    t->GetEntry(e);
    if (fabs(vz) > ana::vzcut) continue;
    pho_object pho(cluster_pt, cluster_e, cluster_eta, cluster_phi, cluster_showershape[10], cluster_showershape[11],
        cluster_time, cluster_bdt_scores[9], pho_object::get_showershape(cluster_showershape, cluster_pt));
    int iabcd = ana::findabcdBin(pho.iso4, pho.bdt, 0);
    if (iabcd != 0 && iabcd != 2) continue;
    int ptbin = ana::findPtBin(pho.pt);
    int p = ptbin - ana::firstUsedPtBin;
    if (p < 0 || p >= nPt) continue;
    float jpt = jet_pt_calib[ir]/ana::jesNominal[ir];
    if (!(jpt > ana::jet_calib_pt_cut[ir])) continue;
    jet_object jet(jpt, jet_e[ir], jet_eta[ir], jet_phi[ir], jet_emfrac[ir], 0, 0, 0);
    float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(ana::jet_calib_pt_cut[ir]/ana::ptBins[ptbin])+1];
    float xj = jet.pt/pho.pt;
    if (xj < lowbin) continue;
    if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) continue;
    if (jet.deltaPhi(pho) < ana::oppcut) continue;
    for (int s = 0; s < nSel; s++) {
      if (fabs(pho.eta) > etaMax[s]) continue;
      (iabcd == 0 ? hA : hC)[s][p]->Fill(xj);
    }
  }

  // ---------------- purities ----------------
  float pA[nSel][nPt], pAlo[nSel][nPt], pAhi[nSel][nPt], pC[nSel][nPt], pClo[nSel][nPt], pChi[nSel][nPt];
  for (int p = 0; p < nPt; p++) {
    float lo = ana::ptBinsUsed[p], hi = ana::ptBinsUsed[p+1];
    pA[0][p] = ana::getPurity(lo, hi, "nominal", ir);  pAlo[0][p] = ana::getPurityErrorLow(lo, hi, "nominal", ir);  pAhi[0][p] = ana::getPurityErrorHigh(lo, hi, "nominal", ir);
    pC[0][p] = ana::getPurityC(lo, hi, "nominal", ir); pClo[0][p] = ana::getPurityCErrorLow(lo, hi, "nominal", ir); pChi[0][p] = ana::getPurityCErrorHigh(lo, hi, "nominal", ir);
  }
  string outdir = ana::path("claude_checks/purity/pdfs");
  TFile * fe = TFile::Open((outdir + "/purity_vs_eta.root").c_str(), "read");
  if (!fe || fe->IsZombie()) { cout << "Run purity_vs_eta.C first." << endl; return; }
  TH1D * hd[4], * hm[4], * fp[4];
  for (int k = 0; k < 4; k++) {
    hd[k] = (TH1D*)fe->Get(Form("r%d_eta0/data_abcd_r%d_eta0_%d", ir, ir, k))->Clone(Form("d07_%d", k));
    hd[k]->Add((TH1D*)fe->Get(Form("r%d_eta1/data_abcd_r%d_eta1_%d", ir, ir, k)));
    hm[k] = (TH1D*)fe->Get(Form("r%d_eta0/mc_abcd_tm_r%d_eta0_%d", ir, ir, k))->Clone(Form("m07_%d", k));
    hm[k]->Add((TH1D*)fe->Get(Form("r%d_eta1/mc_abcd_tm_r%d_eta1_%d", ir, ir, k)));
    hd[k]->SetDirectory(0); hm[k]->SetDirectory(0);
  }
  for (int k = 0; k < 4; k++) { fp[k] = (TH1D*)hm[k]->Clone(Form("fp07_%d", k)); fp[k]->Divide(hm[k], hm[0]); }
  TFile * fout = TFile::Open((outdir + "/draw_xj_eta07.root").c_str(), "recreate");
  fout->mkdir("purity_eta07")->cd(); // combine_hists writes into the current directory
  TGraphAsymmErrors * odC = nullptr;
  TGraphAsymmErrors * od = combine_hists(hd, fp, &odC);
  fout->cd();
  for (int p = 0; p < nPt; p++) {
    int ib = ana::firstUsedPtBin + p;
    double x, y;
    od->GetPoint(ib, x, y);  pA[1][p] = y; pAlo[1][p] = od->GetErrorYlow(ib);  pAhi[1][p] = od->GetErrorYhigh(ib);
    odC->GetPoint(ib, x, y); pC[1][p] = y; pClo[1][p] = odC->GetErrorYlow(ib); pChi[1][p] = odC->GetErrorYhigh(ib);
    // sanity: the |eta|<0.7 ABCD counts summed from purity_vs_eta.C must equal the A/C filled here
    printf("check %2.0f-%2.0f: A %0.f vs %0.f, C %0.f vs %0.f\n", ana::ptBinsUsed[p], ana::ptBinsUsed[p+1],
           hd[0]->GetBinContent(ib+1), hA[1][p]->Integral(), hd[2]->GetBinContent(ib+1), hC[1][p]->Integral());
  }

  // ---------------- correct + numbers ----------------
  TH1D * hCorr[nSel][nPt];
  printf("\nR = 0.4, reco level. Binned <xJ> (errors x1000).\n");
  printf("%-9s | %-26s | %-9s %-9s %-9s | %-9s %-9s\n", "pT", "selection", "N_A", "N_C", "P_A (P_C)", "<A>", "<corr>");
  for (int p = 0; p < nPt; p++) for (int s = 0; s < nSel; s++) {
    hCorr[s][p] = unfold_utility::purityCorrect(hA[s][p], hC[s][p], pA[s][p], pAlo[s][p], pAhi[s][p],
        pC[s][p], pClo[s][p], pChi[s][p], Form("hCorr_s%d_pt%d", s, p), nullptr, true);
    double e1, e2, mA = binMean(hA[s][p], &e1), mC = binMean(hCorr[s][p], &e2);
    printf("%3.0f-%-5.0f | %-26s | %-9.0f %-9.0f %.3f (%.3f) | %.3f(%2.0f) %.3f(%2.0f)\n", ana::ptBinsUsed[p], ana::ptBinsUsed[p+1],
           s ? "|eta| < 0.7" : "|eta| < 1.1 (nominal)", hA[s][p]->Integral(), hC[s][p]->Integral(), pA[s][p], pC[s][p], mA, 1000*e1, mC, 1000*e2);
  }

  // ---------------- drawing ----------------
  drawer d("pythia", "nominal"); // ctor args only pick files for other helpers; labels come from drawAll
  string pdf = outdir + "/draw_xj_eta07.pdf";
  TCanvas * c = new TCanvas("c", "", 1500, 800);
  c->SaveAs((pdf+"[").c_str());
  for (int page = 0; page < 2; page++) { // 0: purity corrected, 1: uncorrected region A
    c->Clear();
    for (int p = 0; p < nPt; p++) {
      c->cd();
      TPad * top = new TPad(Form("etop%d_%d", page, p), "", p/3.0, 0.35, (p+1)/3.0, 1.0);
      TPad * bot = new TPad(Form("ebot%d_%d", page, p), "", p/3.0, 0.0, (p+1)/3.0, 0.35);
      top->SetLeftMargin(0.16); top->SetRightMargin(0.03); top->SetBottomMargin(0.02); top->SetTopMargin(0.05);
      bot->SetLeftMargin(0.16); bot->SetRightMargin(0.03); bot->SetTopMargin(0.02); bot->SetBottomMargin(0.3);
      top->Draw(); bot->Draw();
      top->cd();
      TH1D * s0 = shape(page == 0 ? hCorr[0][p] : hA[0][p], Form("s0_%d_%d", page, p));
      TH1D * s1 = shape(page == 0 ? hCorr[1][p] : hA[1][p], Form("s1_%d_%d", page, p));
      double ymax = std::max(s0->GetMaximum(), s1->GetMaximum());
      s0->SetLineColor(kBlack); s0->SetMarkerColor(kBlack); s0->SetMarkerStyle(20);
      s1->SetLineColor(kRed+1); s1->SetMarkerColor(kRed+1); s1->SetMarkerStyle(24);
      s0->GetXaxis()->SetRangeUser(0, 2); s0->GetXaxis()->SetLabelSize(0);
      s0->GetYaxis()->SetTitleSize(0.055); s0->GetYaxis()->SetLabelSize(0.045); s0->GetYaxis()->SetTitleOffset(1.3);
      s0->SetMinimum(0); s0->SetMaximum(1.7*ymax);
      s0->Draw("e"); s1->Draw("e same");
      TLegend * l = new TLegend(0.6, 0.6, 0.97, 0.76); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.045);
      l->AddEntry(s0, selLab[0], "lp");
      l->AddEntry(s1, selLab[1], "lp");
      l->Draw();
      d.drawAll({"p+p Run24 Data"}, {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ana::ptBinsUsed[p], ana::ptBinsUsed[p+1]),
                Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir]),
                page == 0 ? "purity corrected, reco level" : "region A, uncorrected, reco level",
                Form("P_{A} = %.2f (|#eta|<1.1), %.2f (|#eta|<0.7)", pA[0][p], pA[1][p])}, .2, .88, 13, gPad->GetWh()*0.8);
      bot->cd();
      TH1D * r = (TH1D*)s1->Clone(Form("r_%d_%d", page, p)); r->Divide(s0); sanitize(r);
      // |eta|<0.7 is a subset of |eta|<1.1: the ratio's bin errors (naive division) overestimate
      // its statistical scatter - read it as the size of the shape change only
      r->SetMinimum(0.5); r->SetMaximum(1.5);
      r->GetYaxis()->SetTitle("|#eta|<0.7 / nominal"); r->GetYaxis()->SetNdivisions(505);
      r->GetYaxis()->SetTitleSize(0.09); r->GetYaxis()->SetLabelSize(0.08); r->GetYaxis()->SetTitleOffset(0.8);
      r->GetXaxis()->SetTitleSize(0.11); r->GetXaxis()->SetLabelSize(0.09);
      r->GetXaxis()->SetRangeUser(0, 2);
      r->Draw("e");
      TLine * one = new TLine(0, 1, 2, 1); one->SetLineStyle(2); one->Draw();
    }
    c->SaveAs(pdf.c_str());
  }
  c->SaveAs((pdf+"]").c_str());
  fout->cd();
  for (int s = 0; s < nSel; s++) for (int p = 0; p < nPt; p++) { hA[s][p]->Write(); hC[s][p]->Write(); if (hCorr[s][p]) hCorr[s][p]->Write(); }
  fout->Close();
}
