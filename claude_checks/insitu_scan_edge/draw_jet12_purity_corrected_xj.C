#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/pho_object.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/jet_object.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/reweight_utility.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/unfold_utility.h"
// Reuse the production purity solver (leakage-corrected ABCD quadratic + bootstrap) itself
// rather than a copy of it - combine_hists() is the function macros/puritymaker.C runs.
#include "/home/samson72/sphnx/gammajet_unfold/macros/puritymaker.C"
// Explicit load; run interpreted, never with ACLiC "+" (sibling libgammajet.so collision).
R__LOAD_LIBRARY(libgammajet_unfold.so);

// PPG18 review round 2, issue 3 companion to draw_jet12_sig_bkg_xj.C: a closure test of
// the full purity-correction chain on the Pythia8 QCD "Jet12_long" sample, R = 0.2/0.3/0.4.
//
//   uncorrected  region-A paired xJ (what the Jet12 in-situ reference / Fig. 10 uses)
//   corrected    unfold_utility::purityCorrect(A, C, P_A, P_C) - the same call the
//                unfolding input and draw_insitu_xj.C use - with P_A, P_C derived FROM
//                THIS SAMPLE by the production method (puritymaker.C::combine_hists):
//                  - ABCD counts: Jet12_long paired clusters per photon-pT bin
//                    (= hclusterpt_abcd%i_%i, filled here from the tree)
//                  - leakage fractions f^X = N_sig^X/N_sig^A: Photon5+10+20 truth-matched
//                    templates (hclusterpt_abcd_truthmatched%i_%i from the committed
//                    hists/Photon*_pythia_nominal_unfolding.root, stitched by drawer::get
//                    exactly as puritymaker.C does)
//                  - central value = bootstrap median, errors = 16/84% quantiles
//   truth-tagged region-A reco photons within dR < 0.1 of the tree's truth photon - the
//                purity method's own signal definition (unfolder.cc hclusterpt_abcd_
//                truthmatched fill). Still RECO-level xJ (reco jet and photon pT): the
//                true-photon subset of region A, i.e. what the subtraction should leave.
//
// If the method closes, corrected ~ truth-tagged. The true Jet12 purity (matched/all in A) is
// printed next to the ABCD value. Selection/kinematics are identical to
// draw_jet12_sig_bkg_xj.C (unfolder.cc nominal MC branch). The committed
// Jet12_long_pythia_nominal_unfolding.root hclusterpt_abcd counts are printed next to the
// ones filled here as a cross-check of the selection (small differences are expected from
// the photon-smearing random sequence).

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
  // mean of a histogram from its bin contents (purityCorrect's output is a bin-by-bin
  // subtraction, so TH1's fill-time stats don't describe it)
  double binMean(TH1D * h, double * err = nullptr) {
    double s = 0, sx = 0, sx2 = 0, ve = 0;
    for (int b = 1; b <= h->GetNbinsX(); b++) { double c = h->GetBinContent(b), x = h->GetBinCenter(b); s += c; sx += c*x; sx2 += c*x*x; }
    if (s <= 0) return 0;
    double m = sx/s;
    for (int b = 1; b <= h->GetNbinsX(); b++) ve += pow((h->GetBinCenter(b)-m)*h->GetBinError(b)/s, 2);
    if (err) *err = sqrt(ve);
    return m;
  }
}

void draw_jet12_purity_corrected_xj()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();

  const int irs[3] = {0, 1, 2};
  const int nR = 3, nPt = ana::nPtBinsUsed;
  // ABCD cluster-pT counts per R (ana::ptBins binning, as hclusterpt_abcd)
  TH1D * hABCD[nR][4];
  // region A / C xJ per reported pT bin, and the truth references in A
  TH1D * hA[nR][nPt], * hC[nR][nPt], * hSig[nR][nPt];
  for (int r = 0; r < nR; r++) {
    for (int j = 0; j < 4; j++)
      hABCD[r][j] = new TH1D(Form("jet12_clusterpt_abcd%d_%d", irs[r], j), ";p_{T}^{lead cluster};Counts", ana::nPtBins, ana::ptBins);
    for (int p = 0; p < nPt; p++) {
      auto mk = [&](const char * tag) { return new TH1D(Form("h%s_r%d_pt%d", tag, irs[r], p), ";x_{J#gamma};(1/N) dN/dx_{J#gamma}", ana::nUnfoldXjBins, ana::unfoldXjBins); };
      hA[r][p] = mk("A"); hC[r][p] = mk("C"); hSig[r][p] = mk("Sig");
    }
  }

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
    if (iabcd < 0) continue;
    int ptbin = ana::findPtBin(pho.pt);
    if (ptbin < 0) continue;
    int ipt = ptbin - ana::firstUsedPtBin; // reported-bin index, may be out of [0,nPt)

    bool matched = false;
    if (truth_cluster_pt > 0) {
      pho_object tru(truth_cluster_pt, truth_cluster_e, truth_cluster_eta, truth_cluster_phi,
          truth_cluster_iso3, truth_cluster_iso4, 0, 0.99, 2);
      matched = pho.deltaR(tru) < 0.1;
    }

    for (int r = 0; r < nR; r++) {
      int ir = irs[r];
      if (!(truth_jet_pt[ir] > 0 && truth_jet_pt[ir] < 100)) continue; // Jet12_long keep window
      if (jet_pt_smear_truth[ir] <= 0) continue;
      jet_object jet(jet_pt_smear_truth[ir], jet_e[ir], jet_eta[ir], jet_phi[ir], jet_emfrac[ir], 0, 0, 0);
      // unfolder::check_pair
      float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(ana::jet_calib_pt_cut[ir]/ana::ptBins[ptbin])+1];
      float xj = jet.pt/pho.pt;
      if (xj < lowbin) continue;
      if (fabs(pho.eta) > ana::etacut) continue;
      if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) continue;
      if (jet.deltaPhi(pho) < ana::oppcut) continue;

      hABCD[r][iabcd]->Fill(pho.pt, w);
      if (ipt < 0 || ipt >= nPt) continue;
      if (iabcd == 0) {
        hA[r][ipt]->Fill(xj, w);
        if (matched) hSig[r][ipt]->Fill(xj, w);
      }
      if (iabcd == 2) hC[r][ipt]->Fill(xj, w);
    }
  }

  // ---------------- purity: production method ----------------
  drawer d("pythia", "nominal");
  TFile * fcommitted = TFile::Open("/home/samson72/sphnx/gammajet_unfold/hists/Jet12_long_pythia_nominal_unfolding.root", "read");
  string outdir = "/home/samson72/sphnx/gammajet_unfold/claude_checks/insitu_scan_edge/pdfs";
  gSystem->mkdir(outdir.c_str(), true);
  TFile * fout = TFile::Open((outdir + "/draw_jet12_purity_corrected_xj.root").c_str(), "recreate");
  float pA[nR][nPt], pAlo[nR][nPt], pAhi[nR][nPt], pC[nR][nPt], pClo[nR][nPt], pChi[nR][nPt];
  for (int r = 0; r < nR; r++) {
    int ir = irs[r];
    TH1D * fp[4], * hp[4];
    for (int j = 0; j < 4; j++) {
      hp[j] = d.get(Form("hclusterpt_abcd_truthmatched%i_%i", ir, j), 1);
      fp[j] = (TH1D*)hp[j]->Clone(Form("fp_r%d_%d", ir, j));
      fp[j]->Divide(hp[j], hp[0]);
    }
    fout->mkdir(ana::rnames[ir])->cd(); // combine_hists writes into the current directory
    TGraphAsymmErrors * odC = nullptr;
    cout << "\n--- combine_hists bootstrap quantiles, R = " << ana::JetRs[ir] << " ---" << endl;
    TGraphAsymmErrors * od = combine_hists(hABCD[r], fp, &odC);
    fout->cd();
    printf("R = %.1f  ABCD counts here vs committed hclusterpt_abcd (Jet12_long_pythia_nominal):\n", ana::JetRs[ir]);
    for (int p = 0; p < nPt; p++) {
      int ib = ana::firstUsedPtBin + p;
      printf("  %3.0f-%-3.0f GeV:", ana::ptBins[ib], ana::ptBins[ib+1]);
      for (int j = 0; j < 4; j++) {
        TH1D * hc = fcommitted ? (TH1D*)fcommitted->Get(Form("hclusterpt_abcd%i_%i", ir, j)) : nullptr;
        printf("  %c %7.1f (%7.1f)", 'A'+j, hABCD[r][j]->GetBinContent(ib+1), hc ? hc->GetBinContent(ib+1) : -1.0);
      }
      printf("   leakage b,c,d = %.3f %.3f %.3f\n", fp[1]->GetBinContent(ib+1), fp[2]->GetBinContent(ib+1), fp[3]->GetBinContent(ib+1));
      double x, y;
      od->GetPoint(ib, x, y); pA[r][p] = y; pAlo[r][p] = od->GetErrorYlow(ib); pAhi[r][p] = od->GetErrorYhigh(ib);
      odC->GetPoint(ib, x, y); pC[r][p] = y; pClo[r][p] = odC->GetErrorYlow(ib); pChi[r][p] = odC->GetErrorYhigh(ib);
    }
  }

  // ---------------- corrected distributions + numbers ----------------
  TH1D * hCorr[nR][nPt];
  printf("\n%-5s %-9s | %6s %6s %6s | %8s %8s %8s | %s\n", "R", "pT bin", "P_A", "P_A tr", "P_C",
         "<A>", "<corr>", "<tagged>", "corr/tagged - 1 (mean)");
  for (int r = 0; r < nR; r++) {
    int ir = irs[r];
    for (int p = 0; p < nPt; p++) {
      hCorr[r][p] = unfold_utility::purityCorrect(hA[r][p], hC[r][p], pA[r][p], pAlo[r][p], pAhi[r][p],
          pC[r][p], pClo[r][p], pChi[r][p], Form("hCorr_r%d_pt%d", ir, p), nullptr, true);
      if (!hCorr[r][p]) { hCorr[r][p] = (TH1D*)hA[r][p]->Clone(Form("hCorr_r%d_pt%d", ir, p)); hCorr[r][p]->Reset(); }
      double eA, eCo, eS;
      double mA = binMean(hA[r][p], &eA), mCo = binMean(hCorr[r][p], &eCo), mS = binMean(hSig[r][p], &eS);
      double pTrue = hSig[r][p]->Integral()/std::max(1e-9, hA[r][p]->Integral());
      printf("%-5.1f %3.0f-%-5.0f | %6.3f %6.3f %6.3f | %.3f(%.0f) %.3f(%.0f) %.3f(%.0f) | %+.1f%% +- %.1f%%\n",
             ana::JetRs[ir], ana::ptBinsUsed[p], ana::ptBinsUsed[p+1], pA[r][p], pTrue, pC[r][p],
             mA, 1000*eA, mCo, 1000*eCo, mS, 1000*eS,
             100*(mCo/mS-1), 100*hypot(eCo, eS)/mS);
    }
  }
  printf("(mean errors in parentheses, x1000; corrected-mean error includes the purity term as stored by purityCorrect)\n");

  // ---------------- drawing ----------------
  string pdf = outdir + "/draw_jet12_purity_corrected_xj.pdf";
  TCanvas * c = new TCanvas("c", "", 1500, 800);
  c->SaveAs((pdf+"[").c_str());
  for (int r = 0; r < nR; r++) {
    int ir = irs[r];
    c->Clear();
    for (int p = 0; p < nPt; p++) {
      c->cd();
      TPad * top = new TPad(Form("ptop%d_%d",r,p), "", p/3.0, 0.35, (p+1)/3.0, 1.0);
      TPad * bot = new TPad(Form("pbot%d_%d",r,p), "", p/3.0, 0.0, (p+1)/3.0, 0.35);
      top->SetLeftMargin(0.16); top->SetRightMargin(0.03); top->SetBottomMargin(0.02); top->SetTopMargin(0.05);
      bot->SetLeftMargin(0.16); bot->SetRightMargin(0.03); bot->SetTopMargin(0.02); bot->SetBottomMargin(0.3);
      top->Draw(); bot->Draw();
      top->cd();
      TH1D * sA = shape(hA[r][p], Form("sA_%d_%d", r, p));
      TH1D * sCo = shape(hCorr[r][p], Form("sCo_%d_%d", r, p));
      TH1D * sS = shape(hSig[r][p], Form("sS_%d_%d", r, p));
      double ymax = std::max({sA->GetMaximum(), sCo->GetMaximum(), sS->GetMaximum()});
      sS->SetLineColor(kGreen+2); sS->SetFillColor(kGreen-9); sS->SetFillStyle(1001); sS->SetMarkerSize(0);
      sA->SetLineColor(kGray+2); sA->SetMarkerColor(kGray+2); sA->SetMarkerStyle(24);
      sCo->SetLineColor(kBlack); sCo->SetMarkerColor(kBlack); sCo->SetMarkerStyle(20);
      sS->GetXaxis()->SetRangeUser(0, 2); sS->GetXaxis()->SetLabelSize(0);
      sS->GetYaxis()->SetTitleSize(0.055); sS->GetYaxis()->SetLabelSize(0.045); sS->GetYaxis()->SetTitleOffset(1.3);
      sS->SetMinimum(0); sS->SetMaximum(1.7*ymax);
      sS->Draw("e2"); sA->Draw("e same"); sCo->Draw("e same");
      TLegend * l = new TLegend(0.62, 0.5, 0.97, 0.78); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.042);
      l->AddEntry(sA, "A, uncorrected", "lp");
      l->AddEntry(sCo, "A, corrected", "lp");
      l->AddEntry(sS, "truth-tagged #gamma", "f");
      l->Draw();
      d.drawAll({"Pythia8 QCD (Jet12_long) MC"}, {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ana::ptBinsUsed[p], ana::ptBinsUsed[p+1]),
                Form("Jet R=%.1f", ana::JetRs[ir]), Form("P_{A} = %.2f, P_{C} = %.2f (ABCD)", pA[r][p], pC[r][p])}, .2, .88, 13, gPad->GetWh()*0.8);
      bot->cd();
      TH1D * rCo = (TH1D*)sCo->Clone(Form("rCo_%d_%d", r, p)); rCo->Divide(sS);
      TH1D * rA  = (TH1D*)sA->Clone(Form("rA_%d_%d", r, p));  rA->Divide(sS);
      sanitize(rCo); sanitize(rA);
      rCo->SetMinimum(0); rCo->SetMaximum(2.5);
      rCo->GetYaxis()->SetTitle("ratio to tagged"); rCo->GetYaxis()->SetNdivisions(505);
      rCo->GetYaxis()->SetTitleSize(0.1); rCo->GetYaxis()->SetLabelSize(0.08); rCo->GetYaxis()->SetTitleOffset(0.7);
      rCo->GetXaxis()->SetTitleSize(0.11); rCo->GetXaxis()->SetLabelSize(0.09);
      rCo->GetXaxis()->SetRangeUser(0, 2);
      rCo->Draw("e"); rA->Draw("e same");
      TLine * one = new TLine(0, 1, 2, 1); one->SetLineStyle(2); one->Draw();
    }
    c->SaveAs(pdf.c_str());
  }
  c->SaveAs((pdf+"]").c_str());

  fout->cd();
  for (int r = 0; r < nR; r++) {
    for (int j = 0; j < 4; j++) hABCD[r][j]->Write();
    for (int p = 0; p < nPt; p++) { hA[r][p]->Write(); hC[r][p]->Write(); hSig[r][p]->Write(); hCorr[r][p]->Write(); }
  }
  fout->Close();
}
