#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/pho_object.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/jet_object.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/reweight_utility.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/insitu_utility.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/unfold_utility.h"
// Reuse the production purity solver itself (leakage-corrected ABCD quadratic + bootstrap).
#include "/home/samson72/sphnx/gammajet_unfold/macros/puritymaker.C"
// Explicit load; run interpreted, never with ACLiC "+" (sibling libgammajet.so collision).
R__LOAD_LIBRARY(libgammajet_unfold.so);

// PPG18 review round 2, issue 10: in-situ JES factor p_a with the nominal photon acceptance
// (|eta^gamma| < 1.1) vs a PPG12-style |eta^gamma| < 0.7 cut, all seven jet radii.
//
// The production insitu trees (insitu/inputs/*_insitu.root) carry no photon eta, so every
// scan input is rebuilt here from the trees, with the cut applied consistently everywhere:
//   data A/C      trees/gammajet_Data.root, jet pT = jet_pt_calib (raw, no in-situ factor),
//                 check_pair cuts WITHOUT the xJ floor (the scan applies the floor on the
//                 scaled xJ itself; no pre-cut, so the scan can go below 0.90)
//   MC reference  Photon5+10+20 region A, unfolder.cc MC branch (EMR-smeared photon pT,
//                 jet_pt_smear_truth, treeuser.h truth windows, mcWeight x scalemap), then
//                 insitu_utility::referenceMeans' weighted mean above lowXjFloor
//   purity        nominal: ana::getPurity/getPurityC (committed); |eta| < 0.7:
//                 puritymaker.C::combine_hists on data ABCD counts + truth-matched Photon MC
//                 leakage templates, both with the cut (paired as in unfolder.cc, i.e. with
//                 the nominal floor at the analysis jet scale jet_pt_calib/jesNominal)
//   chi2          grid_insitu.C's purity-corrected chi2, verbatim
// Validation (nominal): rebuilt data counts with raw xJ >= scanLow x floor and the MC
// reference means must equal insitu_utility::cacheDataEvents / referenceMeans on the
// committed insitu inputs.

namespace {
  const int nR = ana::nJetR;
  const int nPt = ana::nPtBinsUsed;
  const int nSel = 2;
  const double etaMax[nSel] = {ana::etacut, 0.7};
  struct Ev { float pho_pt, jet_pt; int ptbin; };
  const float paLo = 0.80, paHi = 1.05; const int nPa = 2500;
}

void insitu_eta07()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();

  vector<Ev> evA[nSel][nR], evC[nSel][nR];
  // MC reference accumulators [sel][r][pt]
  double sw[nSel][nR][nPt] = {}, sw2[nSel][nR][nPt] = {}, swx[nSel][nR][nPt] = {}, swx2[nSel][nR][nPt] = {};
  // purity inputs for the |eta|<0.7 selection [r][abcd]
  TH1D * hD[nR][4], * hM[nR][4];
  for (int ir = 0; ir < nR; ir++) for (int k = 0; k < 4; k++) {
    hD[ir][k] = new TH1D(Form("d07_r%d_%d", ir, k), ";p_{T};Counts", ana::nPtBins, ana::ptBins);
    hM[ir][k] = new TH1D(Form("m07_r%d_%d", ir, k), ";p_{T};Counts", ana::nPtBins, ana::ptBins);
  }
  float lowXj[nR][nPt];
  for (int ir = 0; ir < nR; ir++) for (int p = 0; p < nPt; p++) lowXj[ir][p] = insitu_utility::lowXjFloor(ir, ana::ptBinsUsed[p]);

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
    for (Long64_t e = 0; e < t->GetEntries(); e++) {
      t->GetEntry(e);
      if (fabs(vz) > ana::vzcut) continue;
      if (isMC && !(truth_cluster_pt > truthLo && truth_cluster_pt < truthHi)) continue;
      float w = isMC ? scale*rw.GetWeight(vz, cluster_pt) : 1.0;
      float recoPt = isMC ? cluster_pt + rand.Gaus(0, ana::emResolutionSigma(truth_cluster_pt, ana::emrNominal)*truth_cluster_pt) : cluster_pt;
      pho_object pho(recoPt, cluster_e, cluster_eta, cluster_phi, cluster_showershape[10], cluster_showershape[11],
          cluster_time, cluster_bdt_scores[9], pho_object::get_showershape(cluster_showershape, recoPt));
      if (!(pho.pt >= ana::ptBins[0] && pho.pt < ana::ptBins[ana::nPtBins])) continue;
      if (fabs(pho.eta) > ana::etacut) continue;
      int iabcd = ana::findabcdBin(pho.iso4, pho.bdt, 0);
      if (iabcd < 0) continue;
      bool matched = false;
      if (isMC) {
        pho_object tru(truth_cluster_pt, truth_cluster_e, truth_cluster_eta, truth_cluster_phi,
            truth_cluster_iso3, truth_cluster_iso4, 0, 0.99, 2);
        matched = pho.deltaR(tru) < 0.1;
      }
      int ptbin = ana::findPtBin(pho.pt);
      int p = ptbin - ana::firstUsedPtBin;
      bool central = fabs(pho.eta) < 0.7;
      for (int ir = 0; ir < nR; ir++) {
        float jraw = isMC ? jet_pt_smear_truth[ir] : jet_pt_calib[ir];
        if (!(jraw > 0)) continue;
        jet_object jet(jraw, jet_e[ir], jet_eta[ir], jet_phi[ir], jet_emfrac[ir], 0, 0, 0);
        if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) continue;
        if (jet.deltaPhi(pho) < ana::oppcut) continue;
        // (1) purity inputs, |eta|<0.7, paired as in unfolder.cc (analysis jet scale, nominal floor)
        if (central) {
          float jana = isMC ? jraw : jraw/ana::jesNominal[ir];
          float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(ana::jet_calib_pt_cut[ir]/ana::ptBins[ptbin])+1];
          if (jana > ana::jet_calib_pt_cut[ir] && jana/pho.pt >= lowbin) {
            if (!isMC) hD[ir][iabcd]->Fill(pho.pt);
            else if (matched) hM[ir][iabcd]->Fill(pho.pt, w);
          }
        }
        // (2) scan inputs: reported pT bins, regions A/C (data) or A (MC reference)
        if (p < 0 || p >= nPt) continue;
        for (int s = 0; s < nSel; s++) {
          if (fabs(pho.eta) > etaMax[s]) continue;
          if (!isMC) {
            if (iabcd == 0) evA[s][ir].push_back({pho.pt, jraw, p});
            if (iabcd == 2) evC[s][ir].push_back({pho.pt, jraw, p});
          } else if (iabcd == 0) {
            double x = jraw/pho.pt;
            if (x < lowXj[ir][p]) continue;           // referenceMeans' floor
            if (x < insitu_utility::scanLow*lowXj[ir][p]) continue;
            sw[s][ir][p] += w; sw2[s][ir][p] += w*w; swx[s][ir][p] += w*x; swx2[s][ir][p] += w*x*x;
          }
        }
      }
    }
    f->Close();
    cout << fname << " done" << endl;
  };

  const string tdir = "/home/samson72/sphnx/gammajet_unfold/trees/";
  process((tdir + "gammajet_Data.root").c_str(), false, 1, 0, 0);
  process((tdir + "gammajet_pythia_Photon5.root").c_str(),  true, 146359.3,  0, 12);
  process((tdir + "gammajet_pythia_Photon10.root").c_str(), true, 6944.675, 12, 24);
  process((tdir + "gammajet_pythia_Photon20.root").c_str(), true, 130.4461, 24, 100);

  // ---------------- reference means ----------------
  float refMean[nSel][nR][nPt], refErr[nSel][nR][nPt];
  for (int s = 0; s < nSel; s++) for (int ir = 0; ir < nR; ir++) for (int p = 0; p < nPt; p++) {
    double m = sw[s][ir][p] > 0 ? swx[s][ir][p]/sw[s][ir][p] : 0;
    double var = sw[s][ir][p] > 0 ? swx2[s][ir][p]/sw[s][ir][p] - m*m : 0;
    double neff = sw2[s][ir][p] > 0 ? sw[s][ir][p]*sw[s][ir][p]/sw2[s][ir][p] : 0;
    refMean[s][ir][p] = m; refErr[s][ir][p] = neff > 0 ? sqrt(std::max(var, 0.)/neff) : 0;
  }

  // ---------------- validation (nominal) ----------------
  map<int,double> photon_scale = {{5,146359.3},{10,6944.675},{20,130.4461}};
  const char * inDir = "/home/samson72/sphnx/gammajet_unfold/insitu/inputs";
  printf("\nValidation (nominal selection) vs committed insitu inputs:\n");
  for (int ir = 0; ir < nR; ir++) {
    string dataFile = insitu_utility::insituFilename(inDir, "Data", "", "nominal");
    size_t pA = insitu_utility::cacheDataEvents(dataFile.c_str(), 0, ir).size(), pC = insitu_utility::cacheDataEvents(dataFile.c_str(), 2, ir).size();
    long nA = 0, nC = 0;
    for (auto & ev : evA[0][ir]) if (ev.jet_pt/ev.pho_pt >= insitu_utility::scanLow*lowXj[ir][ev.ptbin]) nA++;
    for (auto & ev : evC[0][ir]) if (ev.jet_pt/ev.pho_pt >= insitu_utility::scanLow*lowXj[ir][ev.ptbin]) nC++;
    float rm[nPt], re[nPt];
    insitu_utility::referenceMeans({
        {insitu_utility::insituFilename(inDir, "Photon5",  "pythia", "nominal"), photon_scale[5]},
        {insitu_utility::insituFilename(inDir, "Photon10", "pythia", "nominal"), photon_scale[10]},
        {insitu_utility::insituFilename(inDir, "Photon20", "pythia", "nominal"), photon_scale[20]},
      }, 0, ir, rm, re, lowXj[ir]);
    printf("  R=%.1f  data A %ld (%zu)  C %ld (%zu) | MC ref <xJ> %.4f %.4f %.4f (%.4f %.4f %.4f)\n", ana::JetRs[ir], nA, pA, nC, pC,
           refMean[0][ir][0], refMean[0][ir][1], refMean[0][ir][2], rm[0], rm[1], rm[2]);
  }

  // ---------------- purities ----------------
  float P[nSel][nR][nPt], PC[nSel][nR][nPt];
  drawer d("pythia", "nominal");
  string outdir = "/home/samson72/sphnx/gammajet_unfold/claude_checks/purity/pdfs";
  TFile * fout = TFile::Open((outdir + "/insitu_eta07.root").c_str(), "recreate");
  for (int ir = 0; ir < nR; ir++) {
    for (int p = 0; p < nPt; p++) {
      P[0][ir][p]  = ana::getPurity(ana::ptBinsUsed[p], ana::ptBinsUsed[p+1], "nominal", ir);
      PC[0][ir][p] = ana::getPurityC(ana::ptBinsUsed[p], ana::ptBinsUsed[p+1], "nominal", ir);
    }
    TH1D * fp[4];
    for (int k = 0; k < 4; k++) { fp[k] = (TH1D*)hM[ir][k]->Clone(Form("fp07_r%d_%d", ir, k)); fp[k]->Divide(hM[ir][k], hM[ir][0]); }
    fout->mkdir(Form("purity_eta07_r%d", ir))->cd();
    TGraphAsymmErrors * odC = nullptr;
    TGraphAsymmErrors * od = combine_hists(hD[ir], fp, &odC);
    fout->cd();
    for (int p = 0; p < nPt; p++) {
      double x, y; int ib = ana::firstUsedPtBin + p;
      od->GetPoint(ib, x, y); P[1][ir][p] = y;
      odC->GetPoint(ib, x, y); PC[1][ir][p] = y;
    }
  }

  // ---------------- scan ----------------
  auto chi2At = [&](int s, int ir, float pa) {
    vector<double> sA(nPt,0), sA2(nPt,0), sC(nPt,0), sC2(nPt,0); vector<int> nA(nPt,0), nC(nPt,0);
    for (auto & ev : evA[s][ir]) { float x = (ev.jet_pt/pa)/ev.pho_pt; if (x < lowXj[ir][ev.ptbin]) continue; sA[ev.ptbin] += x; sA2[ev.ptbin] += x*x; nA[ev.ptbin]++; }
    for (auto & ev : evC[s][ir]) { float x = (ev.jet_pt/pa)/ev.pho_pt; if (x < lowXj[ir][ev.ptbin]) continue; sC[ev.ptbin] += x; sC2[ev.ptbin] += x*x; nC[ev.ptbin]++; }
    double chi = 0;
    for (int p = 0; p < nPt; p++) {
      double rmv = refMean[s][ir][p], rme = refErr[s][ir][p];
      if (rmv <= 0 || nA[p] == 0 || nC[p] == 0) continue;
      float cA, cC; unfold_utility::purityCorrectCoeffs(P[s][ir][p], PC[s][ir][p], nA[p], nC[p], cA, cC);
      double Ncorr = cA*nA[p] - cC*nC[p];
      if (Ncorr <= 0) continue;
      double mean = (cA*sA[p] - cC*sC[p])/Ncorr;
      double var = (cA*sA2[p] - cC*sC2[p])/Ncorr - mean*mean;
      double err = sqrt(std::max(var, 0.)/Ncorr);
      double diff = 1 - mean/rmv;
      double errt = sqrt(err*err/(rmv*rmv) + mean*mean*rme*rme/pow(rmv,4));
      if (errt > 0) chi += diff*diff/(errt*errt);
    }
    return chi;
  };
  struct Fit { double pa, lo, hi, chi; };
  auto fitRange = [&](TGraph * g, double a, double b) {
    Fit f{0, 0, 0, 1e30}; int im = -1;
    for (int i = 0; i < g->GetN(); i++) { double x = g->GetX()[i]; if (x < a - 1e-6 || x > b + 1e-6) continue; if (g->GetY()[i] < f.chi) { f.chi = g->GetY()[i]; f.pa = x; im = i; } }
    f.lo = f.hi = 0;
    for (int i = im; i >= 0; i--) if (g->GetY()[i] > f.chi + 1) { f.lo = f.pa - g->GetX()[i]; break; }
    for (int i = im; i < g->GetN(); i++) if (g->GetY()[i] > f.chi + 1) { f.hi = g->GetX()[i] - f.pa; break; }
    return f;
  };
  TGraph * g[nSel][nR];
  Fit fr[nSel][nR], frIn[nSel][nR];
  for (int s = 0; s < nSel; s++) for (int ir = 0; ir < nR; ir++) {
    g[s][ir] = new TGraph();
    for (int ia = 0; ia <= nPa; ia++) { float pa = paLo + ia*(paHi-paLo)/nPa; g[s][ir]->SetPoint(ia, pa, chi2At(s, ir, pa)); }
    fr[s][ir] = fitRange(g[s][ir], paLo, paHi);
    frIn[s][ir] = fitRange(g[s][ir], insitu_utility::scanLow, insitu_utility::scanHigh);
  }

  printf("\nIn-situ p_a (purity-corrected chi2, Delta chi2 = 1 errors). 'range' = production window [%.2f, %.2f].\n",
         insitu_utility::scanLow, insitu_utility::scanHigh);
  printf("%-4s | %-34s | %-34s | %s\n", "R", "|eta^gamma| < 1.1 (nominal)", "|eta^gamma| < 0.7", "shift (full)");
  for (int ir = 0; ir < nR; ir++) {
    printf("%-4.1f | %.3f -%.3f +%.3f (chi2 %5.1f) [range %.3f] | %.3f -%.3f +%.3f (chi2 %5.1f) [range %.3f] | %+.3f   (ana.h jesNominal %.4f)\n",
           ana::JetRs[ir], fr[0][ir].pa, fr[0][ir].lo, fr[0][ir].hi, fr[0][ir].chi, frIn[0][ir].pa,
           fr[1][ir].pa, fr[1][ir].lo, fr[1][ir].hi, fr[1][ir].chi, frIn[1][ir].pa, fr[1][ir].pa - fr[0][ir].pa, ana::jesNominal[ir]);
  }
  printf("\nPurities P_A (nominal -> |eta|<0.7):\n");
  for (int ir = 0; ir < nR; ir++) printf("  R=%.1f  %.3f->%.3f  %.3f->%.3f  %.3f->%.3f\n", ana::JetRs[ir],
      P[0][ir][0], P[1][ir][0], P[0][ir][1], P[1][ir][1], P[0][ir][2], P[1][ir][2]);
  printf("MC reference <xJ> (nominal -> |eta|<0.7):\n");
  for (int ir = 0; ir < nR; ir++) printf("  R=%.1f  %.4f->%.4f  %.4f->%.4f  %.4f->%.4f\n", ana::JetRs[ir],
      refMean[0][ir][0], refMean[1][ir][0], refMean[0][ir][1], refMean[1][ir][1], refMean[0][ir][2], refMean[1][ir][2]);

  // ---------------- drawing ----------------
  string pdf = outdir + "/insitu_eta07.pdf";
  TCanvas * c = new TCanvas("c", "", 800, 600);
  c->SetLeftMargin(0.14); c->SetBottomMargin(0.13); c->SetRightMargin(0.04); c->SetTopMargin(0.05);
  TH1D * frm = new TH1D("frm", ";jet radius R;in-situ JES factor p_{a}", 1, 0.1, 0.9);
  frm->SetMinimum(0.82); frm->SetMaximum(1.02); frm->GetYaxis()->SetTitleOffset(1.3); frm->Draw();
  TBox * edge = new TBox(0.1, 0.82, 0.9, insitu_utility::scanLow); edge->SetFillColor(kGray); edge->SetFillStyle(1001); edge->Draw();
  frm->Draw("axis same");
  TGraphAsymmErrors * gp[nSel];
  int col[nSel] = {kBlack, kRed+1}; int mk[nSel] = {20, 24};
  for (int s = 0; s < nSel; s++) {
    gp[s] = new TGraphAsymmErrors();
    for (int ir = 0; ir < nR; ir++) { gp[s]->SetPoint(ir, ana::JetRs[ir] + (s ? 0.012 : -0.012), fr[s][ir].pa); gp[s]->SetPointError(ir, 0, 0, fr[s][ir].lo, fr[s][ir].hi); }
    gp[s]->SetMarkerStyle(mk[s]); gp[s]->SetMarkerColor(col[s]); gp[s]->SetLineColor(col[s]); gp[s]->SetMarkerSize(1.3);
    gp[s]->Draw("p same");
  }
  TLegend * l = new TLegend(0.55, 0.2, 0.95, 0.36); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.035);
  l->AddEntry(gp[0], "|#eta^{#gamma}| < 1.1 (nominal)", "p"); l->AddEntry(gp[1], "|#eta^{#gamma}| < 0.7", "p");
  l->AddEntry(edge, "below production scan range", "f"); l->Draw();
  d.drawAll({"p+p Run24 Data vs Pythia8 #gamma+jet"}, {"purity-corrected mean-x_{J} fit", "stat. errors (#Delta#chi^{2}=1)", "inputs rebuilt from trees, extended scan"}, .18, .9, 14, gPad->GetWh()*0.8);
  c->SaveAs(pdf.c_str());
  fout->cd();
  for (int s = 0; s < nSel; s++) for (int ir = 0; ir < nR; ir++) g[s][ir]->Write(Form("chi2_sel%d_r%d", s, ir));
  gp[0]->Write("pa_nominal"); gp[1]->Write("pa_eta07");
  fout->Close();
}
