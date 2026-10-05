#include "../../src/ana.h"
#include "../../src/drawer.h"
#include "../../src/pho_object.h"
#include "../../src/jet_object.h"
#include "../../src/reweight_utility.h"
// Explicit load; run interpreted, never with ACLiC "+" (sibling libgammajet.so collision).
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Three-jet veto study: Delta R between the recoil (leading) jet and the stored "third
// jet", for region-A events that pass the NOMINAL selection, all jet radii, Data vs
// Pythia8 gamma+jet.
//
// Selection = unfolder.cc's nominal reco branch (same as purity/purity_vs_eta.C):
//   |vz| < ana::vzcut; photon pT in the used bins (15-35 GeV); region A of
//   ana::findabcdBin(iso4, bdt, 0); unfolder::check_pair (xJ floor, |eta^gamma| <
//   ana::photonEtaCut, |eta^jet| < etacut - R, dphi > oppcut); jet pT > jet_calib_pt_cut.
//   Data jet pT = jet_pt_calib / ana::jesNominal; MC jet pT = jet_pt_smear_truth, photon
//   pT EMR-smeared (ana::emrNominal), weights = drawer.h scalemap x Reweighter (vz, pT),
//   truth-photon windows from treeuser.h (Photon5/10/20 stitching).
//
// Treemaking definitions (local copy ~/sphnx/claude/files/CaloAna.cc, may differ from SDCC):
//   hasthirdjet = npassingjets > 2, counted over ALL jets (including the photon's own jet,
//     before the timing cut, no eta cut) passing raw pT > 3 GeV OR any calibrated /
//     smeared pT > 5 GeV (CaloAna.cc:336-344, 413).
//   thirdjet_dr = DeltaR(leading jet, the jet that was the leader just before the final
//     leader was found) (CaloAna.cc:391-416). If the final leader was the first jet
//     encountered, the third-jet position is never set (thirdjet_pt = 0) and thirdjet_dr
//     is the distance to (eta, phi) = (0, 0) - meaningless. Only thirdjet_pt > 0 entries
//     are histogrammed; the rest are counted and reported as "no third jet recorded".
//   thirdjet_pt is the raw (uncalibrated) pT.
//
// Output: claude_checks/threejet/pdfs/draw_thirdjet_dr.pdf (+ .root)

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
}

void draw_thirdjet_dr()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  const int nR = ana::nJetR;

  // [0] Data, [1] Pythia8
  TH1D * hdr[2][nR];
  double nSel[2][nR] = {}, nVeto[2][nR] = {}, nNoRec[2][nR] = {};
  double nSel2[2][nR] = {}, nVeto2[2][nR] = {}; // sum of w^2 for errors
  for (int s = 0; s < 2; s++) for (int ir = 0; ir < nR; ir++)
    hdr[s][ir] = new TH1D(Form("hdr_%s_%s", s ? "pythia" : "data", ana::rnames[ir]),
        ";#DeltaR(leading jet, third jet);(1/N) dN/d#DeltaR", 45, 0, 4.5);

  auto process = [&](const char * fname, bool isMC, double scale, float truthLo, float truthHi) {
    int s = isMC ? 1 : 0;
    TFile * f = TFile::Open(fname, "read");
    TTree * t = (TTree*)f->Get("towerntup");
    float vz, cluster_pt, cluster_e, cluster_eta, cluster_phi, cluster_time;
    float cluster_showershape[12], cluster_bdt_scores[11];
    float truth_cluster_pt = 0;
    float jet_pt_calib[7], jet_pt_smear_truth[7], jet_e[7], jet_eta[7], jet_phi[7], jet_emfrac[7];
    float thirdjet_pt[7], thirdjet_dr[7];
    bool hasthirdjet[7];
    t->SetBranchStatus("*", 0);
    vector<const char*> br = {"vz","cluster_pt","cluster_e","cluster_eta","cluster_phi","cluster_time","cluster_showershape",
                              "cluster_bdt_scores","jet_e","jet_eta","jet_phi","jet_emfrac","hasthirdjet","thirdjet_pt","thirdjet_dr"};
    br.push_back(isMC ? "jet_pt_smear_truth" : "jet_pt_calib");
    if (isMC) br.push_back("truth_cluster_pt");
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
    t->SetBranchAddress("hasthirdjet", hasthirdjet);
    t->SetBranchAddress("thirdjet_pt", thirdjet_pt);
    t->SetBranchAddress("thirdjet_dr", thirdjet_dr);
    if (isMC) {
      t->SetBranchAddress("jet_pt_smear_truth", jet_pt_smear_truth);
      t->SetBranchAddress("truth_cluster_pt", &truth_cluster_pt);
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
      int ptbin = ana::findPtBin(pho.pt);
      if (ptbin < ana::firstUsedPtBin || ptbin >= ana::firstUsedPtBin + ana::nPtBinsUsed) continue;
      if (ana::findabcdBin(pho.iso4, pho.bdt, 0) != 0) continue; // region A
      for (int ir = 0; ir < nR; ir++) {
        float jpt = isMC ? jet_pt_smear_truth[ir] : jet_pt_calib[ir]/ana::jesNominal[ir];
        if (!(jpt > ana::jet_calib_pt_cut[ir])) continue;
        jet_object jet(jpt, jet_e[ir], jet_eta[ir], jet_phi[ir], jet_emfrac[ir], 0, 0, 0);
        // unfolder::check_pair
        float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(ana::jet_calib_pt_cut[ir]/ana::ptBins[ptbin])+1];
        if (jet.pt/pho.pt < lowbin) continue;
        if (fabs(pho.eta) > ana::photonEtaCut) continue;
        if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) continue;
        if (jet.deltaPhi(pho) < ana::oppcut) continue;
        nSel[s][ir] += w;  nSel2[s][ir] += w*w;
        if (!hasthirdjet[ir]) continue;
        nVeto[s][ir] += w; nVeto2[s][ir] += w*w;
        if (thirdjet_pt[ir] > 0) hdr[s][ir]->Fill(thirdjet_dr[ir], w);
        else nNoRec[s][ir] += w;
      }
    }
    f->Close();
  };

  const string tdir = ana::path("trees/");
  process((tdir + "gammajet_Data.root").c_str(), false, 1, 0, 0);
  cout << "data done" << endl;
  // drawer.h scalemap[isphoton=1] for sim="pythia", treeuser.h truth-photon windows
  process((tdir + "gammajet_pythia_Photon5.root").c_str(),  true, 146359.3,  0, 12);
  process((tdir + "gammajet_pythia_Photon10.root").c_str(), true, 6944.675, 12, 24);
  process((tdir + "gammajet_pythia_Photon20.root").c_str(), true, 130.4461, 24, 100);
  cout << "MC done" << endl;

  // Vetoed fraction per radius (binomial error; for weighted MC the effective-count form).
  auto frac = [&](int s, int ir, double & f, double & ef) {
    f = nSel[s][ir] > 0 ? nVeto[s][ir]/nSel[s][ir] : 0;
    double neff = nSel2[s][ir] > 0 ? nSel[s][ir]*nSel[s][ir]/nSel2[s][ir] : 0;
    ef = neff > 0 ? sqrt(f*(1-f)/neff) : 0;
  };
  printf("\n  R   | vetoed fraction Data | Pythia8          | vetoed w/o recorded 3rd jet: Data  Pythia8\n");
  TGraphErrors * gf[2];
  for (int s = 0; s < 2; s++) gf[s] = new TGraphErrors();
  for (int ir = 0; ir < nR; ir++) {
    double f0, e0, f1, e1;
    frac(0, ir, f0, e0);
    frac(1, ir, f1, e1);
    gf[0]->SetPoint(ir, ana::JetRs[ir] - 0.008, f0); gf[0]->SetPointError(ir, 0, e0);
    gf[1]->SetPoint(ir, ana::JetRs[ir] + 0.008, f1); gf[1]->SetPointError(ir, 0, e1);
    printf("  %.1f | %.3f +- %.3f        | %.3f +- %.3f    | %.2f  %.2f\n", ana::JetRs[ir], f0, e0, f1, e1,
        nVeto[0][ir] > 0 ? nNoRec[0][ir]/nVeto[0][ir] : 0., nVeto[1][ir] > 0 ? nNoRec[1][ir]/nVeto[1][ir] : 0.);
  }

  // ---------------- drawing ----------------
  string outdir = ana::path("claude_checks/threejet/pdfs");
  gSystem->mkdir(outdir.c_str(), true);
  TFile * fout = TFile::Open((outdir + "/draw_thirdjet_dr.root").c_str(), "recreate");
  drawer d("pythia", "nominal");

  TCanvas * c = new TCanvas("c", "", 2000, 1000);
  c->Divide(4, 2, 0.001, 0.001);
  for (int ir = 0; ir < nR; ir++) {
    c->cd(ir+1);
    gPad->SetTicks();
    gPad->SetLeftMargin(.15);
    gPad->SetRightMargin(.04);
    gPad->SetTopMargin(.05);
    gPad->SetBottomMargin(.13);
    TH1D * sD = shape(hdr[0][ir], Form("%s_shape", hdr[0][ir]->GetName()));
    TH1D * sM = shape(hdr[1][ir], Form("%s_shape", hdr[1][ir]->GetName()));
    sD->SetLineColor(kBlack);   sD->SetMarkerColor(kBlack);   sD->SetMarkerStyle(20); sD->SetLineWidth(2);
    sM->SetLineColor(kAzure+2); sM->SetMarkerColor(kAzure+2); sM->SetMarkerStyle(21); sM->SetLineWidth(2);
    sM->SetMinimum(0);
    sM->SetMaximum(1.8*std::max(sD->GetMaximum(), sM->GetMaximum()));
    sM->GetXaxis()->SetTitleSize(0.05);
    sM->GetYaxis()->SetTitleSize(0.05);
    sM->GetXaxis()->SetLabelSize(0.045);
    sM->GetYaxis()->SetLabelSize(0.045);
    sM->GetYaxis()->SetTitleOffset(1.3);
    sM->Draw("hist e");
    sD->Draw("p e same");
    TLine * lr = new TLine(ana::JetRs[ir], 0, ana::JetRs[ir], sM->GetMaximum()*0.6);
    lr->SetLineStyle(2);
    lr->SetLineColor(kGray+2);
    lr->Draw("same");
    d.drawAll({"p+p Run24 Data"}, {"15 GeV < p_{T}^{#gamma} < 35 GeV", Form("Jet R=%.1f", ana::JetRs[ir]),
        "Region A, nominal selection"}, .2, .9, 16, gPad->GetWh()*gPad->GetHNDC());
    double f0, e0, f1, e1;
    frac(0, ir, f0, e0);
    frac(1, ir, f1, e1);
    TLegend * l = new TLegend(.50, .50, .95, .68);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(43);
    l->SetTextSize(15);
    l->AddEntry(sD, Form("Data (vetoed %.0f%%)", 100*f0), "p");
    l->AddEntry(sM, Form("Pythia8 #gamma+jet (vetoed %.0f%%)", 100*f1), "l");
    l->AddEntry(lr, "#DeltaR = R", "l");
    l->Draw();
    fout->cd();
    hdr[0][ir]->Write();
    hdr[1][ir]->Write();
  }

  // 8th pad: vetoed fraction vs R
  c->cd(8);
  gPad->SetTicks();
  gPad->SetLeftMargin(.15);
  gPad->SetRightMargin(.04);
  gPad->SetTopMargin(.05);
  gPad->SetBottomMargin(.13);
  double fmax = 0;
  for (int s = 0; s < 2; s++) for (int i = 0; i < gf[s]->GetN(); i++) fmax = std::max(fmax, gf[s]->GetY()[i]);
  TH1F * fr = gPad->DrawFrame(0.1, 0, 0.9, 1.2);
  fr->GetXaxis()->SetTitle("Jet R");
  fr->GetYaxis()->SetTitle("Fraction of selected events vetoed");
  fr->GetXaxis()->SetTitleSize(0.05);
  fr->GetYaxis()->SetTitleSize(0.05);
  fr->GetXaxis()->SetLabelSize(0.045);
  fr->GetYaxis()->SetLabelSize(0.045);
  fr->GetYaxis()->SetTitleOffset(1.3);
  gf[0]->SetMarkerStyle(20); gf[0]->SetMarkerColor(kBlack);   gf[0]->SetLineColor(kBlack);
  gf[1]->SetMarkerStyle(21); gf[1]->SetMarkerColor(kAzure+2); gf[1]->SetLineColor(kAzure+2);
  gf[0]->Draw("P same");
  gf[1]->Draw("P same");
  d.drawAll({"p+p Run24 Data"}, {"15 GeV < p_{T}^{#gamma} < 35 GeV", "Region A, nominal selection"},
      .2, .9, 16, gPad->GetWh()*gPad->GetHNDC());
  TLegend * lf = new TLegend(.55, .55, .95, .68);
  lf->SetBorderSize(0);
  lf->SetFillStyle(0);
  lf->SetTextFont(43);
  lf->SetTextSize(15);
  lf->AddEntry(gf[0], "Data", "p");
  lf->AddEntry(gf[1], "Pythia8 #gamma+jet", "p");
  lf->Draw();
  fout->cd();
  gf[0]->Write("gVetoFrac_data");
  gf[1]->Write("gVetoFrac_pythia");

  c->SaveAs((outdir + "/draw_thirdjet_dr.pdf").c_str());
  fout->Close();
}
