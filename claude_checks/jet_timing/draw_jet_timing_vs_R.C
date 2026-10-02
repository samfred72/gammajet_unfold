#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/pho_object.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/jet_object.h"
// Explicit load; run interpreted, never with ACLiC "+" (sibling libgammajet.so collision).
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Radius dependence of the jet timing window 0 < t_MBD - t_jet < 4 ns (issue 1), all seven
// jet radii in one pass over the Sep-1 data tree. Same selection as draw_jet_timing.C
// (region A, paired per unfolder::check_pair incl. the xJ floor, 15 < pT^gamma < 35 GeV,
// jet pT = jet_pt_calib/jesNominal), evaluated separately for every R - the treemaking
// window is applied per radius, each radius picking its own leading in-window jet.
//
//   p1  Delta t_jet distribution per R (normalized), window drawn
//   p2  mean and RMS of in-window Delta t_jet vs R, and the fraction of jets in the first
//       0.25 ns of the window (the pile-up against the 0 ns edge)
//   p3  Gaussian-model fraction outside [0,4] ns vs jet pT, one curve per R
//       (same truncated-Gaussian model as draw_jet_timing.C p8 - size of effect only)
//   p4  mean jet EM fraction vs R, and mean Delta t_jet vs EM fraction per R
//
// Regions A+B+C+D are used for the pT-sliced loss (jet timing is a jet property, not a
// photon-ID one - more statistics), region A for everything else, as in draw_jet_timing.C.

namespace {
  const double tLow = ana::tlowcut, tHigh = ana::thighcut;
  void sanitize(TH1 * h) {
    for (int b = 0; b <= h->GetNcells(); b++)
      if (!std::isfinite(h->GetBinContent(b)) || !std::isfinite(h->GetBinError(b))) { h->SetBinContent(b, 0); h->SetBinError(b, 0); }
  }
  struct LossEst { double frac, err; };
  LossEst gausLoss(TH1 * h, const char * name) {
    LossEst r{0, 0};
    if (h->GetEntries() < 50) return r;
    TF1 * g = new TF1(name, "gaus", tLow, tHigh);
    g->SetParameters(h->GetMaximum(), h->GetMean(), h->GetRMS());
    h->Fit(g, "QRN0L");
    double m = g->GetParameter(1), s = fabs(g->GetParameter(2));
    auto inside = [&](double mm, double ss) { return 0.5*(TMath::Erf((tHigh-mm)/(sqrt(2)*ss)) - TMath::Erf((tLow-mm)/(sqrt(2)*ss))); };
    r.frac = 1 - inside(m, s);
    double fm = 1 - inside(m + g->GetParError(1), s), fs = 1 - inside(m, s + g->GetParError(2));
    r.err = sqrt(pow(fm-r.frac,2) + pow(fs-r.frac,2));
    return r;
  }
}

void draw_jet_timing_vs_R()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();
  const int nR = ana::nJetR;

  TFile * fin = TFile::Open("/home/samson72/sphnx/gammajet_unfold/trees/gammajet_Data.root", "read");
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

  const double phoLow = 15, phoHigh = 35;
  const int nJpt = 5;
  const double jptEdges[nJpt+1] = {5, 9, 13, 18, 24, 35};
  const int nEmf = 5;
  const double emfEdges[nEmf+1] = {0, 0.3, 0.5, 0.7, 0.85, 1.0};
  TH1D * hDt[nR], * hDtPt[nR][nJpt];
  TProfile * pEmf[nR], * pDtEmf[nR];
  long nA[nR], nEdge[nR];
  for (int ir = 0; ir < nR; ir++) {
    hDt[ir] = new TH1D(Form("hDt_r%d", ir), ";t_{MBD} - t_{jet} [ns];normalized", 60, -1, 5);
    for (int j = 0; j < nJpt; j++) hDtPt[ir][j] = new TH1D(Form("hDtPt_r%d_%d", ir, j), "", 60, -1, 5);
    pEmf[ir] = new TProfile(Form("pEmf_r%d", ir), "", 1, 0, 1);
    pDtEmf[ir] = new TProfile(Form("pDtEmf_r%d", ir), ";jet EM fraction;#LT t_{MBD} - t_{jet} #GT [ns]", nEmf, emfEdges);
    nA[ir] = nEdge[ir] = 0;
  }

  Long64_t n = t->GetEntries();
  for (Long64_t e = 0; e < n; e++) {
    t->GetEntry(e);
    if (fabs(vz) > ana::vzcut) continue;
    if (cluster_pt < phoLow || cluster_pt >= phoHigh) continue;
    pho_object pho(cluster_pt, cluster_e, cluster_eta, cluster_phi,
        cluster_showershape[10], cluster_showershape[11], cluster_time,
        cluster_bdt_scores[9], pho_object::get_showershape(cluster_showershape, cluster_pt));
    int iabcd = ana::findabcdBin(pho.iso4, pho.bdt, 0);
    if (iabcd < 0) continue;
    if (fabs(pho.eta) > ana::etacut) continue;
    int ptbin = ana::findPtBin(pho.pt);
    for (int ir = 0; ir < nR; ir++) {
      if (jet_pt_calib[ir] <= 0) continue;
      float jetPt = jet_pt_calib[ir] / ana::jesNominal[ir];
      jet_object jet(jetPt, jet_e[ir], jet_eta[ir], jet_phi[ir], jet_emfrac[ir], 0, 0, jet_time[ir]);
      if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) continue;
      float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(ana::jet_calib_pt_cut[ir]/ana::ptBins[ptbin])+1];
      if (jetPt/pho.pt < lowbin) continue;
      if (jet.deltaPhi(pho) < ana::oppcut) continue;
      float dt = mbd_time - jet_time[ir];
      int ijp = -1; for (int j = 0; j < nJpt; j++) if (jetPt >= jptEdges[j] && jetPt < jptEdges[j+1]) ijp = j;
      if (ijp >= 0) hDtPt[ir][ijp]->Fill(dt);
      if (iabcd != 0) continue;
      hDt[ir]->Fill(dt);
      pEmf[ir]->Fill(0.5, jet_emfrac[ir]);
      pDtEmf[ir]->Fill(jet_emfrac[ir], dt);
      nA[ir]++;
      if (dt < 0.25) nEdge[ir]++;
    }
  }

  drawer d("pythia", "nominal"); // ctor args only pick files for other helpers; labels come from drawAll
  // ---------------- numbers ----------------
  TGraphErrors * gMean = new TGraphErrors(), * gRms = new TGraphErrors(), * gEdge = new TGraphErrors(), * gEmf = new TGraphErrors();
  TGraphErrors * gLoss[nR];
  printf("\nRegion A paired, %.0f-%.0f GeV. Delta t = t_MBD - t_jet (in-window entries only).\n", phoLow, phoHigh);
  printf("%-4s %7s %7s %6s %9s %8s | Gaussian-model fraction outside [0,4] ns by jet pT (A+B+C+D):", "R", "N", "<dt>", "RMS", "first0.25", "<EMfrac>");
  for (int j = 0; j < nJpt; j++) printf(" %2.0f-%-2.0f", jptEdges[j], jptEdges[j+1]);
  printf("\n");
  for (int ir = 0; ir < nR; ir++) {
    double R = ana::JetRs[ir];
    double fEdge = nA[ir] ? (double)nEdge[ir]/nA[ir] : 0;
    gMean->SetPoint(ir, R, hDt[ir]->GetMean()); gMean->SetPointError(ir, 0, hDt[ir]->GetMeanError());
    gRms->SetPoint(ir, R, hDt[ir]->GetRMS());   gRms->SetPointError(ir, 0, hDt[ir]->GetRMSError());
    gEdge->SetPoint(ir, R, fEdge);             gEdge->SetPointError(ir, 0, nA[ir] ? sqrt(fEdge*(1-fEdge)/nA[ir]) : 0);
    gEmf->SetPoint(ir, R, pEmf[ir]->GetBinContent(1)); gEmf->SetPointError(ir, 0, pEmf[ir]->GetBinError(1));
    printf("%-4.1f %7ld %7.3f %6.3f %9.3f %8.3f |", R, nA[ir], hDt[ir]->GetMean(), hDt[ir]->GetRMS(), fEdge, pEmf[ir]->GetBinContent(1));
    gLoss[ir] = new TGraphErrors();
    for (int j = 0; j < nJpt; j++) {
      LossEst le = gausLoss(hDtPt[ir][j], Form("g_r%d_%d", ir, j));
      gLoss[ir]->SetPoint(j, 0.5*(jptEdges[j]+jptEdges[j+1]) + (ir-3)*0.25, le.frac);
      gLoss[ir]->SetPointError(j, 0, le.err);
      printf(" %.3f", le.frac);
      fflush(stdout);
    }
    printf("\n");
  }
  printf("\nMean in-window Delta t_jet by jet EM fraction (region A):\n%-4s", "R");
  for (int k = 0; k < nEmf; k++) printf("  %.2f-%.2f", emfEdges[k], emfEdges[k+1]);
  printf("\n");
  for (int ir = 0; ir < nR; ir++) {
    printf("%-4.1f", ana::JetRs[ir]);
    for (int k = 0; k < nEmf; k++) printf("  %9.2f", pDtEmf[ir]->GetBinContent(k+1));
    printf("\n");
  }

  // ---------------- drawing ----------------
  string outdir = "/home/samson72/sphnx/gammajet_unfold/claude_checks/jet_timing/pdfs";
  gSystem->mkdir(outdir.c_str(), true);
  string pdf = outdir + "/draw_jet_timing_vs_R.pdf";
  TCanvas * c = new TCanvas("c", "", 700, 600);
  auto pad = [&]() { c->Clear(); c->SetLeftMargin(0.15); c->SetBottomMargin(0.13); c->SetRightMargin(0.05); c->SetTopMargin(0.05); };
  auto legend = [](double x1, double y1, double x2, double y2) { TLegend * l = new TLegend(x1,y1,x2,y2); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.03); return l; };
  int cols[7] = {kRed+1, kOrange+1, kGreen+2, kBlack, kAzure+7, kBlue+1, kMagenta+1};
  string phoLabel = Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", phoLow, phoHigh);
  c->SaveAs((pdf+"[").c_str());

  // p1
  pad();
  double ymax = 0;
  for (int ir = 0; ir < nR; ir++) { if (hDt[ir]->Integral() > 0) hDt[ir]->Scale(1.0/hDt[ir]->Integral(), "width"); sanitize(hDt[ir]); ymax = std::max(ymax, hDt[ir]->GetMaximum()); }
  TLegend * l1 = legend(0.68, 0.5, 0.93, 0.9);
  for (int ir = 0; ir < nR; ir++) {
    hDt[ir]->SetLineColor(cols[ir]); hDt[ir]->SetLineWidth(2);
    hDt[ir]->SetMaximum(1.6*ymax); hDt[ir]->SetMinimum(0);
    hDt[ir]->Draw(ir ? "hist same" : "hist");
    l1->AddEntry(hDt[ir], Form("R = %.1f", ana::JetRs[ir]), "l");
  }
  for (double x : {tLow, tHigh}) { TLine * l = new TLine(x, 0, x, 1.6*ymax); l->SetLineStyle(2); l->SetLineColor(kGray+2); l->Draw(); }
  l1->Draw();
  d.drawAll({"p+p Run24 Data"}, {phoLabel, "Region A, paired", "dashed: treemaking window (data only)"}, .18, .9, 14, gPad->GetWh()*0.8);
  c->SaveAs(pdf.c_str());

  // p2: mean / RMS / edge fraction vs R
  pad();
  { TH1D * fr = new TH1D("frR", ";jet radius R;[ns]  or  fraction #times 10", 1, 0.1, 0.9);
    fr->SetMinimum(0); fr->SetMaximum(3); fr->Draw();
    TGraphErrors * gEdge10 = new TGraphErrors();
    for (int i = 0; i < gEdge->GetN(); i++) { gEdge10->SetPoint(i, gEdge->GetX()[i], 10*gEdge->GetY()[i]); gEdge10->SetPointError(i, 0, 10*gEdge->GetEY()[i]); }
    gMean->SetMarkerStyle(20); gRms->SetMarkerStyle(24); gEdge10->SetMarkerStyle(21); gEdge10->SetMarkerColor(kRed+1); gEdge10->SetLineColor(kRed+1);
    gMean->Draw("p same"); gRms->Draw("p same"); gEdge10->Draw("p same");
    TLegend * l = legend(0.55, 0.72, 0.95, 0.9);
    l->AddEntry(gMean, "#LT t_{MBD} - t_{jet} #GT [ns]", "p"); l->AddEntry(gRms, "RMS [ns]", "p");
    l->AddEntry(gEdge10, "fraction in first 0.25 ns (#times 10)", "p"); l->Draw(); }
  d.drawAll({"p+p Run24 Data"}, {phoLabel, "Region A, paired", "in-window entries only"}, .18, .9, 14, gPad->GetWh()*0.8);
  c->SaveAs(pdf.c_str());

  // p3: loss vs jet pT per R
  pad();
  { TH1D * fr = new TH1D("frL", ";corrected jet p_{T} [GeV];fraction outside [0,4] ns (Gaussian model)", 1, jptEdges[0], jptEdges[nJpt]);
    fr->SetMinimum(0); fr->SetMaximum(1.0); fr->Draw();
    TLegend * l = legend(0.68, 0.5, 0.93, 0.9);
    for (int ir = 0; ir < nR; ir++) {
      gLoss[ir]->SetMarkerStyle(20); gLoss[ir]->SetMarkerColor(cols[ir]); gLoss[ir]->SetLineColor(cols[ir]);
      gLoss[ir]->Draw("pl same"); l->AddEntry(gLoss[ir], Form("R = %.1f", ana::JetRs[ir]), "p");
    }
    l->Draw(); }
  d.drawAll({"p+p Run24 Data"}, {phoLabel, "Regions A+B+C+D, paired", "Gaussian fit inside window", "valid only if edge entries are in-time jets", "R #geq 0.6, 5-9 GeV: fit unstable"}, .18, .9, 14, gPad->GetWh()*0.8);
  c->SaveAs(pdf.c_str());

  // p4: EM fraction
  pad();
  { TH1D * fr = new TH1D("frE", ";jet EM fraction;#LT t_{MBD} - t_{jet} #GT [ns]", 1, 0, 1);
    fr->SetMinimum(0); fr->SetMaximum(4.5); fr->Draw();
    TLegend * l = legend(0.68, 0.5, 0.93, 0.9);
    for (int ir = 0; ir < nR; ir++) {
      sanitize(pDtEmf[ir]);
      pDtEmf[ir]->SetMarkerStyle(20); pDtEmf[ir]->SetMarkerColor(cols[ir]); pDtEmf[ir]->SetLineColor(cols[ir]);
      pDtEmf[ir]->Draw("e same");
      l->AddEntry(pDtEmf[ir], Form("R = %.1f (#LTEMf#GT %.2f)", ana::JetRs[ir], pEmf[ir]->GetBinContent(1)), "p");
    }
    l->Draw(); }
  d.drawAll({"p+p Run24 Data"}, {phoLabel, "Region A, paired", "in-window entries only"}, .18, .9, 14, gPad->GetWh()*0.8);
  c->SaveAs(pdf.c_str());
  c->SaveAs((pdf+"]").c_str());
}
