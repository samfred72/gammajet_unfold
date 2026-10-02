#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"

// Explicit load - the sibling gammajet project's libgammajet.so has same-named classes.
// Run interpreted (root -b -q draw_timingana_tcal.C), never with ACLiC "+".
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Standard vs Dading Chen's corrected tower times (CaloTowerTimeCalibration sidecar), from
// TimingAna trees made with timecalib = true (SDCC timingana/macros/condor). Every tower has
// both times: *_time_std = standard TOWERINFO_CALIB time (official per-tower mean already
// subtracted), *_time = corrected (sample/group alignment, sector and tower offsets, slew).
//
// Towers used: E > towerEmin (his slew fit is derived above 0.5 GeV), finite corrected time,
// EMCal raw towers excluding zero-suppressed ones (status bit 5, t = 0 by construction).
// The same tower set enters the standard and corrected versions, so differences come from
// the calibration alone. Jet-level times are energy-weighted over those towers.
//
//   p1  tower time distributions per calorimeter, standard vs corrected
//   p2  one-sample-shifted fraction (|t| > 9 ns) vs tower energy, standard vs corrected
//   p3  median tower time vs tower energy (slew), standard vs corrected
//   p4  median t_EMCal - t_HCal vs jet pT, standard vs corrected (MBD-free)
//   p5  median jet time vs jet pT and vs EM fraction, standard vs corrected
//   p6  jet-time width (half 68% interval) vs jet pT: t_jet and t_MBD - t_jet, std vs corrected
//   p7  jet time vs t_MBD, standard and corrected (correlation)

namespace {
  const int ir = 2; // R = 0.4, must match the TimingAna run
  const char *inglob = "/home/samson72/sphnx/gammajet_unfold/claude_checks/jet_timing/timingana/tcal/*.root";
  const char *outpdf = "/home/samson72/sphnx/gammajet_unfold/claude_checks/jet_timing/pdfs/draw_timingana_tcal_r04.pdf";
  const float vzCut = 60;
  const float ptMin = 5;
  const float towerEmin = 0.5; // GeV
  const float shiftCut = 9;    // ns; one ADC sample is 17.6 ns

  void sanitize(TH1 *h) {
    for (int i = 0; i <= h->GetNcells(); i++)
      if (!std::isfinite(h->GetBinContent(i)) || !std::isfinite(h->GetBinError(i))) { h->SetBinContent(i, 0); h->SetBinError(i, 0); }
  }
  void margins(TVirtualPad *p, bool z = false) {
    p->SetLeftMargin(.15); p->SetBottomMargin(.13); p->SetTopMargin(.05); p->SetRightMargin(z ? .15 : .05);
  }
  // per-x-bin quantile of a TH2 as a TGraphErrors (error = half the 16-84% interval / sqrt(n))
  TGraphErrors *quantileGraph(TH2 *h, double q, int minN = 20, bool halfWidth = false) {
    TGraphErrors *g = new TGraphErrors();
    for (int b = 1; b <= h->GetNbinsX(); b++) {
      TH1D *p = h->ProjectionY(Form("%s_py%d", h->GetName(), b), b, b);
      if (p->GetEntries() < minN) { delete p; continue; }
      double qs[3] = {0.16, q, 0.84}, v[3];
      p->GetQuantiles(3, v, qs);
      double w = 0.5 * (v[2] - v[0]);
      int n = g->GetN();
      g->SetPoint(n, h->GetXaxis()->GetBinCenter(b), halfWidth ? w : v[1]);
      g->SetPointError(n, 0, halfWidth ? w / std::sqrt(2. * p->GetEntries()) : w / std::sqrt(p->GetEntries()));
      delete p;
    }
    return g;
  }
}

void draw_timingana_tcal(const char *files = inglob)
{
  gStyle->SetOptStat(0);
  TChain *T = new TChain("T");
  int nfiles = T->Add(files);
  if (nfiles == 0 || T->GetEntries() == 0) { std::cerr << "no input in " << files << std::endl; return; }

  float vz, mbd, pt, eta, emfrac;
  std::vector<int> *tw_calo = nullptr, *em_status = nullptr;
  std::vector<float> *tw_e = nullptr, *tw_t = nullptr, *tw_ts = nullptr, *em_e = nullptr, *em_t = nullptr, *em_ts = nullptr;
  std::vector<bool> *em_good = nullptr;
  T->SetBranchAddress("vz", &vz);
  T->SetBranchAddress("mbd_time", &mbd);
  T->SetBranchAddress("jet_pt_calib", &pt);
  T->SetBranchAddress("jet_eta", &eta);
  T->SetBranchAddress("jet_emfrac", &emfrac);
  T->SetBranchAddress("tw_calo", &tw_calo);
  T->SetBranchAddress("tw_e", &tw_e);
  T->SetBranchAddress("tw_time", &tw_t);
  T->SetBranchAddress("tw_time_std", &tw_ts);
  T->SetBranchAddress("em_e", &em_e);
  T->SetBranchAddress("em_time", &em_t);
  T->SetBranchAddress("em_time_std", &em_ts);
  T->SetBranchAddress("em_isgood", &em_good);
  T->SetBranchAddress("em_status", &em_status);

  const char *cname[3] = {"EMCal", "iHCal", "oHCal"};
  const char *vname[2] = {"standard", "corrected"};
  int vcol[2] = {kRed + 1, kBlue + 1};
  int ccol[3] = {kBlue, kSpring - 1, kMagenta + 1};
  const int nte = 20; double teb[nte + 1];
  for (int i = 0; i <= nte; i++) teb[i] = towerEmin * std::pow(40.0, i / (double) nte); // 0.5-20 GeV
  const int nptb = 9;
  double ptb[nptb + 1] = {5, 7, 9, 11, 13, 15, 18, 22, 28, 40};

  TH1D *hT[3][2];
  TH2D *hTvsE[3][2];
  TProfile *pShift[3][2];
  TH2D *hDiffPt[2], *hJetPt[2], *hJetEf[2], *hJetMbdPt[2], *hJetVsMbd[2];
  for (int v = 0; v < 2; v++) {
    for (int c = 0; c < 3; c++) {
      hT[c][v] = new TH1D(Form("hT%d%d", c, v), ";tower t [ns];normalized counts", 240, -30, 30);
      hTvsE[c][v] = new TH2D(Form("hTvsE%d%d", c, v), ";tower energy [GeV];tower t [ns]", nte, teb, 480, -30, 30);
      pShift[c][v] = new TProfile(Form("pShift%d%d", c, v), Form(";tower energy [GeV];fraction with |t| > %.0f ns", shiftCut), nte, teb);
    }
    hDiffPt[v] = new TH2D(Form("hDiffPt%d", v), ";p_{T}^{jet} [GeV];t_{EMCal} - t_{HCal} [ns]", nptb, ptb, 480, -30, 30);
    hJetPt[v] = new TH2D(Form("hJetPt%d", v), ";p_{T}^{jet} [GeV];t_{jet} [ns]", nptb, ptb, 480, -30, 30);
    hJetEf[v] = new TH2D(Form("hJetEf%d", v), ";jet EM fraction;t_{jet} [ns]", 10, 0, 1, 480, -30, 30);
    hJetMbdPt[v] = new TH2D(Form("hJetMbdPt%d", v), ";p_{T}^{jet} [GeV];t_{MBD} - t_{jet} [ns]", nptb, ptb, 480, -30, 30);
    hJetVsMbd[v] = new TH2D(Form("hJetVsMbd%d", v), ";t_{MBD} [ns];t_{jet} [ns]", 120, -15, 15, 120, -15, 15);
  }

  long nsel = 0;
  for (Long64_t i = 0; i < T->GetEntries(); i++) {
    T->GetEntry(i);
    if (!std::isfinite(mbd) || !std::isfinite(vz) || std::fabs(vz) > vzCut) continue;
    if (std::fabs(eta) > 1.1 - ana::JetRs[ir]) continue;
    double jpt = pt / ana::jesNominal[ir];
    if (jpt < ptMin) continue;
    nsel++;
    // [calo group: 0 = EMCal, 1 = HCal][version] energy-weighted sums
    double tw[2][2] = {{0, 0}, {0, 0}}, te[2] = {0, 0};
    auto fillTower = [&](int c, float e, float tc, float ts) {
      if (e <= towerEmin || !std::isfinite(tc) || !std::isfinite(ts)) return;
      float t[2] = {ts, tc};
      for (int v = 0; v < 2; v++) {
        hT[c][v]->Fill(t[v]);
        hTvsE[c][v]->Fill(e, t[v]);
        pShift[c][v]->Fill(e, std::fabs(t[v]) > shiftCut ? 1 : 0);
        tw[c > 0][v] += e * t[v];
      }
      te[c > 0] += e;
    };
    for (size_t j = 0; j < em_e->size(); j++)
      if (em_good->at(j) && !((em_status->at(j) >> 5) & 1)) fillTower(0, em_e->at(j), em_t->at(j), em_ts->at(j));
    for (size_t j = 0; j < tw_e->size(); j++)
      if (tw_calo->at(j) > 0) fillTower(tw_calo->at(j), tw_e->at(j), tw_t->at(j), tw_ts->at(j));
    for (int v = 0; v < 2; v++) {
      if (te[0] > 0 && te[1] > 0) hDiffPt[v]->Fill(jpt, tw[0][v] / te[0] - tw[1][v] / te[1]);
      if (te[0] + te[1] > 0) {
        double tj = (tw[0][v] + tw[1][v]) / (te[0] + te[1]);
        hJetPt[v]->Fill(jpt, tj);
        hJetEf[v]->Fill(emfrac, tj);
        hJetMbdPt[v]->Fill(jpt, mbd - tj);
        hJetVsMbd[v]->Fill(mbd, tj);
      }
    }
  }
  std::cout << "files: " << nfiles << ", selected jets: " << nsel << " of " << T->GetEntries() << std::endl;

  drawer d("pythia", "nominal"); // ctor args only pick files for other helpers
  std::vector<std::string> samples = {"p+p Run24 Data"};
  std::vector<std::string> base = {Form("p_{T}^{jet} > %.0f GeV, |#eta^{jet}| < %.1f", ptMin, 1.1 - ana::JetRs[ir]),
                                   Form("Jet R=%.1f", ana::JetRs[ir]),
                                   Form("towers E > %.1f GeV, EMCal no ZS", towerEmin)};
  auto label = [&](std::vector<std::string> extra, float x = .18, float y = .9) {
    std::vector<std::string> f = base;
    f.insert(f.end(), extra.begin(), extra.end());
    d.drawAll(samples, f, x, y, 14, gPad->GetWh() * 0.8);
  };
  auto legend = [](double x1, double y1, double x2, double y2) {
    TLegend *l = new TLegend(x1, y1, x2, y2); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(.028); return l;
  };
  auto styleG = [](TGraph *g, int col, int mk) { g->SetLineColor(col); g->SetMarkerColor(col); g->SetMarkerStyle(mk); g->SetMarkerSize(.9); };

  TCanvas *c = new TCanvas("c", "", 800, 700);
  c->Print(Form("%s[", outpdf));

  // p1: per-calorimeter tower time, std vs corrected
  c->Clear(); c->Divide(2, 2);
  for (int k = 0; k < 3; k++) {
    c->cd(k + 1); margins(gPad); gPad->SetLogy();
    TLegend *l = legend(.5, .78, .95, .93);
    double ymax = 0;
    for (int v = 0; v < 2; v++) { sanitize(hT[k][v]); if (hT[k][v]->Integral() > 0) hT[k][v]->Scale(1. / hT[k][v]->Integral()); ymax = std::max(ymax, hT[k][v]->GetMaximum()); }
    for (int v = 0; v < 2; v++) {
      hT[k][v]->SetLineColor(vcol[v]); hT[k][v]->SetLineWidth(2);
      hT[k][v]->GetYaxis()->SetRangeUser(1e-5, ymax * 30);
      hT[k][v]->Draw(v ? "hist same" : "hist");
      double fsh = hT[k][v]->Integral(1, hT[k][v]->FindBin(-shiftCut)) + hT[k][v]->Integral(hT[k][v]->FindBin(shiftCut), hT[k][v]->GetNbinsX());
      l->AddEntry(hT[k][v], Form("%s %s: RMS %.2f, |t|>%.0f: %.3f", cname[k], vname[v], hT[k][v]->GetRMS(), shiftCut, fsh), "l");
    }
    l->Draw();
  }
  c->cd(4); label({}, .1, .85);
  c->Print(outpdf);

  // p2: shifted fraction vs tower energy
  c->Clear(); margins(gPad); gPad->SetLogx();
  TLegend *l2 = legend(.5, .7, .93, .93);
  for (int k = 0; k < 3; k++) for (int v = 0; v < 2; v++) {
    sanitize(pShift[k][v]);
    pShift[k][v]->SetLineColor(ccol[k]); pShift[k][v]->SetMarkerColor(ccol[k]);
    pShift[k][v]->SetMarkerStyle(v ? 20 : 24); pShift[k][v]->SetLineStyle(v ? 1 : 2);
    pShift[k][v]->GetYaxis()->SetRangeUser(0, 0.6);
    pShift[k][v]->Draw(k + v ? "same" : "");
    l2->AddEntry(pShift[k][v], Form("%s %s", cname[k], vname[v]), "lp");
  }
  l2->Draw();
  label({"open: standard, filled: corrected"});
  c->Print(outpdf);

  // p3: median tower time vs energy (slew)
  c->Clear(); margins(gPad); gPad->SetLogx();
  TMultiGraph *mg3 = new TMultiGraph(); mg3->SetTitle(";tower energy [GeV];median tower t [ns]");
  TLegend *l3 = legend(.5, .7, .93, .93);
  for (int k = 0; k < 3; k++) for (int v = 0; v < 2; v++) {
    TGraphErrors *g = quantileGraph(hTvsE[k][v], 0.5);
    styleG(g, ccol[k], v ? 20 : 24); g->SetLineStyle(v ? 1 : 2);
    mg3->Add(g, "pl"); l3->AddEntry(g, Form("%s %s", cname[k], vname[v]), "lp");
  }
  mg3->Draw("a"); mg3->GetYaxis()->SetRangeUser(-8, 6);
  l3->Draw();
  label({"open: standard, filled: corrected"});
  c->Print(outpdf);

  // p4: t_EMCal - t_HCal vs jet pT
  c->Clear(); margins(gPad);
  TMultiGraph *mg4 = new TMultiGraph(); mg4->SetTitle(";p_{T}^{jet} [GeV];median t_{EMCal} - t_{HCal} [ns]");
  TLegend *l4 = legend(.55, .8, .93, .93);
  for (int v = 0; v < 2; v++) {
    TGraphErrors *g = quantileGraph(hDiffPt[v], 0.5); styleG(g, vcol[v], 20);
    mg4->Add(g, "pl"); l4->AddEntry(g, vname[v], "lp");
  }
  mg4->Draw("a"); mg4->GetYaxis()->SetRangeUser(-4, 6);
  l4->Draw();
  label({"no MBD time involved"});
  c->Print(outpdf);

  // p5: median jet time vs pT and vs EM fraction
  c->Clear(); c->Divide(2, 1);
  TH2D *h5[2][2] = {{hJetPt[0], hJetPt[1]}, {hJetEf[0], hJetEf[1]}};
  const char *t5[2] = {";p_{T}^{jet} [GeV];median t_{jet} [ns]", ";jet EM fraction;median t_{jet} [ns]"};
  for (int p = 0; p < 2; p++) {
    c->cd(p + 1); margins(gPad);
    TMultiGraph *mg = new TMultiGraph(); mg->SetTitle(t5[p]);
    TLegend *l = legend(.5, .8, .93, .93);
    for (int v = 0; v < 2; v++) {
      TGraphErrors *g = quantileGraph(h5[p][v], 0.5); styleG(g, vcol[v], 20);
      mg->Add(g, "pl"); l->AddEntry(g, vname[v], "lp");
    }
    mg->Draw("a"); mg->GetYaxis()->SetRangeUser(-4, 6);
    l->Draw();
    if (p == 0) label({}, .2, .4);
  }
  c->Print(outpdf);

  // p6: width vs pT
  c->Clear(); margins(gPad);
  TMultiGraph *mg6 = new TMultiGraph(); mg6->SetTitle(";p_{T}^{jet} [GeV];half 16-84% width [ns]");
  TLegend *l6 = legend(.45, .75, .93, .93);
  for (int v = 0; v < 2; v++) {
    TGraphErrors *g1 = quantileGraph(hJetPt[v], 0.5, 20, true); styleG(g1, vcol[v], 20);
    TGraphErrors *g2 = quantileGraph(hJetMbdPt[v], 0.5, 20, true); styleG(g2, vcol[v], 24); g2->SetLineStyle(2);
    mg6->Add(g1, "pl"); mg6->Add(g2, "pl");
    l6->AddEntry(g1, Form("t_{jet}, %s", vname[v]), "lp");
    l6->AddEntry(g2, Form("t_{MBD} - t_{jet}, %s", vname[v]), "lp");
  }
  mg6->Draw("a"); mg6->GetYaxis()->SetRangeUser(0, 5);
  l6->Draw();
  label({}, .18, .35);
  c->Print(outpdf);

  // p7: jet time vs MBD time
  c->Clear(); c->Divide(2, 1);
  for (int v = 0; v < 2; v++) {
    c->cd(v + 1); margins(gPad, true); gPad->SetLogz();
    sanitize(hJetVsMbd[v]);
    hJetVsMbd[v]->Draw("colz");
    d.drawAll(samples, {vname[v], Form("corr. = %.2f", hJetVsMbd[v]->GetCorrelationFactor())}, .2, .9, 14, gPad->GetWh() * 0.8);
  }
  c->Print(outpdf);
  c->Print(Form("%s]", outpdf));

  // numeric summary
  printf("\ntower time RMS / shifted fraction (|t| > %.0f ns), E > %.1f GeV\n", shiftCut, towerEmin);
  for (int k = 0; k < 3; k++) for (int v = 0; v < 2; v++) {
    TH1D *h = hT[k][v]; // already unit-normalized
    double fsh = h->Integral(1, h->FindBin(-shiftCut)) + h->Integral(h->FindBin(shiftCut), h->GetNbinsX());
    printf("  %-6s %-9s RMS %5.2f  shifted %.3f  (n=%.0f)\n", cname[k], vname[v], h->GetRMS(), fsh, hTvsE[k][v]->GetEntries());
  }
  printf("\njet time width (half 16-84%%) vs pT:  pT  | t_jet std  t_jet corr | t_MBD-t_jet std  corr\n");
  TGraphErrors *w[4] = {quantileGraph(hJetPt[0], .5, 20, true), quantileGraph(hJetPt[1], .5, 20, true),
                        quantileGraph(hJetMbdPt[0], .5, 20, true), quantileGraph(hJetMbdPt[1], .5, 20, true)};
  for (int i = 0; i < w[0]->GetN(); i++) {
    printf("  %5.1f |", w[0]->GetX()[i]);
    for (int k = 0; k < 4; k++) printf(" %8.2f", i < w[k]->GetN() ? w[k]->GetY()[i] : -1);
    printf("\n");
  }
  TGraphErrors *dd[2] = {quantileGraph(hDiffPt[0], .5), quantileGraph(hDiffPt[1], .5)};
  printf("\nmedian t_EMCal - t_HCal vs pT (std / corr):");
  for (int i = 0; i < dd[0]->GetN(); i++) printf("  %.0f: %.2f / %.2f", dd[0]->GetX()[i], dd[0]->GetY()[i], i < dd[1]->GetN() ? dd[1]->GetY()[i] : -99);
  printf("\njet vs MBD time correlation (std / corr): %.2f / %.2f\n", hJetVsMbd[0]->GetCorrelationFactor(), hJetVsMbd[1]->GetCorrelationFactor());
}
