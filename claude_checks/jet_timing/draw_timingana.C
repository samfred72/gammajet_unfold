#include "../../src/ana.h"
#include "../../src/drawer.h"

// Explicit load - the sibling gammajet project's libgammajet.so has same-named classes.
// Run interpreted (root -b -q draw_timingana.C), never with ACLiC "+".
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Jet-timing study on the TimingAna tree (SDCC gammajet/timingana, one DST): the leading
// jet with NO timing requirement, with its time split by calorimeter.
//
// "EMCal" below is the energy-weighted time of the RAW EMCal towers under the jet's
// retowers, excluding zero-suppressed towers (status bit 5): ZS towers have no waveform
// fit and carry t = 0 exactly, but RetowerCEMC folds them into the retower time, so the
// retower-based EMCal time (what CaloAna's jet_time uses) is pulled toward t = 0. Page 2
// shows the retower version too.
//
//   p1  Delta t_EMCal vs Delta t_HCal per jet (both calorimeters present)
//   p2  mean Delta t vs jet pT: EMCal, iHCal, oHCal, HCal, whole jet (CaloAna definition)
//   p3  Delta t_EMCal vs jet pT and Delta t_HCal vs jet pT (2D)
//   p4  t_EMCal - t_HCal vs jet pT (independent of the MBD time)
//   p5  mean Delta t vs jet EM fraction, whole jet and per calorimeter
//   p6  tower Delta t vs tower energy per calorimeter (amplitude dependence)
//   p7  raw MBD time: south, north, and the t0-corrected mean (t_MBD above)
//   p8  raw jet-level times (no MBD subtraction): whole jet and per calorimeter, with t_MBD
//   p9  raw tower times (E > 0.1 GeV, no MBD subtraction) per calorimeter, with t_MBD
//
// Delta t = t_MBD - t, the convention of the treemaking window (0 < Delta t < 4 ns).
// Jet-level times are energy-weighted over towers with E > 0.1 GeV, as in CaloAna.
// Jet pT is jet_pt_calib / ana::jesNominal, as in unfolder.cc.

namespace {
  const int ir = 2; // R = 0.4, must match the TimingAna run
  const char *infile = ana::path("claude_checks/jet_timing/timingana/timing_DST_JETCALO_run2pp_ana521_2025p007_v001-00047289-00000_r04.root");
  const char *outpdf = ana::path("claude_checks/jet_timing/pdfs/draw_timingana_r04.pdf");
  const float vzCut = 60;
  const float ptMin = 5;

  // NaN/Inf would break ROOT's axis auto-ranging (blank canvas).
  void sanitize(TH1 *h) {
    for (int i = 0; i <= h->GetNcells(); i++)
      if (!std::isfinite(h->GetBinContent(i)) || !std::isfinite(h->GetBinError(i))) { h->SetBinContent(i, 0); h->SetBinError(i, 0); }
  }
  bool ok(float t) { return std::isfinite(t) && t > -900; }
  void margins(TVirtualPad *p, bool z = false) {
    p->SetLeftMargin(.15); p->SetBottomMargin(.13); p->SetTopMargin(.05); p->SetRightMargin(z ? .15 : .05);
  }
  void window(double x1, double x2, bool vertical) {
    for (double t : {ana::tlowcut, ana::thighcut}) {
      TLine *l = vertical ? new TLine(t, gPad->GetUymin(), t, gPad->GetUymax()) : new TLine(x1, t, x2, t);
      l->SetLineStyle(2); l->SetLineColor(kRed); l->Draw();
    }
  }
}

void draw_timingana()
{
  gStyle->SetOptStat(0);
  TFile *f = TFile::Open(infile);
  if (!f || f->IsZombie()) { std::cerr << "cannot open " << infile << std::endl; return; }
  TTree *T = (TTree *) f->Get("T");

  float mbds, mbdn;
  float vz, mbd, pt, eta, emfrac, tj, tem, tih, toh, th;
  std::vector<int> *tw_calo = nullptr; std::vector<float> *tw_e = nullptr, *tw_time = nullptr;
  std::vector<float> *em_e = nullptr, *em_time = nullptr; std::vector<bool> *em_good = nullptr;
  std::vector<int> *em_status = nullptr;
  T->SetBranchAddress("vz", &vz);
  T->SetBranchAddress("mbd_time", &mbd);
  T->SetBranchAddress("mbd_time_south", &mbds);
  T->SetBranchAddress("mbd_time_north", &mbdn);
  T->SetBranchAddress("jet_pt_calib", &pt);
  T->SetBranchAddress("jet_eta", &eta);
  T->SetBranchAddress("jet_emfrac", &emfrac);
  T->SetBranchAddress("jet_time", &tj);
  T->SetBranchAddress("jet_time_em", &tem);
  T->SetBranchAddress("jet_time_ih", &tih);
  T->SetBranchAddress("jet_time_oh", &toh);
  T->SetBranchAddress("jet_time_hcal", &th);
  T->SetBranchAddress("tw_calo", &tw_calo);
  T->SetBranchAddress("tw_e", &tw_e);
  T->SetBranchAddress("tw_time", &tw_time);
  T->SetBranchAddress("em_e", &em_e);
  T->SetBranchAddress("em_time", &em_time);
  T->SetBranchAddress("em_isgood", &em_good);
  T->SetBranchAddress("em_status", &em_status);

  const int nptb = 9;
  double ptb[nptb + 1] = {5, 7, 9, 11, 13, 15, 18, 22, 28, 40};
  TH2D *hEmVsH = new TH2D("hEmVsH", ";#Deltat_{HCal} [ns];#Deltat_{EMCal} [ns]", 120, -15, 15, 120, -15, 15);
  TH2D *hEmVsPt = new TH2D("hEmVsPt", ";p_{T}^{jet} [GeV];#Deltat_{EMCal} [ns]", nptb, ptb, 120, -15, 15);
  TH2D *hHVsPt = new TH2D("hHVsPt", ";p_{T}^{jet} [GeV];#Deltat_{HCal} [ns]", nptb, ptb, 120, -15, 15);
  TH2D *hDiffVsPt = new TH2D("hDiffVsPt", ";p_{T}^{jet} [GeV];t_{EMCal} - t_{HCal} [ns]", nptb, ptb, 120, -15, 15);
  const int np = 6;
  const char *pname[np] = {"jet (CaloAna def.)", "EMCal (raw, no ZS)", "iHCal", "oHCal", "HCal (iH+oH)", "EMCal retowers"};
  int pcol[np] = {kBlack, kBlue, kSpring - 1, kMagenta + 1, kOrange + 7, kAzure + 8};
  TProfile *pPt[np], *pEf[np];
  for (int k = 0; k < np; k++) {
    pPt[k] = new TProfile(Form("pPt%d", k), ";p_{T}^{jet} [GeV];#LT#Deltat#GT [ns]", nptb, ptb);
    pEf[k] = new TProfile(Form("pEf%d", k), ";jet EM fraction;#LT#Deltat#GT [ns]", 10, 0, 1);
  }
  TProfile *pDiff = new TProfile("pDiff", "", nptb, ptb);
  // raw (absolute) time distributions, 0.25 ns bins
  const char *rawax = ";t [ns];normalized counts";
  TH1D *hMbd[3] = {new TH1D("hMbdS", rawax, 160, -20, 20), new TH1D("hMbdN", rawax, 160, -20, 20), new TH1D("hMbd", rawax, 160, -20, 20)};
  const char *mbdname[3] = {"MBD south", "MBD north", "t_{MBD} = (S+N)/2 - t_{0}"};
  TH1D *hJetRaw[np], *hTwRaw[3];
  for (int k = 0; k < np; k++) hJetRaw[k] = new TH1D(Form("hJetRaw%d", k), rawax, 160, -20, 20);
  for (int k = 0; k < 3; k++) hTwRaw[k] = new TH1D(Form("hTwRaw%d", k), rawax, 160, -20, 20);
  const char *tname[3] = {"EMCal (raw, no ZS)", "iHCal", "oHCal"};
  // log-spaced tower energy axis from 0.1 to 20 GeV
  const int nte = 24; double teb[nte + 1];
  for (int i = 0; i <= nte; i++) teb[i] = 0.1 * std::pow(200.0, i / (double) nte);
  TProfile *pTw[3];
  for (int k = 0; k < 3; k++) pTw[k] = new TProfile(Form("pTw%d", k), ";tower energy [GeV];#LT#Deltat_{tower}#GT [ns]", nte, teb);

  long nsel = 0;
  for (Long64_t i = 0; i < T->GetEntries(); i++) {
    T->GetEntry(i);
    if (!std::isfinite(mbd) || !std::isfinite(vz) || std::fabs(vz) > vzCut) continue;
    if (std::fabs(eta) > 1.1 - ana::JetRs[ir]) continue;
    double jpt = pt / ana::jesNominal[ir];
    if (jpt < ptMin) continue;
    nsel++;
    double tw = 0, te = 0;
    for (size_t j = 0; j < em_e->size(); j++) {
      bool zs = (em_status->at(j) >> 5) & 1;
      if (!em_good->at(j) || zs || em_e->at(j) <= 0.1) continue;
      tw += em_e->at(j) * em_time->at(j); te += em_e->at(j);
      pTw[0]->Fill(em_e->at(j), mbd - em_time->at(j));
      hTwRaw[0]->Fill(em_time->at(j));
    }
    float temr = tem; // retower-based EMCal time, for comparison only
    tem = te > 0 ? tw / te : -999;
    float t[np] = {tj, tem, tih, toh, th, temr};
    for (int k = 0; k < np; k++) if (ok(t[k])) hJetRaw[k]->Fill(t[k]);
    hMbd[0]->Fill(mbds); hMbd[1]->Fill(mbdn); hMbd[2]->Fill(mbd);
    for (int k = 0; k < np; k++) if (ok(t[k])) { pPt[k]->Fill(jpt, mbd - t[k]); pEf[k]->Fill(emfrac, mbd - t[k]); }
    if (ok(tem)) hEmVsPt->Fill(jpt, mbd - tem);
    if (ok(th)) hHVsPt->Fill(jpt, mbd - th);
    if (ok(tem) && ok(th)) {
      hEmVsH->Fill(mbd - th, mbd - tem);
      hDiffVsPt->Fill(jpt, tem - th);
      pDiff->Fill(jpt, tem - th);
    }
    for (size_t j = 0; j < tw_e->size(); j++)
      if (tw_calo->at(j) > 0 && tw_e->at(j) > 0.1) {
        pTw[tw_calo->at(j)]->Fill(tw_e->at(j), mbd - tw_time->at(j));
        hTwRaw[tw_calo->at(j)]->Fill(tw_time->at(j));
      }
  }
  std::cout << "selected jets: " << nsel << " of " << T->GetEntries() << std::endl;

  for (TH1 *h : std::initializer_list<TH1 *>{hEmVsH, hEmVsPt, hHVsPt, hDiffVsPt, pDiff}) sanitize(h);
  for (int k = 0; k < np; k++) { sanitize(pPt[k]); sanitize(pEf[k]); }
  for (int k = 0; k < 3; k++) sanitize(pTw[k]);

  drawer d("pythia", "nominal"); // ctor args only pick files for other helpers
  std::vector<std::string> samples = {"p+p Run24 Data"};
  std::string rLabel = Form("Jet R=%.1f", ana::JetRs[ir]);
  std::string ptLabel = Form("p_{T}^{jet} > %.0f GeV, |#eta^{jet}| < %.1f", ptMin, 1.1 - ana::JetRs[ir]);
  std::string selLabel = "Leading jet, no timing cut";
  auto label = [&](std::vector<std::string> extra, float x = .18, float y = .9) {
    std::vector<std::string> feats = {selLabel, ptLabel, rLabel};
    feats.insert(feats.end(), extra.begin(), extra.end());
    d.drawAll(samples, feats, x, y, 14, gPad->GetWh() * 0.8);
  };

  TCanvas *c = new TCanvas("c", "", 800, 700);
  c->Print(Form("%s[", outpdf));

  // p1
  margins(gPad, true); gPad->SetLogz();
  hEmVsH->Draw("colz");
  window(-15, 15, true); window(-15, 15, false);
  label({Form("corr. = %.2f", hEmVsH->GetCorrelationFactor()), "dashed: treemaking window"});
  c->Print(outpdf);

  // p2
  c->Clear(); margins(gPad); gPad->SetLogz(0);
  TLegend *lg = new TLegend(.55, .62, .93, .92); lg->SetBorderSize(0); lg->SetFillStyle(0);
  for (int k = 0; k < np; k++) {
    pPt[k]->SetLineColor(pcol[k]); pPt[k]->SetMarkerColor(pcol[k]); pPt[k]->SetMarkerStyle(20); pPt[k]->SetMarkerSize(.8);
    pPt[k]->GetYaxis()->SetRangeUser(-4, 8);
    pPt[k]->Draw(k ? "same" : "");
    lg->AddEntry(pPt[k], pname[k], "lp");
  }
  window(ptb[0], ptb[nptb], false);
  lg->Draw();
  label({"dashed: treemaking window"});
  c->Print(outpdf);

  // p3
  c->Clear(); c->Divide(2, 1);
  c->cd(1); margins(gPad, true); gPad->SetLogz(); hEmVsPt->Draw("colz"); window(ptb[0], ptb[nptb], false); label({}, .2);
  c->cd(2); margins(gPad, true); gPad->SetLogz(); hHVsPt->Draw("colz"); window(ptb[0], ptb[nptb], false);
  c->Print(outpdf);

  // p4
  c->Clear(); margins(gPad, true); gPad->SetLogz();
  hDiffVsPt->Draw("colz");
  pDiff->SetLineColor(kRed); pDiff->SetMarkerColor(kRed); pDiff->SetMarkerStyle(20); pDiff->Draw("same");
  label({"red: mean", "no MBD time involved"});
  c->Print(outpdf);

  // p5
  c->Clear(); margins(gPad); gPad->SetLogz(0);
  TLegend *lg2 = new TLegend(.6, .65, .93, .92); lg2->SetBorderSize(0); lg2->SetFillStyle(0);
  for (int k = 0; k < np - 1; k++) { // retower comparison only on p2
    pEf[k]->SetLineColor(pcol[k]); pEf[k]->SetMarkerColor(pcol[k]); pEf[k]->SetMarkerStyle(20); pEf[k]->SetMarkerSize(.8);
    pEf[k]->GetYaxis()->SetRangeUser(-4, 8);
    pEf[k]->Draw(k ? "same" : "");
    lg2->AddEntry(pEf[k], pname[k], "lp");
  }
  window(0, 1, false);
  lg2->Draw();
  label({"dashed: treemaking window"});
  c->Print(outpdf);

  // p6
  c->Clear(); margins(gPad); gPad->SetLogx();
  TLegend *lg3 = new TLegend(.6, .75, .93, .92); lg3->SetBorderSize(0); lg3->SetFillStyle(0);
  for (int k = 0; k < 3; k++) {
    pTw[k]->SetLineColor(pcol[k + 1]); pTw[k]->SetMarkerColor(pcol[k + 1]); pTw[k]->SetMarkerStyle(20); pTw[k]->SetMarkerSize(.8);
    pTw[k]->GetYaxis()->SetRangeUser(-6, 10);
    pTw[k]->Draw(k ? "same" : "");
    lg3->AddEntry(pTw[k], tname[k], "lp");
  }
  window(teb[0], teb[nte], false);
  lg3->Draw();
  label({"towers in leading jet, E > 0.1 GeV"});
  c->Print(outpdf);

  // p7-p9: raw time distributions, unit area, log y
  auto overlay = [&](std::vector<TH1D *> hs, std::vector<std::string> names, std::vector<int> cols, std::vector<std::string> extra) {
    c->Clear(); margins(gPad); gPad->SetLogx(0); gPad->SetLogy();
    TLegend *l = new TLegend(.6, .92 - .055 * hs.size(), .93, .92); l->SetBorderSize(0); l->SetFillStyle(0);
    double ymax = 0;
    for (auto h : hs) { sanitize(h); if (h->Integral() > 0) h->Scale(1. / h->Integral()); ymax = std::max(ymax, h->GetMaximum()); }
    for (size_t k = 0; k < hs.size(); k++) {
      hs[k]->SetLineColor(cols[k]); hs[k]->SetLineWidth(2);
      hs[k]->GetYaxis()->SetRangeUser(1e-4, ymax * 20);
      hs[k]->Draw(k ? "hist same" : "hist");
      l->AddEntry(hs[k], Form("%s: #mu=%.2f, RMS=%.2f", names[k].c_str(), hs[k]->GetMean(), hs[k]->GetRMS()), "l");
    }
    l->SetTextSize(.025);
    l->Draw();
    label(extra);
    c->Print(outpdf);
  };
  for (int k = 0; k < 3; k++) sanitize(hMbd[k]);
  overlay({hMbd[0], hMbd[1], hMbd[2]}, {mbdname[0], mbdname[1], mbdname[2]}, {kRed, kBlue, kBlack}, {"raw times, no subtraction"});
  overlay({hJetRaw[0], hJetRaw[1], hJetRaw[5], hJetRaw[2], hJetRaw[3], hJetRaw[4], hMbd[2]},
          {"jet (CaloAna def.)", "EMCal (raw, no ZS)", "EMCal retowers", "iHCal", "oHCal", "HCal (iH+oH)", "t_{MBD}"},
          {kBlack, kBlue, kAzure + 8, kSpring - 1, kMagenta + 1, kOrange + 7, kRed}, {"jet-level times, no MBD subtraction"});
  overlay({hTwRaw[0], hTwRaw[1], hTwRaw[2], hMbd[2]}, {"EMCal towers (no ZS)", "iHCal towers", "oHCal towers", "t_{MBD}"},
          {kBlue, kSpring - 1, kMagenta + 1, kRed}, {"towers in leading jet, E > 0.1 GeV", "no MBD subtraction"});

  c->Print(Form("%s]", outpdf));

  // numeric summary (goes to draw_timingana_r04.log)
  printf("\nmean Delta t [ns] (entries) vs jet pT\n%-10s", "pT");
  for (int k = 0; k < np; k++) printf("%20s", pname[k]);
  for (int b = 1; b <= nptb; b++) {
    printf("\n%4.0f-%-5.0f", ptb[b - 1], ptb[b]);
    for (int k = 0; k < np; k++) printf("%13.2f(%5.0f)", pPt[k]->GetBinContent(b), pPt[k]->GetBinEntries(b));
  }
  printf("\n\nmean Delta t [ns] (entries) vs jet EM fraction\n%-10s", "EMfrac");
  for (int k = 0; k < np; k++) printf("%20s", pname[k]);
  for (int b = 1; b <= 10; b++) {
    printf("\n%.1f-%.1f   ", (b - 1) / 10., b / 10.);
    for (int k = 0; k < np; k++) printf("%13.2f(%5.0f)", pEf[k]->GetBinContent(b), pEf[k]->GetBinEntries(b));
  }
  printf("\n\nt_EMCal - t_HCal mean vs pT:");
  for (int b = 1; b <= nptb; b++) printf(" %.2f", pDiff->GetBinContent(b));
  printf("\nEMCal vs HCal correlation (p1): %.2f\n", hEmVsH->GetCorrelationFactor());
  printf("\ntower <Delta t> vs tower E:\n%-14s%20s%20s%20s", "E [GeV]", tname[0], tname[1], tname[2]);
  for (int b = 1; b <= nte; b++) {
    printf("\n%5.2f-%-7.2f", teb[b - 1], teb[b]);
    for (int k = 0; k < 3; k++) printf("%13.2f(%5.0f)", pTw[k]->GetBinContent(b), pTw[k]->GetBinEntries(b));
  }
  printf("\n\nraw time distributions: mean / RMS [ns]\n");
  for (int k = 0; k < 3; k++) printf("%-28s %6.2f / %5.2f\n", mbdname[k], hMbd[k]->GetMean(), hMbd[k]->GetRMS());
  for (int k = 0; k < np; k++) printf("jet-level %-18s %6.2f / %5.2f\n", pname[k], hJetRaw[k]->GetMean(), hJetRaw[k]->GetRMS());
  for (int k = 0; k < 3; k++) printf("towers %-21s %6.2f / %5.2f\n", tname[k], hTwRaw[k]->GetMean(), hTwRaw[k]->GetRMS());
}
