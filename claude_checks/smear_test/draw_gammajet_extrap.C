#include "../../src/ana.h"
#include "../../src/insitu_utility.h"
// Explicit load - see drawing/draw_final_result.C. Run interpreted, never with ACLiC "+".
R__LOAD_LIBRARY(libgammajet_unfold.so);

// gamma+jet counterpart of draw_jet12_extrap.C / draw_jet12_xj.C: the effect of the low-pT
// extrapolation of the JER smearing width on the in-situ JES input (Pythia8 Photon5/10/20, R=0.4).
//
// Width functions (nominal template h_jer_smear_r04_pileup_EMfracJES_nominal, as fill_jet12_extrap.C):
//   tmpl   Interpolate(pt) - the treemaker's (flat below the first bin centre, 5.19 GeV)
//   lin15  template above 15 GeV, its tangent at 15 GeV continued below
//   lin5   template above 5.19 GeV, its first segment continued below
//   calib  no extra smearing (reference)
//
// Inputs: gammajet_inputs/<mode>/Photon{5,10,20}_pythia_nominal_insitu.root, the in-situ trees written by
// unfolder.cc on branch jerextrap with GAMMAJET_JER_EXTRAP=<mode> (all modes with the same reweighting
// file). There the MC jet_pt_smear_truth is re-smeared with the stored deviate recovered as
// z = (smear - calib)/(pt_ref w_tmpl(pt_ref)): pt_ref = pt_calib below 5 GeV (exact, the treemaker only
// truth-matches jets above it), else the leading truth recoil jet's pT if within 0.75R, else pt_calib
// (approximation). gammajet_inputs/bookkeeping.txt has those three counts per mode, for R=0.4 region-A in-situ pairs.
// gammajet_inputs/pa_values.txt has the nominal in-situ p_a (one scan, 5 GeV jet cut) per mode and radius,
// from logs/<mode>/11a_insitu_scan.log of the jerextrap worktree; tmpl is ana::jesNominal (converged).
//
// The MC selection is the in-situ reference's: region A, the used photon-pT bins, the x_J floor, weights
// = cross section x the in-situ tree's weight (insitu_utility::referenceMeans).
//
//   p1  width functions
//   p2  recoil-jet pT, all used photon-pT bins, with ratio to tmpl
//   p3  x_J per photon-pT bin, with ratio to tmpl
//   p4  MC reference <x_J> vs photon pT, and difference from tmpl
//   p5  in-situ p_a vs R per mode
//   p6  pt_ref bookkeeping

namespace {
  const int nMode = 4;
  const char * modes[nMode]  = {"calib", "tmpl", "lin15", "lin5"};
  const char * labels[nMode] = {"no extra smearing", "Default (flat below 5.2 GeV)", "linear below 15 GeV", "linear below 5.2 GeV"};
  const int cols[nMode]      = {kGray+2, kBlack, kBlue+1, kRed+1};
  const int ir = 2; // R = 0.4
  map<int,double> photon_scale = {{5,146359.3},{10,6944.675},{20,130.4461}}; // as grid_insitu.C

  // NaN/Inf poisons ROOT's auto-ranging (CLAUDE.md) - zero such bins before drawing.
  void sanitize(TH1 * h) {
    for (int b = 0; b <= h->GetNcells(); b++)
      if (!std::isfinite(h->GetBinContent(b)) || !std::isfinite(h->GetBinError(b))) { h->SetBinContent(b, 0); h->SetBinError(b, 0); }
  }
  vector<pair<string,double>> samples(const string & mode) {
    string dir = ana::path("claude_checks/smear_test/gammajet_inputs/" + mode);
    vector<pair<string,double>> s;
    for (int p : {5, 10, 20}) s.push_back({dir + Form("/Photon%d_pythia_nominal_insitu.root", p), photon_scale[p]});
    return s;
  }
  // Recoil-jet pT over all used photon-pT bins, same filters as insitu_utility::referenceMeans.
  TH1D * jetPtHist(const string & mode, const float lowXj[]) {
    TH1D * h = new TH1D(Form("hjet_%s", mode.c_str()), ";p_{T}^{jet} [GeV];normalized", 40, 0, 40);
    for (auto & s : samples(mode)) {
      TFile * f = TFile::Open(s.first.c_str(), "READ");
      if (!f || f->IsZombie()) { cout << "WARNING: missing " << s.first << endl; continue; }
      TTree * t = (TTree*)f->Get("insitutree");
      Float_t pho_pt, jet_pt, w; Int_t abcd, evIr;
      t->SetBranchAddress("pho_pt", &pho_pt); t->SetBranchAddress("jet_pt", &jet_pt);
      t->SetBranchAddress("abcd", &abcd); t->SetBranchAddress("weight", &w); t->SetBranchAddress("ir", &evIr);
      for (Long64_t e = 0; e < t->GetEntries(); e++) {
        t->GetEntry(e);
        if (abcd != 0 || evIr != ir) continue;
        int ipt = ana::findPtBin(pho_pt) - ana::firstUsedPtBin;
        if (ipt < 0 || ipt >= ana::nPtBinsUsed) continue;
        if (jet_pt/pho_pt < lowXj[ipt]) continue;
        h->Fill(jet_pt, s.second*w);
      }
      f->Close();
    }
    return h;
  }
  void unitNorm(TH1 * h) { if (h->Integral() > 0) h->Scale(1.0/h->Integral(), "width"); }
  TH1D * ratioTo(TH1D * num, TH1D * den, const char * name) {
    TH1D * r = (TH1D*)num->Clone(name);
    r->Divide(den);
    sanitize(r);
    return r;
  }
  TLegend * legend(double x1, double y1, double x2, double y2) {
    TLegend * l = new TLegend(x1, y1, x2, y2); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.035); return l;
  }
  void label(vector<string> extra, double x = .2, double y = .88, double csize = 600) {
    vector<string> feats = {"15 GeV < p_{T}^{#gamma} < 35 GeV", "Jet R=0.4, Region A, reco level"};
    for (auto & e : extra) feats.push_back(e);
    insitu_utility::drawSPhenixLabel({"Pythia8 #gamma+jet MC"}, feats, x, y, 14, csize);
  }
}

void draw_gammajet_extrap() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();

  float lowXj[ana::nPtBinsUsed];
  for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) lowXj[ipt] = insitu_utility::lowXjFloor(ir, ana::ptBinsUsed[ipt]);

  string outdir = ana::path("claude_checks/smear_test/pdfs");
  gSystem->mkdir(outdir.c_str(), true);
  string pdf = outdir + "/gammajet_extrap_R04.pdf";
  TCanvas * c = new TCanvas("c", "", 800, 700);
  c->SaveAs((pdf + "[").c_str());

  // ---- p1: width functions ----
  {
    TFile f(ana::path("claude_checks/smear_test/jerband_smearing_templates.root"), "read");
    TH1D * hNom = (TH1D*)f.Get("h_jer_smear_r04_pileup_EMfracJES_nominal");
    hNom->SetDirectory(0);
    f.Close();
    const double w15 = hNom->Interpolate(15), s15 = (hNom->Interpolate(15.01) - hNom->Interpolate(14.99)) / 0.02;
    const double c1 = hNom->GetBinCenter(1), w5 = hNom->GetBinContent(1);
    const double s5 = (hNom->GetBinContent(2) - hNom->GetBinContent(1)) / (hNom->GetBinCenter(2) - c1);
    auto width = [&](int m, double pt) {
      if (m == 2) return pt >= 15 ? hNom->Interpolate(pt) : w15 + s15*(pt - 15);
      if (m == 3) return pt >= c1 ? hNom->Interpolate(pt) : w5 + s5*(pt - c1);
      return hNom->Interpolate(pt);
    };
    c->Clear(); c->SetLeftMargin(0.14); c->SetBottomMargin(0.12); c->SetRightMargin(0.04); c->SetTopMargin(0.05);
    TH1D * frame = new TH1D("frame", ";p_{T}^{ref} [GeV];extra JER width #sigma/p_{T}", 100, 0, 30);
    frame->SetMinimum(0); frame->SetMaximum(1.2); frame->Draw("axis");
    TLegend * l = legend(.5, .62, .93, .78);
    for (int m = 1; m < nMode; m++) {
      TGraph * g = new TGraph();
      for (int i = 0; i <= 300; i++) { double pt = 0.1*i; g->SetPoint(i, pt, width(m, pt)); }
      g->SetLineColor(cols[m]); g->SetLineWidth(3); g->SetLineStyle(m == 1 ? 1 : 2);
      g->Draw("l same");
      l->AddEntry(g, labels[m], "l");
    }
    l->Draw();
    TLine * l5 = new TLine(5, 0, 5, 1.2); l5->SetLineStyle(3); l5->Draw();
    insitu_utility::drawSPhenixLabel({"Pythia8 #gamma+jet MC"}, {"h_jer_smear_r04_pileup_EMfracJES_nominal", "dotted: 5 GeV jet cut"}, .2, .88, 14, c->GetWh());
    c->SaveAs(pdf.c_str());
  }

  // ---- p2: recoil-jet pT ----
  TH1D * hJet[nMode];
  for (int m = 0; m < nMode; m++) { hJet[m] = jetPtHist(modes[m], lowXj); unitNorm(hJet[m]); sanitize(hJet[m]); }
  {
    c->Clear();
    TPad * top = new TPad("top", "", 0, 0.32, 1, 1); top->SetBottomMargin(0.02); top->SetLeftMargin(0.14); top->SetRightMargin(0.04); top->SetTopMargin(0.05); top->SetLogy(); top->Draw();
    TPad * bot = new TPad("bot", "", 0, 0, 1, 0.32); bot->SetTopMargin(0.02); bot->SetBottomMargin(0.3); bot->SetLeftMargin(0.14); bot->SetRightMargin(0.04); bot->Draw();
    top->cd();
    hJet[1]->SetMaximum(50*hJet[1]->GetMaximum()); hJet[1]->SetMinimum(1e-4);
    hJet[1]->GetXaxis()->SetLabelSize(0);
    TLegend * l = legend(.55, .55, .93, .75);
    for (int m : {1, 0, 2, 3}) {
      hJet[m]->SetLineColor(cols[m]); hJet[m]->SetLineWidth(2);
      hJet[m]->Draw(m == 1 ? "hist" : "hist same");
      l->AddEntry(hJet[m], labels[m], "l");
    }
    l->Draw();
    label({}, .2, .88, top->GetWh());
    bot->cd();
    TH1D * frame = (TH1D*)hJet[1]->Clone("frameJ"); frame->Reset();
    frame->GetYaxis()->SetTitle("ratio to Default"); frame->SetMinimum(0.6); frame->SetMaximum(1.4);
    frame->GetXaxis()->SetLabelSize(0.1); frame->GetXaxis()->SetTitleSize(0.11); frame->GetYaxis()->SetLabelSize(0.08); frame->GetYaxis()->SetTitleSize(0.08); frame->GetYaxis()->SetTitleOffset(0.7);
    frame->Draw("axis");
    for (int m : {0, 2, 3}) { TH1D * r = ratioTo(hJet[m], hJet[1], Form("rJ%d", m)); r->SetLineColor(cols[m]); r->SetMarkerColor(cols[m]); r->SetMarkerStyle(20); r->SetMarkerSize(0.6); r->Draw("same e"); }
    TLine * one = new TLine(0, 1, 40, 1); one->SetLineStyle(2); one->Draw();
    c->SaveAs(pdf.c_str());
  }

  // ---- p3: x_J per photon-pT bin ----
  vector<TH1D*> hXj[nMode];
  for (int m = 0; m < nMode; m++) {
    hXj[m] = insitu_utility::buildMCXjByPtBin(samples(modes[m]), 0, ir, Form("hxj_%s", modes[m]), lowXj);
    for (auto h : hXj[m]) { unitNorm(h); sanitize(h); }
  }
  {
    c->Clear();
    c->Divide(ana::nPtBinsUsed, 1, 0.001, 0.001);
    for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) {
      c->cd(ipt + 1);
      TPad * top = new TPad(Form("t%d", ipt), "", 0, 0.32, 1, 1); top->SetBottomMargin(0.02); top->SetLeftMargin(0.17); top->SetRightMargin(0.03); top->SetTopMargin(0.05); top->Draw();
      TPad * bot = new TPad(Form("b%d", ipt), "", 0, 0, 1, 0.32); bot->SetTopMargin(0.02); bot->SetBottomMargin(0.3); bot->SetLeftMargin(0.17); bot->SetRightMargin(0.03); bot->Draw();
      top->cd();
      hXj[1][ipt]->SetMaximum(1.6*hXj[1][ipt]->GetMaximum()); hXj[1][ipt]->SetMinimum(0);
      hXj[1][ipt]->GetYaxis()->SetTitle("normalized");
      hXj[1][ipt]->GetXaxis()->SetLabelSize(0);
      for (int m : {1, 0, 2, 3}) { hXj[m][ipt]->SetLineColor(cols[m]); hXj[m][ipt]->SetLineWidth(2); hXj[m][ipt]->Draw(m == 1 ? "hist" : "hist same"); }
      insitu_utility::drawSPhenixLabel({"Pythia8 #gamma+jet MC"},
          {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1]), "Jet R=0.4, Region A"}, .22, .88, 10, top->GetWh());
      if (ipt == ana::nPtBinsUsed - 1) {
        TLegend * l = legend(.35, .45, .97, .65); l->SetTextSize(0.045);
        for (int m : {1, 0, 2, 3}) l->AddEntry(hXj[m][ipt], labels[m], "l");
        l->Draw();
      }
      bot->cd();
      TH1D * frame = (TH1D*)hXj[1][ipt]->Clone(Form("fx%d", ipt)); frame->Reset();
      frame->GetYaxis()->SetTitle("ratio to Default"); frame->SetMinimum(0.5); frame->SetMaximum(1.5);
      frame->GetXaxis()->SetLabelSize(0.1); frame->GetXaxis()->SetTitleSize(0.11); frame->GetYaxis()->SetLabelSize(0.08); frame->GetYaxis()->SetTitleSize(0.08); frame->GetXaxis()->SetNdivisions(505);
      frame->Draw("axis");
      for (int m : {0, 2, 3}) { TH1D * r = ratioTo(hXj[m][ipt], hXj[1][ipt], Form("rx%d_%d", m, ipt)); r->SetLineColor(cols[m]); r->SetMarkerColor(cols[m]); r->SetMarkerStyle(20); r->SetMarkerSize(0.5); r->Draw("same e"); }
      TLine * one = new TLine(frame->GetXaxis()->GetXmin(), 1, frame->GetXaxis()->GetXmax(), 1); one->SetLineStyle(2); one->Draw();
    }
    c->SaveAs(pdf.c_str());
  }

  // ---- p4: reference <x_J> ----
  float mean[nMode][ana::nPtBinsUsed], err[nMode][ana::nPtBinsUsed];
  for (int m = 0; m < nMode; m++) insitu_utility::referenceMeans(samples(modes[m]), 0, ir, mean[m], err[m], lowXj);
  printf("MC reference <xJ> (R=0.4, region A)\n%-6s", "");
  for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) printf("   %2.0f-%2.0f GeV        ", ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1]);
  printf("\n");
  for (int m = 0; m < nMode; m++) {
    printf("%-6s", modes[m]);
    for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) printf("  %.4f+-%.4f (%+.4f)", mean[m][ipt], err[m][ipt], mean[m][ipt] - mean[1][ipt]);
    printf("\n");
  }
  {
    c->Clear();
    TPad * top = new TPad("top4", "", 0, 0.38, 1, 1); top->SetBottomMargin(0.02); top->SetLeftMargin(0.14); top->SetRightMargin(0.04); top->SetTopMargin(0.05); top->Draw();
    TPad * bot = new TPad("bot4", "", 0, 0, 1, 0.38); bot->SetTopMargin(0.02); bot->SetBottomMargin(0.25); bot->SetLeftMargin(0.14); bot->SetRightMargin(0.04); bot->Draw();
    top->cd();
    TH1D * fr = new TH1D("fr4", ";p_{T}^{#gamma} [GeV];MC #LT x_{J#gamma} #GT", 1, ana::ptBinsUsed[0], ana::ptBinsUsed[ana::nPtBinsUsed]);
    double lo = 1e9, hi = -1e9;
    for (int m = 0; m < nMode; m++) for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) { lo = std::min<double>(lo, mean[m][ipt]); hi = std::max<double>(hi, mean[m][ipt]); }
    fr->SetMinimum(lo - 0.1*(hi - lo)); fr->SetMaximum(hi + 0.8*(hi - lo)); fr->GetXaxis()->SetLabelSize(0); fr->Draw("axis");
    TLegend * l = legend(.55, .55, .93, .75);
    bot->cd();
    TH1D * frb = new TH1D("frb4", ";p_{T}^{#gamma} [GeV];#Delta#LT x_{J#gamma}#GT vs Default", 1, ana::ptBinsUsed[0], ana::ptBinsUsed[ana::nPtBinsUsed]);
    double dmax = 0.005;
    for (int m = 0; m < nMode; m++) for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) dmax = std::max<double>(dmax, fabs(mean[m][ipt] - mean[1][ipt]) + err[m][ipt]);
    frb->SetMinimum(-1.2*dmax); frb->SetMaximum(1.2*dmax);
    frb->GetXaxis()->SetLabelSize(0.08); frb->GetXaxis()->SetTitleSize(0.09); frb->GetYaxis()->SetLabelSize(0.07); frb->GetYaxis()->SetTitleSize(0.07); frb->GetYaxis()->SetTitleOffset(0.9);
    frb->Draw("axis");
    TLine * zero = new TLine(ana::ptBinsUsed[0], 0, ana::ptBinsUsed[ana::nPtBinsUsed], 0); zero->SetLineStyle(2); zero->Draw();
    for (int m = 0; m < nMode; m++) {
      TGraphErrors * g = new TGraphErrors(), * d = new TGraphErrors();
      for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) {
        double x = 0.5*(ana::ptBinsUsed[ipt] + ana::ptBinsUsed[ipt+1]) + 0.4*(m - 1.5), ex = 0;
        g->SetPoint(ipt, x, mean[m][ipt]); g->SetPointError(ipt, ex, err[m][ipt]);
        d->SetPoint(ipt, x, mean[m][ipt] - mean[1][ipt]); d->SetPointError(ipt, ex, err[m][ipt]);
      }
      for (auto gg : {g, d}) { gg->SetMarkerStyle(20); gg->SetMarkerColor(cols[m]); gg->SetLineColor(cols[m]); }
      top->cd(); g->Draw("p same"); l->AddEntry(g, labels[m], "p");
      if (m != 1) { bot->cd(); d->Draw("p same"); }
    }
    top->cd(); l->Draw();
    label({"error bars: MC statistics (Kish)"}, .2, .88, top->GetWh());
    c->SaveAs(pdf.c_str());
  }

  // ---- p5: in-situ p_a per mode ----
  {
    map<string, vector<double>> pa, eLo, eHi;
    pa["tmpl"].assign(ana::jesNominal, ana::jesNominal + ana::nJetR);
    eLo["tmpl"].assign(ana::jesStatErrLow, ana::jesStatErrLow + ana::nJetR);
    eHi["tmpl"].assign(ana::jesStatErrHigh, ana::jesStatErrHigh + ana::nJetR);
    std::ifstream in(ana::path("claude_checks/smear_test/gammajet_inputs/pa_values.txt"));
    string mode;
    while (in >> mode) {
      for (int r = 0; r < ana::nJetR; r++) { double a, lo, hi; in >> a >> lo >> hi; pa[mode].push_back(a); eLo[mode].push_back(lo); eHi[mode].push_back(hi); }
    }
    c->Clear(); c->SetLeftMargin(0.14); c->SetBottomMargin(0.12); c->SetRightMargin(0.04); c->SetTopMargin(0.05);
    TH1D * fr = new TH1D("fr5", ";Jet R;in-situ p_{a} (nominal)", 1, 0.15, 0.85);
    fr->SetMinimum(0.86); fr->SetMaximum(0.98); fr->Draw("axis");
    TLegend * l = legend(.5, .7, .93, .84);
    printf("in-situ p_a (nominal)\n");
    for (int m = 1; m < nMode; m++) {
      if (!pa.count(modes[m])) { cout << "WARNING: no p_a for " << modes[m] << endl; continue; }
      TGraphAsymmErrors * g = new TGraphAsymmErrors();
      printf("%-6s", modes[m]);
      for (int r = 0; r < ana::nJetR; r++) {
        g->SetPoint(r, ana::JetRs[r] + 0.012*(m - 2), pa[modes[m]][r]);
        g->SetPointError(r, 0, 0, eLo[modes[m]][r], eHi[modes[m]][r]);
        printf("  %.4f (%+.4f)", pa[modes[m]][r], pa[modes[m]][r] - pa["tmpl"][r]);
      }
      printf("\n");
      g->SetMarkerStyle(20); g->SetMarkerColor(cols[m]); g->SetLineColor(cols[m]);
      g->Draw("p same");
      l->AddEntry(g, labels[m], "p");
    }
    l->Draw();
    insitu_utility::drawSPhenixLabel({"p+p Run24 Data"}, {"Pythia8 #gamma+jet MC", "Purity-corrected, stat. errors",
        "Default: converged; linear: one scan"}, .18, .9, 14, c->GetWh());
    c->SaveAs(pdf.c_str());
  }

  // ---- p6: pt_ref bookkeeping ----
  {
    c->Clear();
    TLatex tx; tx.SetNDC(); tx.SetTextSize(0.03);
    tx.DrawLatex(0.08, 0.92, "Reference p_{T} for the re-smearing: R=0.4 region-A in-situ pairs, Photon5+10+20 (unweighted)");
    std::ifstream in(ana::path("claude_checks/smear_test/gammajet_inputs/bookkeeping.txt"));
    string mode; double y = 0.84;
    tx.DrawLatex(0.08, y, "mode        p_{T,calib} < 5 GeV (exact)     truth-matched (#DeltaR < 0.75R)     p_{T,calib} fallback");
    long n0, n1, n2;
    while (in >> mode >> n0 >> n1 >> n2) {
      y -= 0.06;
      double tot = n0 + n1 + n2;
      tx.DrawLatex(0.08, y, Form("%-8s     %ld (%.1f%%)          %ld (%.1f%%)          %ld (%.1f%%)", mode.c_str(), n0, 100*n0/tot, n1, 100*n1/tot, n2, 100*n2/tot));
    }
    c->SaveAs(pdf.c_str());
  }
  c->SaveAs((pdf + "]").c_str());
  printf("Wrote %s\n", pdf.c_str());
}
