#include "../../src/ana.h"
#include "../../src/drawer.h"
#include "../../src/insitu_utility.h"
#include "../../src/unfold_utility.h"
// Explicit load; run interpreted, never with ACLiC "+" (sibling libgammajet.so collision).
R__LOAD_LIBRARY(libgammajet_unfold.so);

// PPG18 review round 2, issue 3: data counterpart of draw_jet12_purity_corrected_xj.C.
// Run that macro first - this one reads its output (Jet12_long ABCD purities, region-A
// uncorrected / corrected / truth-tagged xJ) and puts Data next to it, R = 0.2/0.3/0.4.
//
// Data side is built with exactly the in-situ scan's ingredients (grid_insitu.C):
//   insitu_utility::cacheDataEvents(Data_nominal_insitu.root, A or C, ir) at the RAW
//   jet scale (pa = 1, i.e. jet_pt_calib with no in-situ factor), the xJ floor
//   insitu_utility::lowXjFloor, the committed data purities ana::getPurity/getPurityC
//   (hists/purity_nominal.root, per R), and insitu_utility::purityCorrectByPtBin.
// The in-situ reference is also drawn: Pythia8 gamma+jet (Photon5+10+20, same cross-
// section weights as grid_insitu.C) region A via insitu_utility::buildMCXjByPtBin.
//
// Pages:
//   1       purity vs photon pT per R: data P_A, P_C; Jet12 ABCD P_A, P_C; Jet12 true P_A
//   2..4    per R: top = Data uncorrected / purity corrected / gamma+jet MC reference;
//           bottom = corrected/uncorrected shape ratio for Data and for Jet12 (how much
//           the subtraction moves each bin in each sample)
//   log     means (binned) and the implied in-situ scale <Data>/<MC> before and after
//           the correction - what the mean-xJ scan is fitting.

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

void draw_data_vs_jet12_purity_xj()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();

  const int irs[3] = {0, 1, 2};
  const int nR = 3, nPt = ana::nPtBinsUsed;
  const char * inDir = ana::path("insitu/inputs");
  const char * systag = "nominal";
  map<int,double> photon_scale = {{5,146359.3},{10,6944.675},{20,130.4461}}; // = grid_insitu.C
  string outdir = ana::path("claude_checks/insitu_scan_edge/pdfs");
  TFile * fj = TFile::Open((outdir + "/draw_jet12_purity_corrected_xj.root").c_str(), "read");
  if (!fj || fj->IsZombie()) { cout << "Run draw_jet12_purity_corrected_xj.C first." << endl; return; }

  TH1D * dA[nR][nPt], * dCorr[nR][nPt], * mcRef[nR][nPt];
  TH1D * jA[nR][nPt], * jCorr[nR][nPt], * jSig[nR][nPt];
  float pA[nR][nPt], pAlo[nR][nPt], pAhi[nR][nPt], pC[nR][nPt], pClo[nR][nPt], pChi[nR][nPt];
  double jpA[nR][nPt], jpAlo[nR][nPt], jpAhi[nR][nPt], jpC[nR][nPt], jpClo[nR][nPt], jpChi[nR][nPt], jpTrue[nR][nPt];

  for (int r = 0; r < nR; r++) {
    int ir = irs[r];
    float lowXj[nPt];
    for (int p = 0; p < nPt; p++) lowXj[p] = insitu_utility::lowXjFloor(ir, ana::ptBinsUsed[p]);
    string dataFile = insitu_utility::insituFilename(inDir, "Data", "", systag);
    vector<DataEvent> evA = insitu_utility::cacheDataEvents(dataFile.c_str(), 0, ir);
    vector<DataEvent> evC = insitu_utility::cacheDataEvents(dataFile.c_str(), 2, ir);
    vector<TH1D*> hA = insitu_utility::buildXjByPtBin(evA, 1.0, nPt, Form("dataA_r%d", ir), lowXj);
    vector<TH1D*> hC = insitu_utility::buildXjByPtBin(evC, 1.0, nPt, Form("dataC_r%d", ir), lowXj);
    for (int p = 0; p < nPt; p++) {
      float lo = ana::ptBinsUsed[p], hi = ana::ptBinsUsed[p+1];
      pA[r][p] = ana::getPurity(lo, hi, systag, ir);
      pAlo[r][p] = ana::getPurityErrorLow(lo, hi, systag, ir);
      pAhi[r][p] = ana::getPurityErrorHigh(lo, hi, systag, ir);
      pC[r][p] = ana::getPurityC(lo, hi, systag, ir);
      pClo[r][p] = ana::getPurityCErrorLow(lo, hi, systag, ir);
      pChi[r][p] = ana::getPurityCErrorHigh(lo, hi, systag, ir);
    }
    vector<TH1D*> hCorr = insitu_utility::purityCorrectByPtBin(hA, hC, nPt, pA[r], pAlo[r], pAhi[r], pC[r], pClo[r], pChi[r],
                                                               Form("dataCorr_r%d", ir));
    vector<TH1D*> hMC = insitu_utility::buildMCXjByPtBin({
        {insitu_utility::insituFilename(inDir, "Photon5",  "pythia", systag), photon_scale[5]},
        {insitu_utility::insituFilename(inDir, "Photon10", "pythia", systag), photon_scale[10]},
        {insitu_utility::insituFilename(inDir, "Photon20", "pythia", systag), photon_scale[20]},
      }, 0, ir, Form("mcRefA_r%d", ir), lowXj);

    TGraphAsymmErrors * gP  = (TGraphAsymmErrors*)fj->Get(Form("%s/combined", ana::rnames[ir]));
    TGraphAsymmErrors * gPC = (TGraphAsymmErrors*)fj->Get(Form("%s/combined_C", ana::rnames[ir]));
    for (int p = 0; p < nPt; p++) {
      dA[r][p] = hA[p]; dCorr[r][p] = hCorr[p]; mcRef[r][p] = hMC[p];
      jA[r][p]    = (TH1D*)fj->Get(Form("hA_r%d_pt%d", ir, p));
      jCorr[r][p] = (TH1D*)fj->Get(Form("hCorr_r%d_pt%d", ir, p));
      jSig[r][p]  = (TH1D*)fj->Get(Form("hSig_r%d_pt%d", ir, p));
      int ib = ana::firstUsedPtBin + p;
      double x;
      gP->GetPoint(ib, x, jpA[r][p]);  jpAlo[r][p] = gP->GetErrorYlow(ib);  jpAhi[r][p] = gP->GetErrorYhigh(ib);
      gPC->GetPoint(ib, x, jpC[r][p]); jpClo[r][p] = gPC->GetErrorYlow(ib); jpChi[r][p] = gPC->GetErrorYhigh(ib);
      jpTrue[r][p] = jSig[r][p]->Integral()/std::max(1e-9, jA[r][p]->Integral());
    }
  }

  // ---------------- numbers ----------------
  printf("\nPurity (bootstrap median, -lo/+hi):\n");
  printf("%-4s %-7s | %-22s %-22s | %-22s %-22s %8s\n", "R", "pT", "Data P_A", "Data P_C", "Jet12 P_A (ABCD)", "Jet12 P_C (ABCD)", "Jet12 true");
  for (int r = 0; r < nR; r++) for (int p = 0; p < nPt; p++)
    printf("%-4.1f %2.0f-%-4.0f | %.3f -%.3f +%.3f     %.3f -%.3f +%.3f     | %.3f -%.3f +%.3f     %.3f -%.3f +%.3f     %8.3f\n",
           ana::JetRs[irs[r]], ana::ptBinsUsed[p], ana::ptBinsUsed[p+1],
           pA[r][p], pAlo[r][p], pAhi[r][p], pC[r][p], pClo[r][p], pChi[r][p],
           jpA[r][p], jpAlo[r][p], jpAhi[r][p], jpC[r][p], jpClo[r][p], jpChi[r][p], jpTrue[r][p]);

  printf("\nBinned <xJ> (errors x1000). Data at raw jet scale (pa = 1). 'shift' = corrected - uncorrected.\n");
  printf("%-4s %-7s | %9s %9s %7s | %9s %9s %7s %9s | %9s | %6s %6s\n", "R", "pT", "Data A", "Data corr", "shift",
         "J12 A", "J12 corr", "shift", "J12 tag", "gj MC A", "raw/MC", "corr/MC");
  for (int r = 0; r < nR; r++) for (int p = 0; p < nPt; p++) {
    double e1, e2, e3, e4, e5, e6;
    double mdA = binMean(dA[r][p], &e1), mdC = binMean(dCorr[r][p], &e2);
    double mjA = binMean(jA[r][p], &e3), mjC = binMean(jCorr[r][p], &e4), mjS = binMean(jSig[r][p], &e5);
    double mMC = binMean(mcRef[r][p], &e6);
    printf("%-4.1f %2.0f-%-4.0f | %.3f(%2.0f) %.3f(%2.0f) %+.3f | %.3f(%2.0f) %.3f(%2.0f) %+.3f %.3f(%2.0f) | %.3f(%2.0f) | %6.3f %6.3f\n",
           ana::JetRs[irs[r]], ana::ptBinsUsed[p], ana::ptBinsUsed[p+1],
           mdA, 1000*e1, mdC, 1000*e2, mdC-mdA, mjA, 1000*e3, mjC, 1000*e4, mjC-mjA, mjS, 1000*e5, mMC, 1000*e6,
           mMC > 0 ? mdA/mMC : 0, mMC > 0 ? mdC/mMC : 0);
  }
  printf("(raw/MC and corr/MC ~ the single scale factor the mean-xJ scan fits per bin; the scan floor is %.3f)\n", insitu_utility::scanLow);

  // ---------------- drawing ----------------
  drawer d("pythia", "nominal");
  string pdf = outdir + "/draw_data_vs_jet12_purity_xj.pdf";
  TCanvas * c = new TCanvas("c", "", 1500, 800);
  c->SaveAs((pdf+"[").c_str());

  // page 1: purity
  c->Clear(); c->Divide(3, 1);
  for (int r = 0; r < nR; r++) {
    c->cd(r+1); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05);
    TH1D * fr = new TH1D(Form("frP%d", r), ";p_{T}^{#gamma} [GeV];purity", 1, ana::ptBinsUsed[0], ana::ptBinsUsed[nPt]);
    fr->SetMinimum(0); fr->SetMaximum(1.6); fr->GetYaxis()->SetTitleSize(0.05); fr->GetXaxis()->SetTitleSize(0.05); fr->Draw();
    auto mkG = [&](int col, int mk, double off) { TGraphAsymmErrors * g = new TGraphAsymmErrors(); g->SetMarkerColor(col); g->SetLineColor(col); g->SetMarkerStyle(mk); g->SetMarkerSize(1.2); return g; };
    TGraphAsymmErrors * gDA = mkG(kBlack, 20, 0), * gDC = mkG(kBlack, 24, 0), * gJA = mkG(kRed+1, 21, 0), * gJC = mkG(kRed+1, 25, 0), * gJT = mkG(kGreen+2, 29, 0);
    for (int p = 0; p < nPt; p++) {
      double x = 0.5*(ana::ptBinsUsed[p]+ana::ptBinsUsed[p+1]);
      gDA->SetPoint(p, x-0.6, pA[r][p]); gDA->SetPointError(p, 0, 0, pAlo[r][p], pAhi[r][p]);
      gDC->SetPoint(p, x-0.6, pC[r][p]); gDC->SetPointError(p, 0, 0, pClo[r][p], pChi[r][p]);
      gJA->SetPoint(p, x+0.6, jpA[r][p]); gJA->SetPointError(p, 0, 0, jpAlo[r][p], jpAhi[r][p]);
      gJC->SetPoint(p, x+0.6, jpC[r][p]); gJC->SetPointError(p, 0, 0, jpClo[r][p], jpChi[r][p]);
      gJT->SetPoint(p, x+1.2, jpTrue[r][p]);
    }
    for (auto g : {gDA, gDC, gJA, gJC, gJT}) g->Draw("p same");
    TLegend * l = new TLegend(0.18, 0.56, 0.62, 0.78); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.035);
    l->AddEntry(gDA, "Data P_{A}", "p"); l->AddEntry(gDC, "Data P_{C}", "p");
    l->AddEntry(gJA, "Jet12 P_{A} (ABCD)", "p"); l->AddEntry(gJC, "Jet12 P_{C} (ABCD)", "p");
    l->AddEntry(gJT, "Jet12 true P_{A} (tagged/A)", "p");
    l->Draw();
    d.drawAll({"Data vs Pythia8 QCD (Jet12_long)"}, {Form("Jet R=%.1f", ana::JetRs[irs[r]]), "paired, same ABCD method"}, .18, .9, 13, gPad->GetWh()*0.8);
  }
  c->SaveAs(pdf.c_str());

  // pages 2-4: xJ
  for (int r = 0; r < nR; r++) {
    int ir = irs[r];
    c->Clear();
    for (int p = 0; p < nPt; p++) {
      c->cd();
      TPad * top = new TPad(Form("dtop%d_%d",r,p), "", p/3.0, 0.35, (p+1)/3.0, 1.0);
      TPad * bot = new TPad(Form("dbot%d_%d",r,p), "", p/3.0, 0.0, (p+1)/3.0, 0.35);
      top->SetLeftMargin(0.16); top->SetRightMargin(0.03); top->SetBottomMargin(0.02); top->SetTopMargin(0.05);
      bot->SetLeftMargin(0.16); bot->SetRightMargin(0.03); bot->SetTopMargin(0.02); bot->SetBottomMargin(0.3);
      top->Draw(); bot->Draw();
      top->cd();
      TH1D * sA = shape(dA[r][p], Form("sdA_%d_%d", r, p));
      TH1D * sCo = shape(dCorr[r][p], Form("sdC_%d_%d", r, p));
      TH1D * sMC = shape(mcRef[r][p], Form("sMC_%d_%d", r, p));
      double ymax = std::max({sA->GetMaximum(), sCo->GetMaximum(), sMC->GetMaximum()});
      sMC->SetLineColor(kGreen+2); sMC->SetFillColor(kGreen-9); sMC->SetFillStyle(1001); sMC->SetMarkerSize(0);
      sA->SetLineColor(kGray+2); sA->SetMarkerColor(kGray+2); sA->SetMarkerStyle(24);
      sCo->SetLineColor(kBlack); sCo->SetMarkerColor(kBlack); sCo->SetMarkerStyle(20);
      sMC->GetXaxis()->SetRangeUser(0, 2); sMC->GetXaxis()->SetLabelSize(0);
      sMC->GetYaxis()->SetTitle("(1/N) dN/dx_{J#gamma}");
      sMC->GetYaxis()->SetTitleSize(0.055); sMC->GetYaxis()->SetLabelSize(0.045); sMC->GetYaxis()->SetTitleOffset(1.3);
      sMC->SetMinimum(0); sMC->SetMaximum(1.7*ymax);
      sMC->Draw("e2"); sA->Draw("e same"); sCo->Draw("e same");
      TLegend * l = new TLegend(0.6, 0.55, 0.97, 0.78); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.042);
      l->AddEntry(sA, "Data, uncorr.", "lp");
      l->AddEntry(sCo, "Data, corr.", "lp");
      l->AddEntry(sMC, "#gamma+jet MC", "f");
      l->Draw();
      d.drawAll({"p+p Run24 Data"}, {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ana::ptBinsUsed[p], ana::ptBinsUsed[p+1]),
                Form("Jet R=%.1f, raw jet scale", ana::JetRs[ir]), Form("P_{A} = %.2f, P_{C} = %.2f", pA[r][p], pC[r][p])}, .2, .88, 13, gPad->GetWh()*0.8);
      bot->cd();
      TH1D * sjA = shape(jA[r][p], Form("sjA_%d_%d", r, p));
      TH1D * sjC = shape(jCorr[r][p], Form("sjC_%d_%d", r, p));
      TH1D * rD = (TH1D*)sCo->Clone(Form("rD_%d_%d", r, p)); rD->Divide(sA);
      TH1D * rJ = (TH1D*)sjC->Clone(Form("rJ_%d_%d", r, p)); rJ->Divide(sjA);
      sanitize(rD); sanitize(rJ);
      // same events on both sides of each ratio - drawn without errors, it shows only
      // how far the correction moves each bin
      for (TH1D * h : {rD, rJ}) for (int b = 0; b <= h->GetNbinsX()+1; b++) h->SetBinError(b, 0);
      rD->SetMarkerStyle(20); rD->SetMarkerColor(kBlack); rD->SetLineColor(kBlack);
      rJ->SetMarkerStyle(21); rJ->SetMarkerColor(kRed+1); rJ->SetLineColor(kRed+1);
      rD->SetMinimum(0.4); rD->SetMaximum(1.6);
      rD->GetYaxis()->SetTitle("corr. / uncorr."); rD->GetYaxis()->SetNdivisions(505);
      rD->GetYaxis()->SetTitleSize(0.1); rD->GetYaxis()->SetLabelSize(0.08); rD->GetYaxis()->SetTitleOffset(0.7);
      rD->GetXaxis()->SetTitleSize(0.11); rD->GetXaxis()->SetLabelSize(0.09);
      rD->GetXaxis()->SetRangeUser(0, 2);
      rD->Draw("p"); rJ->Draw("p same");
      TLine * one = new TLine(0, 1, 2, 1); one->SetLineStyle(2); one->Draw();
      TLegend * lb = new TLegend(0.2, 0.78, 0.9, 0.95); lb->SetNColumns(2); lb->SetBorderSize(0); lb->SetFillStyle(0); lb->SetTextSize(0.08);
      lb->AddEntry(rD, "Data", "p"); lb->AddEntry(rJ, "Jet12 MC", "p"); lb->Draw();
    }
    c->SaveAs(pdf.c_str());
  }
  c->SaveAs((pdf+"]").c_str());

  TFile * fout = TFile::Open((outdir + "/draw_data_vs_jet12_purity_xj.root").c_str(), "recreate");
  for (int r = 0; r < nR; r++) for (int p = 0; p < nPt; p++) { dA[r][p]->Write(); dCorr[r][p]->Write(); mcRef[r][p]->Write(); }
  fout->Close();
}
