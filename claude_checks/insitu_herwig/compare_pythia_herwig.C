#include "../../src/ana.h"
#include "../../src/insitu_utility.h"
#include <string>
#include <vector>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TStyle.h"
using namespace std;

R__LOAD_LIBRARY(libgammajet_unfold.so);

// Why does the in-situ p_a move by +0.6% (R=0.2) to +4.7% (R=0.8) when the MC reference is Herwig?
// Compares Pythia and Herwig gamma+jet MC (nominal systag) at each stage:
//   page 1: <x_J> Herwig/Pythia vs R, per pT bin, at truth level (response truth projection), at
//           reco level (the in-situ reference: insitu trees, region A, low-x_J floor), and 1/(p_a ratio).
//   pages:  per radius, normalized x_J shapes per pT bin: truth (top), reco + Data at nominal p_a (bottom).
//   then:   MC jet response <reco/truth jet pT> vs truth jet pT, per radius.
//   last:   MC photon response: <reco/truth photon pT> vs truth pT, and reco pT shapes per truth pT bin.
// MC reco jets are the smeared truth jets (unfolder.cc: jet_pt_smear_truth), so a generator difference
// at truth level carries straight into the reference.

const char * insitu_input_dir = ana::path("insitu/inputs");
const char * hists_dir = ana::path("hists");
const char * insitu_output_dir = ana::path("insitu/output");
const int nPt = ana::nPtBinsUsed;
const int nGlob = ana::nUnfoldXjBins + 2; // findUnfoldBin's per-pT-bin stride

// cross-section weights (drawer.h's scalemap)
map<string, map<int,double>> xsec = {
  {"pythia", {{5,146359.3},{10,6944.675},{20,130.4461}}},
  {"herwig", {{5,6.48487e+07},{10,3.62808e+02},{20,5.34010e+01}}},
};

// Sum of hname over Photon5/10/20 unfolding outputs, cross-section weighted.
TH1 * sumUnfold(const string & sim, const char * hname) {
  TH1 * sum = nullptr;
  for (int s : {5, 10, 20}) {
    TFile * f = TFile::Open(Form("%s/Photon%d_%s_nominal_unfolding.root", hists_dir, s, sim.c_str()), "READ");
    if (!f || f->IsZombie()) continue;
    TH1 * h = (TH1*)f->Get(hname);
    if (!h) { f->Close(); continue; }
    h = (TH1*)h->Clone(Form("%s_%s_%d", hname, sim.c_str(), s)); h->SetDirectory(0);
    h->Scale(xsec[sim][s]);
    if (!sum) { sum = (TH1*)h->Clone(Form("%s_%s", hname, sim.c_str())); sum->SetDirectory(0); }
    else sum->Add(h);
    f->Close();
  }
  return sum;
}

// x_J shape of reported pT bin ipt out of a flattened (findUnfoldBin) histogram, normalized to unit area.
TH1D * shapeFromGlobal(TH1 * g, int ipt, const char * name) {
  TH1D * h = new TH1D(name, ";x_{J#gamma};(1/N) dN/dx_{J#gamma}", ana::nUnfoldXjBins, ana::unfoldXjBins);
  int iptAll = ipt + ana::firstUsedPtBin;
  for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
    int b = g->FindBin(iptAll*nGlob + ixj + 1);
    h->SetBinContent(ixj+1, g->GetBinContent(b)); h->SetBinError(ixj+1, g->GetBinError(b));
  }
  if (h->Integral() > 0) h->Scale(1.0/h->Integral(), "width");
  return h;
}

// Reco x_J shape (region A, low-x_J floor) from in-situ trees; Data jets divided by pa.
TH1D * shapeFromInsitu(const string & sim, int ir, int ipt, float floorXj, double pa, const char * name) {
  TH1D * h = new TH1D(name, ";x_{J#gamma};(1/N) dN/dx_{J#gamma}", ana::nUnfoldXjBins, ana::unfoldXjBins);
  h->Sumw2();
  vector<pair<string,double>> files;
  if (sim == "Data") files = {{insitu_utility::insituFilename(insitu_input_dir, "Data", "", "nominal"), 1.0}};
  else for (int s : {5, 10, 20}) files.push_back({insitu_utility::insituFilename(insitu_input_dir, Form("Photon%d", s), sim.c_str(), "nominal"), xsec[sim][s]});
  for (auto & fs : files) {
    TFile * f = TFile::Open(fs.first.c_str(), "READ");
    if (!f || f->IsZombie()) continue;
    TTree * t = (TTree*)f->Get("insitutree");
    Float_t pho_pt, jet_pt, w; Int_t abcd, evIr;
    t->SetBranchAddress("pho_pt", &pho_pt); t->SetBranchAddress("jet_pt", &jet_pt);
    t->SetBranchAddress("abcd", &abcd); t->SetBranchAddress("weight", &w); t->SetBranchAddress("ir", &evIr);
    for (Long64_t e = 0; e < t->GetEntries(); e++) {
      t->GetEntry(e);
      if (abcd != 0 || evIr != ir) continue;
      if (ana::findPtBin(pho_pt) != ipt + ana::firstUsedPtBin) continue;
      double x = jet_pt/pa/pho_pt;
      if (x < floorXj) continue;
      h->Fill(x, fs.second*w);
    }
    f->Close();
  }
  if (h->Integral() > 0) h->Scale(1.0/h->Integral(), "width");
  return h;
}

// Mean and its error of a unit-area density histogram.
void meanOf(TH1D * h, double & m, double & e) {
  double sw = 0, swx = 0, swx2 = 0, se2 = 0;
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    double w = h->GetBinContent(i)*h->GetBinWidth(i), x = h->GetBinCenter(i);
    sw += w; swx += w*x;
  }
  m = sw > 0 ? swx/sw : 0;
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    double ew = h->GetBinError(i)*h->GetBinWidth(i), x = h->GetBinCenter(i);
    se2 += ew*ew*(x - m)*(x - m);
  }
  e = sw > 0 ? sqrt(se2)/sw : 0;
}

double readPa(const string & tag, int ir) {
  TFile * f = TFile::Open(Form("%s/grid_insitu_%s.root", insitu_output_dir, tag.c_str()), "READ");
  if (!f || f->IsZombie()) return -1;
  TTree * t = (TTree*)f->Get(Form("%s/results", ana::rnames[ir]));
  float pa = -1; t->SetBranchAddress("pa_puritycorrected", &pa); t->GetEntry(0); f->Close();
  return pa;
}

void style(TH1 * h, int col, int mk, int ls = 1) {
  h->SetLineColor(col); h->SetMarkerColor(col); h->SetMarkerStyle(mk); h->SetLineStyle(ls); h->SetLineWidth(2); h->SetMarkerSize(0.8);
}

void compare_pythia_herwig() {
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  TH1::AddDirectory(kFALSE);
  const int nR = ana::nJetR;
  const int colP = kBlue+1, colH = kRed+1;
  string pdf = "compare_pythia_herwig.pdf";

  // ---- compute everything ----
  // [ir][ipt]
  vector<vector<TH1D*>> tP(nR, vector<TH1D*>(nPt)), tH = tP, rP = tP, rH = tP, rD = tP;
  vector<vector<double>> mtP(nR, vector<double>(nPt)), mtH = mtP, mrP = mtP, mrH = mtP, mrD = mtP;
  vector<vector<double>> etP(nR, vector<double>(nPt)), etH = etP, erP = etP, erH = etP, erD = etP;
  vector<double> paNom(nR), paHer(nR);
  for (int ir = 0; ir < nR; ir++) {
    paNom[ir] = readPa("nominal", ir); paHer[ir] = readPa("herwig", ir);
    TH1 * gtP = sumUnfold("pythia", Form("htruthxj%d", ir));
    TH1 * gtH = sumUnfold("herwig", Form("htruthxj%d", ir));
    for (int ipt = 0; ipt < nPt; ipt++) {
      float fl = insitu_utility::lowXjFloor(ir, ana::ptBinsUsed[ipt]);
      tP[ir][ipt] = shapeFromGlobal(gtP, ipt, Form("tP_%d_%d", ir, ipt));
      tH[ir][ipt] = shapeFromGlobal(gtH, ipt, Form("tH_%d_%d", ir, ipt));
      rP[ir][ipt] = shapeFromInsitu("pythia", ir, ipt, fl, 1.0, Form("rP_%d_%d", ir, ipt));
      rH[ir][ipt] = shapeFromInsitu("herwig", ir, ipt, fl, 1.0, Form("rH_%d_%d", ir, ipt));
      rD[ir][ipt] = shapeFromInsitu("Data", ir, ipt, fl, paNom[ir], Form("rD_%d_%d", ir, ipt));
      meanOf(tP[ir][ipt], mtP[ir][ipt], etP[ir][ipt]); meanOf(tH[ir][ipt], mtH[ir][ipt], etH[ir][ipt]);
      meanOf(rP[ir][ipt], mrP[ir][ipt], erP[ir][ipt]); meanOf(rH[ir][ipt], mrH[ir][ipt], erH[ir][ipt]);
      meanOf(rD[ir][ipt], mrD[ir][ipt], erD[ir][ipt]);
    }
  }

  printf("\n<x_J> (binned means of the shapes; truth: response truth projection, reco: in-situ trees region A + floor)\n");
  printf("%-5s %-8s %8s %8s %8s | %8s %8s %8s | %8s | p_a H/P\n", "R", "pT", "trP", "trH", "H/P", "recoP", "recoH", "H/P", "Data");
  for (int ir = 0; ir < nR; ir++) for (int ipt = 0; ipt < nPt; ipt++)
    printf("%-5.1f %2.0f-%-5.0f %8.4f %8.4f %8.4f | %8.4f %8.4f %8.4f | %8.4f | %.4f\n", ana::JetRs[ir], ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1],
        mtP[ir][ipt], mtH[ir][ipt], mtH[ir][ipt]/mtP[ir][ipt], mrP[ir][ipt], mrH[ir][ipt], mrH[ir][ipt]/mrP[ir][ipt], mrD[ir][ipt], paHer[ir]/paNom[ir]);

  TCanvas * c = new TCanvas("c", "", 1500, 900);
  c->SaveAs((pdf + "[").c_str());

  // ---- page 1: mean ratios vs R ----
  c->Divide(3, 1, 0.002, 0.002);
  for (int ipt = 0; ipt < nPt; ipt++) {
    c->cd(ipt+1);
    gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.11); gPad->SetTicks(1,1);
    TH1D * fr = new TH1D(Form("frm%d", ipt), ";jet radius R;Herwig / Pythia", 1, 0.1, 0.9);
    fr->SetMinimum(0.90); fr->SetMaximum(1.04);
    fr->GetYaxis()->SetTitleOffset(1.5); fr->Draw("axis");
    TLine * one = new TLine(0.1, 1, 0.9, 1); one->SetLineStyle(2); one->Draw();
    TGraphErrors * gt = new TGraphErrors(), * gr = new TGraphErrors(), * gp = new TGraphErrors();
    for (int ir = 0; ir < nR; ir++) {
      double R = ana::JetRs[ir];
      double q = mtH[ir][ipt]/mtP[ir][ipt];
      gt->SetPoint(ir, R - 0.01, q); gt->SetPointError(ir, 0, q*hypot(etH[ir][ipt]/mtH[ir][ipt], etP[ir][ipt]/mtP[ir][ipt]));
      q = mrH[ir][ipt]/mrP[ir][ipt];
      gr->SetPoint(ir, R + 0.01, q); gr->SetPointError(ir, 0, q*hypot(erH[ir][ipt]/mrH[ir][ipt], erP[ir][ipt]/mrP[ir][ipt]));
      gp->SetPoint(ir, R, paNom[ir]/paHer[ir]);
    }
    gt->SetMarkerStyle(24); gt->SetMarkerColor(kGreen+2); gt->SetLineColor(kGreen+2); gt->Draw("p same");
    gr->SetMarkerStyle(20); gr->SetMarkerColor(kRed+1); gr->SetLineColor(kRed+1); gr->Draw("p same");
    gp->SetMarkerStyle(33); gp->SetMarkerSize(1.6); gp->SetMarkerColor(kBlack); gp->Draw("p same");
    insitu_utility::drawSPhenixLabel({"Pythia8 vs Herwig7 #gamma+jet MC"},
        {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1])}, .22, .9, 18, gPad->GetWh());
    if (ipt == 0) {
      TLegend * l = new TLegend(0.22, 0.14, 0.95, 0.34); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.035);
      l->AddEntry(gt, "truth #LTx_{J}#GT ratio", "p");
      l->AddEntry(gr, "reco #LTx_{J}#GT ratio (in-situ ref., region A)", "p");
      l->AddEntry(gp, "p_{a}^{Pythia} / p_{a}^{Herwig} (all p_{T})", "p");
      l->Draw();
    }
  }
  c->SaveAs(pdf.c_str());

  // ---- per-radius shape pages ----
  for (int ir = 0; ir < nR; ir++) {
    c->Clear();
    c->Divide(3, 2, 0.002, 0.002);
    for (int row = 0; row < 2; row++) for (int ipt = 0; ipt < nPt; ipt++) {
      c->cd(row*3 + ipt + 1);
      gPad->SetLeftMargin(0.14); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.13); gPad->SetTicks(1,1);
      TH1D * a = row == 0 ? tP[ir][ipt] : rP[ir][ipt];
      TH1D * b = row == 0 ? tH[ir][ipt] : rH[ir][ipt];
      style(a, colP, 20); style(b, colH, 21);
      double ymax = std::max(a->GetMaximum(), b->GetMaximum());
      if (row == 1) { style(rD[ir][ipt], kBlack, 24); ymax = std::max(ymax, rD[ir][ipt]->GetMaximum()); }
      a->SetMaximum(1.5*ymax); a->SetMinimum(0);
      a->GetXaxis()->SetTitleSize(0.05); a->GetYaxis()->SetTitleSize(0.05); a->GetYaxis()->SetTitleOffset(1.25);
      a->Draw("hist e"); b->Draw("hist e same");
      if (row == 1) rD[ir][ipt]->Draw("p same");
      double mA = row == 0 ? mtP[ir][ipt] : mrP[ir][ipt], mB = row == 0 ? mtH[ir][ipt] : mrH[ir][ipt];
      TLegend * l = new TLegend(0.52, 0.62, 0.97, 0.80); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.042);
      l->AddEntry(a, Form("Pythia #LTx_{J}#GT=%.3f", mA), "l");
      l->AddEntry(b, Form("Herwig #LTx_{J}#GT=%.3f", mB), "l");
      if (row == 1) l->AddEntry(rD[ir][ipt], Form("Data/p_{a}^{nom} #LTx_{J}#GT=%.3f", mrD[ir][ipt]), "p");
      l->Draw();
      insitu_utility::drawSPhenixLabel({row == 0 ? "truth level" : "reco level, region A (in-situ ref.)"},
          {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1]), Form("Jet R=%.1f", ana::JetRs[ir])},
          .17, .9, 15, gPad->GetWh()*gPad->GetHNDC());
    }
    c->SaveAs(pdf.c_str());
  }

  // ---- jet response ----
  c->Clear();
  c->Divide(4, 2, 0.002, 0.002);
  for (int ir = 0; ir < nR; ir++) {
    c->cd(ir+1);
    gPad->SetLeftMargin(0.16); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.13); gPad->SetTicks(1,1);
    TH2D * hP = (TH2D*)sumUnfold("pythia", Form("hjetresponse%d", ir));
    TH2D * hH = (TH2D*)sumUnfold("herwig", Form("hjetresponse%d", ir));
    if (!hP || !hH) continue;
    // Hresponse: x = reco, y = truth. <reco>/truth-bin-center per truth bin.
    TH1D * qP = hP->ProjectionY(Form("qP%d", ir)); qP->Reset(); TH1D * qH = (TH1D*)qP->Clone(Form("qH%d", ir));
    for (auto pr : {make_pair(hP, qP), make_pair(hH, qH)}) {
      for (int j = 1; j <= pr.first->GetNbinsY(); j++) {
        TH1D * px = pr.first->ProjectionX("_px", j, j);
        if (px->GetEntries() < 5) continue;
        double yc = pr.first->GetYaxis()->GetBinCenter(j);
        pr.second->SetBinContent(j, px->GetMean()/yc); pr.second->SetBinError(j, px->GetMeanError()/yc);
        delete px;
      }
    }
    style(qP, colP, 20); style(qH, colH, 21);
    qP->SetTitle(";truth jet p_{T} [GeV];#LTreco p_{T}#GT / truth p_{T}"); qP->GetXaxis()->SetRangeUser(5, 50);
    qP->SetMinimum(0.8); qP->SetMaximum(1.3); qP->GetYaxis()->SetTitleOffset(1.4);
    qP->GetXaxis()->SetTitleSize(0.05); qP->GetYaxis()->SetTitleSize(0.05);
    qP->Draw("e"); qH->Draw("e same");
    TLine * one = new TLine(qP->GetXaxis()->GetBinLowEdge(qP->GetXaxis()->GetFirst()), 1, 50, 1); one->SetLineStyle(2); one->Draw();
    insitu_utility::drawSPhenixLabel({"MC jet response (matched)"}, {Form("Jet R=%.1f", ana::JetRs[ir])}, .2, .9, 15, gPad->GetWh()*gPad->GetHNDC());
    if (ir == 0) {
      TLegend * l = new TLegend(0.2, 0.18, 0.9, 0.32); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.05);
      l->AddEntry(qP, "Pythia", "p"); l->AddEntry(qH, "Herwig", "p"); l->Draw();
    }
  }
  c->SaveAs(pdf.c_str());

  // ---- photon response (matched; photon-only, so one radius, R=0.4) ----
  c->Clear();
  c->Divide(4, 1, 0.002, 0.002);
  {
    TH2D * hP = (TH2D*)sumUnfold("pythia", "hphoresponse2");
    TH2D * hH = (TH2D*)sumUnfold("herwig", "hphoresponse2");
    c->cd(1);
    gPad->SetLeftMargin(0.17); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.12); gPad->SetTicks(1,1);
    TH1D * qP = hP->ProjectionY("pqP"); qP->Reset(); TH1D * qH = (TH1D*)qP->Clone("pqH");
    for (auto pr : {make_pair(hP, qP), make_pair(hH, qH)}) {
      for (int j = 1; j <= pr.first->GetNbinsY(); j++) {
        TH1D * px = pr.first->ProjectionX("_ppx", j, j);
        if (px->GetEntries() < 5) { delete px; continue; }
        double yc = pr.first->GetYaxis()->GetBinCenter(j);
        pr.second->SetBinContent(j, px->GetMean()/yc); pr.second->SetBinError(j, px->GetMeanError()/yc);
        delete px;
      }
    }
    style(qP, colP, 20); style(qH, colH, 21);
    qP->SetTitle(";truth photon p_{T} [GeV];#LTreco p_{T}#GT / truth p_{T}"); qP->GetXaxis()->SetRangeUser(13, 40);
    qP->SetMinimum(0.85); qP->SetMaximum(1.1); qP->GetYaxis()->SetTitleOffset(1.6);
    qP->Draw("e"); qH->Draw("e same");
    TLine * one = new TLine(13, 1, 40, 1); one->SetLineStyle(2); one->Draw();
    insitu_utility::drawSPhenixLabel({"MC photon response (matched)"}, {"p_{T}^{#gamma} bins as analysis"}, .22, .9, 16, gPad->GetWh());
    TLegend * l = new TLegend(0.22, 0.16, 0.9, 0.28); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.045);
    l->AddEntry(qP, "Pythia", "p"); l->AddEntry(qH, "Herwig", "p"); l->Draw();
    // reco photon pT / truth bin centre shape per reported truth pT bin
    for (int ipt = 0; ipt < nPt; ipt++) {
      c->cd(ipt+2);
      gPad->SetLeftMargin(0.15); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.12); gPad->SetTicks(1,1);
      double lo = ana::ptBinsUsed[ipt], hi = ana::ptBinsUsed[ipt+1];
      int j1 = hP->GetYaxis()->FindBin(lo + 1e-3), j2 = hP->GetYaxis()->FindBin(hi - 1e-3);
      TH1D * a = hP->ProjectionX(Form("prP%d", ipt), j1, j2), * b = hH->ProjectionX(Form("prH%d", ipt), j1, j2);
      a->Scale(1.0/a->Integral()); b->Scale(1.0/b->Integral());
      style(a, colP, 20); style(b, colH, 21);
      a->SetTitle(";reco photon p_{T} [GeV];normalized"); a->GetXaxis()->SetRangeUser(lo - 8, hi + 8);
      a->SetMaximum(1.5*std::max(a->GetMaximum(), b->GetMaximum())); a->SetMinimum(0); a->GetYaxis()->SetTitleOffset(1.5);
      a->Draw("hist e"); b->Draw("hist e same");
      TLegend * l2 = new TLegend(0.5, 0.66, 0.97, 0.80); l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextSize(0.042);
      l2->AddEntry(a, Form("Pythia #LTreco#GT=%.2f", a->GetMean()), "l"); l2->AddEntry(b, Form("Herwig #LTreco#GT=%.2f", b->GetMean()), "l"); l2->Draw();
      insitu_utility::drawSPhenixLabel({"truth photon in"}, {Form("%.0f GeV < p_{T}^{#gamma,truth} < %.0f GeV", lo, hi)}, .19, .9, 16, gPad->GetWh());
    }
  }
  c->SaveAs(pdf.c_str());
  c->SaveAs((pdf + "]").c_str());
  cout << "Wrote " << pdf << endl;
}
