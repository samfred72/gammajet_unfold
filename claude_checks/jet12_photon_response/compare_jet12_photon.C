#include "../../src/ana.h"
#include "../../src/insitu_utility.h"
#include <string>
#include <vector>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
using namespace std;

R__LOAD_LIBRARY(libgammajet_unfold.so);

// Is the Jet12_long sample's detector response consistent with the Pythia Photon5/10/20 samples?
// (The Herwig production had a ~2.5% higher photon response; claude_checks/insitu_herwig.) Nominal
// systag, truth-matched region-A pairs (the response-matrix selection).
//   page 1: photon response <reco/truth pT> vs truth pT, and reco pT shapes per truth pT bin (R=0.4).
//   page 2: jet response <reco/truth pT> vs truth jet pT, per radius.
//   page 3: x_J shapes at R=0.4 per pT bin, truth (top) and reco region A in-situ trees (bottom).
// Photon-sample histograms are cross-section weighted (drawer.h); Jet12_long is one sample. The two
// differ in process mix (Jet12_long photons are QCD 2->2 fragmentation/radiation photons), so x_J
// is expected to differ; the detector responses should not.

const char * hists_dir = ana::path("hists");
const char * insitu_input_dir = ana::path("insitu/inputs");
const int nPt = ana::nPtBinsUsed;
const int G = ana::nUnfoldXjBins + 2;
const int irPho = 2; // photon response from the R=0.4 pairing
map<int,double> photonXsec = {{5,146359.3},{10,6944.675},{20,130.4461}};

TH1 * getPhoton(const char * hname) {
  TH1 * sum = nullptr;
  for (int s : {5, 10, 20}) {
    TFile * f = TFile::Open(Form("%s/Photon%d_pythia_nominal_unfolding.root", hists_dir, s), "READ");
    if (!f || f->IsZombie()) continue;
    TH1 * h = (TH1*)f->Get(hname);
    if (h && h->GetEntries() > 0) {
      h = (TH1*)h->Clone(Form("%s_P%d", hname, s)); h->SetDirectory(0); h->Scale(photonXsec[s]);
      if (!sum) { sum = (TH1*)h->Clone(Form("%s_pho", hname)); sum->SetDirectory(0); } else sum->Add(h);
    }
    f->Close();
  }
  return sum;
}
TH1 * getJet12(const char * hname) {
  TFile * f = TFile::Open(Form("%s/Jet12_long_pythia_nominal_unfolding.root", hists_dir), "READ");
  TH1 * h = (TH1*)f->Get(hname); h = (TH1*)h->Clone(Form("%s_j12", hname)); h->SetDirectory(0); f->Close();
  return h;
}
// <reco>/truth per truth bin of a (x = reco, y = truth) Hresponse; truth bins merged in groups of
// `group` to beat down Jet12_long statistics.
TH1D * responseProfile(TH2 * h, const char * name, double lo, double hi, int group) {
  vector<double> edges;
  int j1 = h->GetYaxis()->FindBin(lo + 1e-3), j2 = h->GetYaxis()->FindBin(hi - 1e-3);
  for (int j = j1; j <= j2 + 1; j += group) edges.push_back(h->GetYaxis()->GetBinLowEdge(std::min(j, j2 + 1)));
  if (edges.back() < h->GetYaxis()->GetBinUpEdge(j2)) edges.push_back(h->GetYaxis()->GetBinUpEdge(j2));
  TH1D * q = new TH1D(name, "", edges.size() - 1, edges.data());
  for (int k = 0; k + 1 < (int)edges.size(); k++) {
    int a = h->GetYaxis()->FindBin(edges[k] + 1e-3), b = h->GetYaxis()->FindBin(edges[k+1] - 1e-3);
    TH1D * px = h->ProjectionX("_px", a, b); TH1D * py = h->ProjectionY("_py", a, b);
    if (px->GetEffectiveEntries() >= 10) {
      // mean reco over mean truth in this truth slice, error from the reco mean
      q->SetBinContent(k+1, px->GetMean()/py->GetMean()); q->SetBinError(k+1, px->GetMeanError()/py->GetMean());
    }
    delete px; delete py;
  }
  return q;
}
// x_J shape of reported pT bin ipt from a flattened histogram, unit area / bin width.
TH1D * shapeFromGlobal(TH1 * g, int ipt, const char * name) {
  TH1D * h = new TH1D(name, ";x_{J#gamma};(1/N) dN/dx_{J#gamma}", ana::nUnfoldXjBins, ana::unfoldXjBins);
  for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
    int b = g->FindBin((ipt + ana::firstUsedPtBin)*G + ixj + 1);
    h->SetBinContent(ixj+1, g->GetBinContent(b)); h->SetBinError(ixj+1, g->GetBinError(b));
  }
  if (h->Integral() > 0) h->Scale(1.0/h->Integral(), "width");
  return h;
}
TH1D * shapeFromInsitu(bool jet12, int ir, int ipt, float floorXj, const char * name) {
  TH1D * h = new TH1D(name, ";x_{J#gamma};(1/N) dN/dx_{J#gamma}", ana::nUnfoldXjBins, ana::unfoldXjBins); h->Sumw2();
  vector<pair<string,double>> files;
  if (jet12) files = {{insitu_utility::insituFilename(insitu_input_dir, "Jet12_long", "pythia", "nominal"), 1.0}};
  else for (int s : {5, 10, 20}) files.push_back({insitu_utility::insituFilename(insitu_input_dir, Form("Photon%d", s), "pythia", "nominal"), photonXsec[s]});
  for (auto & fs : files) {
    TFile * f = TFile::Open(fs.first.c_str(), "READ"); if (!f || f->IsZombie()) continue;
    TTree * t = (TTree*)f->Get("insitutree");
    Float_t pho_pt, jet_pt, w; Int_t abcd, evIr;
    t->SetBranchAddress("pho_pt", &pho_pt); t->SetBranchAddress("jet_pt", &jet_pt);
    t->SetBranchAddress("abcd", &abcd); t->SetBranchAddress("weight", &w); t->SetBranchAddress("ir", &evIr);
    for (Long64_t e = 0; e < t->GetEntries(); e++) {
      t->GetEntry(e);
      if (abcd != 0 || evIr != ir || ana::findPtBin(pho_pt) != ipt + ana::firstUsedPtBin) continue;
      double x = jet_pt/pho_pt; if (x < floorXj) continue;
      h->Fill(x, fs.second*w);
    }
    f->Close();
  }
  if (h->Integral() > 0) h->Scale(1.0/h->Integral(), "width");
  return h;
}
void style(TH1 * h, int col, int mk) { h->SetLineColor(col); h->SetMarkerColor(col); h->SetMarkerStyle(mk); h->SetLineWidth(2); h->SetMarkerSize(0.9); }
void pad(double l = 0.16) { gPad->SetLeftMargin(l); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05); gPad->SetBottomMargin(0.12); gPad->SetTicks(1,1); }

void compare_jet12_photon() {
  gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
  TH1::AddDirectory(kFALSE);
  const int colP = kBlue+1, colJ = kOrange+7;
  string pdf = "compare_jet12_photon.pdf";
  TCanvas * c = new TCanvas("c", "", 1500, 800);
  c->SaveAs((pdf + "[").c_str());

  // ---- page 1: photon response ----
  TH2 * pP = (TH2*)getPhoton(Form("hphoresponse%d", irPho));
  TH2 * pJ = (TH2*)getJet12(Form("hphoresponse%d", irPho));
  printf("Matched region-A photons (R=0.4 pairing): Photon samples %.0f entries, Jet12_long %.0f entries\n", pP->GetEntries(), pJ->GetEntries());
  c->Divide(4, 1, 0.002, 0.002);
  c->cd(1); pad(0.18);
  TH1D * qP = responseProfile(pP, "qP", 13, 40, 1), * qJ = responseProfile(pJ, "qJ", 13, 40, 3);
  style(qP, colP, 20); style(qJ, colJ, 21);
  qP->SetTitle(";truth photon p_{T} [GeV];#LTreco p_{T}#GT / #LTtruth p_{T}#GT");
  qP->SetMinimum(0.85); qP->SetMaximum(1.1); qP->GetYaxis()->SetTitleOffset(1.7);
  qP->Draw("e"); qJ->Draw("e same");
  TLine * one = new TLine(13, 1, 40, 1); one->SetLineStyle(2); one->Draw();
  insitu_utility::drawSPhenixLabel({"Pythia8 MC, photon response"}, {"truth-matched region-A pairs"}, .22, .9, 16, gPad->GetWh());
  TLegend * l = new TLegend(0.22, 0.16, 0.95, 0.30); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.042);
  l->AddEntry(qP, "Photon5/10/20", "p"); l->AddEntry(qJ, "Jet12_long (3 GeV bins)", "p"); l->Draw();
  printf("\nPhoton response <reco>/<truth> per reported truth pT bin:\n");
  for (int ipt = 0; ipt < nPt; ipt++) {
    double lo = ana::ptBinsUsed[ipt], hi = ana::ptBinsUsed[ipt+1];
    TH1D * rP = responseProfile(pP, Form("rP%d", ipt), lo, hi, 100), * rJ = responseProfile(pJ, Form("rJ%d", ipt), lo, hi, 100);
    printf("  %2.0f-%2.0f GeV: Photon %.4f +- %.4f   Jet12_long %.4f +- %.4f   ratio J/P %.4f +- %.4f\n", lo, hi,
        rP->GetBinContent(1), rP->GetBinError(1), rJ->GetBinContent(1), rJ->GetBinError(1),
        rJ->GetBinContent(1)/rP->GetBinContent(1), rJ->GetBinError(1)/rP->GetBinContent(1));
    c->cd(ipt+2); pad(0.15);
    int j1 = pP->GetYaxis()->FindBin(lo + 1e-3), j2 = pP->GetYaxis()->FindBin(hi - 1e-3);
    TH1D * a = pP->ProjectionX(Form("prP%d", ipt), j1, j2), * b = pJ->ProjectionX(Form("prJ%d", ipt), j1, j2);
    a->Scale(1.0/a->Integral()); b->Scale(1.0/b->Integral());
    style(a, colP, 20); style(b, colJ, 21);
    a->SetTitle(";reco photon p_{T} [GeV];normalized"); a->GetXaxis()->SetRangeUser(lo - 8, hi + 8);
    a->SetMaximum(1.5*std::max(a->GetMaximum(), b->GetMaximum())); a->SetMinimum(0); a->GetYaxis()->SetTitleOffset(1.5);
    a->Draw("hist e"); b->Draw("e same");
    TLegend * l2 = new TLegend(0.45, 0.66, 0.97, 0.80); l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextSize(0.04);
    l2->AddEntry(a, Form("Photon #LTreco#GT=%.2f", a->GetMean()), "l"); l2->AddEntry(b, Form("Jet12_long #LTreco#GT=%.2f", b->GetMean()), "p"); l2->Draw();
    insitu_utility::drawSPhenixLabel({"truth photon in"}, {Form("%.0f GeV < p_{T}^{#gamma,truth} < %.0f GeV", lo, hi)}, .19, .9, 16, gPad->GetWh());
  }
  c->SaveAs(pdf.c_str());

  // ---- page 2: jet response per radius ----
  c->Clear(); c->Divide(4, 2, 0.002, 0.002);
  printf("\nJet response <reco>/<truth>, truth jet pT 15-35 GeV:\n");
  for (int ir = 0; ir < ana::nJetR; ir++) {
    c->cd(ir+1); pad();
    TH2 * hP = (TH2*)getPhoton(Form("hjetresponse%d", ir)), * hJ = (TH2*)getJet12(Form("hjetresponse%d", ir));
    TH1D * jP = responseProfile(hP, Form("jP%d", ir), 8, 45, 1), * jJ = responseProfile(hJ, Form("jJ%d", ir), 8, 45, 2);
    style(jP, colP, 20); style(jJ, colJ, 21);
    jP->SetTitle(";truth jet p_{T} [GeV];#LTreco p_{T}#GT / #LTtruth p_{T}#GT");
    jP->SetMinimum(0.8); jP->SetMaximum(1.3); jP->GetYaxis()->SetTitleOffset(1.4);
    jP->Draw("e"); jJ->Draw("e same");
    TLine * o = new TLine(8, 1, 45, 1); o->SetLineStyle(2); o->Draw();
    insitu_utility::drawSPhenixLabel({"jet response (matched)"}, {Form("Jet R=%.1f", ana::JetRs[ir])}, .2, .9, 15, gPad->GetWh()*gPad->GetHNDC());
    TH1D * aP = responseProfile(hP, "aP", 15, 35, 100), * aJ = responseProfile(hJ, "aJ", 15, 35, 100);
    printf("  R=%.1f: Photon %.4f +- %.4f   Jet12_long %.4f +- %.4f\n", ana::JetRs[ir], aP->GetBinContent(1), aP->GetBinError(1), aJ->GetBinContent(1), aJ->GetBinError(1));
    if (ir == 0) { TLegend * lj = new TLegend(0.2, 0.16, 0.9, 0.30); lj->SetBorderSize(0); lj->SetFillStyle(0); lj->SetTextSize(0.05);
      lj->AddEntry(jP, "Photon5/10/20", "p"); lj->AddEntry(jJ, "Jet12_long", "p"); lj->Draw(); }
  }
  c->SaveAs(pdf.c_str());

  // ---- page 3: x_J shapes at R=0.4 ----
  c->Clear(); c->Divide(3, 2, 0.002, 0.002);
  const int ir = 2;
  TH1 * gtP = getPhoton(Form("htruthxj%d", ir)), * gtJ = getJet12(Form("htruthxj%d", ir));
  printf("\n<x_J> at R=0.4 (binned): truth Photon, Jet12_long | reco region A Photon, Jet12_long\n");
  for (int row = 0; row < 2; row++) for (int ipt = 0; ipt < nPt; ipt++) {
    c->cd(row*3 + ipt + 1); pad(0.14);
    float fl = insitu_utility::lowXjFloor(ir, ana::ptBinsUsed[ipt]);
    TH1D * a = row == 0 ? shapeFromGlobal(gtP, ipt, Form("tP%d", ipt)) : shapeFromInsitu(false, ir, ipt, fl, Form("rP%d", ipt));
    TH1D * b = row == 0 ? shapeFromGlobal(gtJ, ipt, Form("tJ%d", ipt)) : shapeFromInsitu(true, ir, ipt, fl, Form("rJ%d", ipt));
    style(a, colP, 20); style(b, colJ, 21);
    a->SetMaximum(1.5*std::max(a->GetMaximum(), b->GetMaximum())); a->SetMinimum(0);
    a->Draw("hist e"); b->Draw("e same");
    TLegend * l3 = new TLegend(0.5, 0.66, 0.97, 0.80); l3->SetBorderSize(0); l3->SetFillStyle(0); l3->SetTextSize(0.042);
    l3->AddEntry(a, Form("Photon #LTx_{J}#GT=%.3f", a->GetMean()), "l"); l3->AddEntry(b, Form("Jet12_long #LTx_{J}#GT=%.3f", b->GetMean()), "p"); l3->Draw();
    insitu_utility::drawSPhenixLabel({row == 0 ? "truth level" : "reco level, region A"},
        {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1]), "Jet R=0.4"}, .17, .9, 15, gPad->GetWh()*gPad->GetHNDC());
    printf("  %s %2.0f-%2.0f: Photon %.4f  Jet12_long %.4f\n", row == 0 ? "truth" : "reco ", ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1], a->GetMean(), b->GetMean());
  }
  c->SaveAs(pdf.c_str());
  c->SaveAs((pdf + "]").c_str());
  cout << "Wrote " << pdf << endl;
}
