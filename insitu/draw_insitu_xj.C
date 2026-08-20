#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
using namespace std;

// ana::findPtBin/getPurity/etc. are implemented in ana.cc, compiled into
// libgammajet_unfold.so - load it explicitly so cling resolves those symbols against
// the real compiled definitions rather than misbinding (see draw_purity_corrected.C).
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Nominal jet radius index (R=0.4), matching every other macro in this directory
// (e.g. grid_insitu.C/grid_insitu_jet12.C's default ir=2) - the insitutree files read
// below were themselves filled at this fixed R by unfolder.cc.
const int ir = 2;

// sPHENIX label block: bold-italic "sPHENIX Internal" title, then one line per sample,
// then one line per feature - same text/font convention and stacking formula as
// drawer::drawAll() (see src/drawer.cc) and grid_insitu.C's own drawSPhenixLabel(),
// reimplemented locally so this self-contained macro doesn't have to construct a full
// drawer (which opens a batch of unrelated unfolding-output files it has no other use for).
void drawSPhenixLabel(vector<string> samples, vector<string> features, float drawx, float drawy, int fontsize, float csize) {
  float titlescale = 1.25;
  float subtitlescale = 1.25;
  float ydiff = fontsize * 0.0017 * 700.0/csize;
  auto drawOne = [&](const char * text, float xp, float yp, int size) {
    TLatex * tex = new TLatex(xp, yp, text);
    tex->SetTextFont(43);
    tex->SetTextSize(size);
    tex->SetTextColor(kBlack);
    tex->SetLineWidth(1);
    tex->SetNDC();
    tex->Draw();
  };
  drawOne("#bf{#it{sPHENIX}} #kern[0.5]{Internal}", drawx, drawy, (int)(fontsize*titlescale));
  for (unsigned i = 0; i < samples.size(); i++) {
    drawOne(samples[i].c_str(), drawx, drawy-ydiff*subtitlescale*(i+1), (int)(fontsize*subtitlescale));
  }
  for (unsigned i = 0; i < features.size(); i++) {
    drawOne(features[i].c_str(), drawx, drawy-ydiff*subtitlescale*samples.size()-ydiff*(i+1)*subtitlescale, fontsize);
  }
}

// Reads the insitutree (pho_pt, jet_pt, abcd) written by unfolder.cc into
// gammajet_unfold/insitu/, and builds x_{J#gamma} = jet_pt/pho_pt histograms for:
//   1. Region A (signal region) in Data
//   2. Region A in Pythia8 gamma+jet MC (Photon5+10+20 combined, cross-section weighted)
//   3. Purity-corrected Data (Region A minus the scaled Region C background)
//
// Physics-level comparison uses only ana::ptBinsUsed (15-35 GeV, ana::firstUsedPtBin
// through +ana::nPtBinsUsed) - same restriction as drawing/draw_purity_corrected.C,
// dropping both the 13-15 GeV migration buffer bin and the 35-100 GeV overflow bin
// (ana::ptBins[0] and [nPtBins-1] respectively - see ana.h's comment).
const int nPtBinsUsed = ana::nPtBinsUsed;
const char * systag = "nominal";
const char * insitu_dir = "/home/samson72/sphnx/gammajet_unfold/insitu";

// Cross-section weights for combining the Photon5/10/20 samples - same numbers as
// drawer.h's scalemap[isphoton=1][sample] for sim="pythia".
map<int,double> photon_scale = {{5,146359.3},{10,6944.675},{20,130.4461}};

// Fills one x_{J} histogram per photon-pT bin (ana::ptBins binning) from a single
// insitutree file, selecting only the requested ABCD region, scaled by `weight`.
// Returns nullptr entries (left as empty histograms) if the file/tree is missing.
vector<TH1D*> fillXjByPtBin(const char * filename, int abcdSelect, double weight, const char * tag) {
  vector<TH1D*> h(ana::nPtBins);
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    h[ipt] = new TH1D(Form("hxj_%s_pt%d", tag, ipt), ";x_{J#gamma};Counts",
        ana::nUnfoldXjBins, ana::unfoldXjBins);
    h[ipt]->Sumw2();
  }

  TFile * f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie()) {
    cout << "WARNING: could not open " << filename << " - leaving histograms empty." << endl;
    return h;
  }
  TTree * t = (TTree*)f->Get("insitutree");

  Float_t pho_pt, jet_pt;
  Int_t abcd;
  t->SetBranchAddress("pho_pt", &pho_pt);
  t->SetBranchAddress("jet_pt", &jet_pt);
  t->SetBranchAddress("abcd", &abcd);

  Long64_t nentries = t->GetEntries();
  for (Long64_t e = 0; e < nentries; e++) {
    t->GetEntry(e);
    if (abcd != abcdSelect) continue;
    int ipt = ana::findPtBin(pho_pt);
    if (ipt < 0) continue;
    h[ipt]->Fill(jet_pt/pho_pt, weight);
  }
  f->Close();
  return h;
}

// Adds a second sample's per-pT-bin histograms into the first, in place.
void addInto(vector<TH1D*> & total, const vector<TH1D*> & add) {
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    total[ipt]->Add(add[ipt]);
  }
}

// signal(xJ) = A(xJ) - (1-P)*(N_A/N_C)*C(xJ), P = ana::getPurity(ptlow,pthigh,systag) -
// same formula as drawing/draw_purity_corrected.C's purityCorrectP(), one purity value
// (and its bootstrap errors) per pT bin, fully correlated across xJ bins in that bin.
TH1D * purityCorrectP(TH1D * A, TH1D * C, float p, float pErrLow, float pErrHigh, const char * name) {
  float NA = A->Integral();
  float NC = C->Integral();
  TH1D * h = (TH1D*)A->Clone(name);
  if (NC <= 0) {
    cout << "WARNING: " << name << " has zero region-C statistics - cannot cross-normalize, using raw region A." << endl;
    return h;
  }
  float scale = (1-p)*(NA/NC);
  for (int i = 1; i <= A->GetNbinsX(); i++) {
    float a  = A->GetBinContent(i);
    float ae = A->GetBinError(i);
    float c  = C->GetBinContent(i);
    float ce = C->GetBinError(i);
    float dPurityLow  = (NA/NC)*c*pErrLow;
    float dPurityHigh = (NA/NC)*c*pErrHigh;
    float content = a - scale*c;
    float errLow  = sqrt(ae*ae + pow(scale*ce,2) + pow(dPurityLow,2));
    float errHigh = sqrt(ae*ae + pow(scale*ce,2) + pow(dPurityHigh,2));
    h->SetBinContent(i, content);
    h->SetBinError(i, std::max(errLow, errHigh));
  }
  return h;
}

// x_{J} bins are non-uniform (ana::unfoldXjBins) - divide by bin width so the
// comparison plot shows a density, not raw counts with an artificial shelf where the
// bin width changes.
TH1D * densityForDisplay(TH1D * h, const char * name) {
  TH1D * hd = (TH1D*)h->Clone(name);
  hd->Scale(1., "width");
  hd->GetYaxis()->SetTitle("Counts / bin width");
  return hd;
}

void draw_insitu_xj() {
  gStyle->SetOptStat(0);

  // -----------------------------
  // 1. Region A in Data
  // -----------------------------
  vector<TH1D*> hA_data = fillXjByPtBin(Form("%s/Data_%s_insitu.root", insitu_dir, systag), 0, 1.0, "A_data");

  // -----------------------------
  // 2. Region A in Pythia (Photon5+10+20, cross-section weighted)
  // -----------------------------
  vector<TH1D*> hA_pythia = fillXjByPtBin(Form("%s/Photon5_pythia_%s_insitu.root", insitu_dir, systag), 0, photon_scale[5], "A_pythia_p5");
  addInto(hA_pythia, fillXjByPtBin(Form("%s/Photon10_pythia_%s_insitu.root", insitu_dir, systag), 0, photon_scale[10], "A_pythia_p10"));
  addInto(hA_pythia, fillXjByPtBin(Form("%s/Photon20_pythia_%s_insitu.root", insitu_dir, systag), 0, photon_scale[20], "A_pythia_p20"));

  // -----------------------------
  // 3. Purity-corrected Data (Region A - scaled Region C)
  // -----------------------------
  vector<TH1D*> hC_data = fillXjByPtBin(Form("%s/Data_%s_insitu.root", insitu_dir, systag), 2, 1.0, "C_data");

  vector<TH1D*> hCorrected(ana::nPtBins);
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin + nPtBinsUsed; ipt++) {
    float ptlow  = ana::ptBins[ipt];
    float pthigh = ana::ptBins[ipt+1];
    float p        = ana::getPurity(ptlow, pthigh, systag);
    float pErrLow  = ana::getPurityErrorLow(ptlow, pthigh, systag);
    float pErrHigh = ana::getPurityErrorHigh(ptlow, pthigh, systag);
    hCorrected[ipt] = purityCorrectP(hA_data[ipt], hC_data[ipt], p, pErrLow, pErrHigh,
        Form("hxjcorrected_data_pt%d", ipt));
  }

  // -----------------------------
  // Combine the used pT bins (15-35 GeV, ana::ptBinsUsed) into one histogram per
  // version - hA_data/hA_pythia/hCorrected are indexed by ana::findPtBin's raw
  // (unrestricted) bin index, so start/iterate from ana::firstUsedPtBin.
  // -----------------------------
  TH1D * hxjA_data       = (TH1D*)hA_data[ana::firstUsedPtBin]->Clone("hxjA_data");
  TH1D * hxjA_pythia     = (TH1D*)hA_pythia[ana::firstUsedPtBin]->Clone("hxjA_pythia");
  TH1D * hxjcorrected_data = (TH1D*)hCorrected[ana::firstUsedPtBin]->Clone("hxjcorrected_data");
  for (int ipt = ana::firstUsedPtBin+1; ipt < ana::firstUsedPtBin + nPtBinsUsed; ipt++) {
    hxjA_data->Add(hA_data[ipt]);
    hxjA_pythia->Add(hA_pythia[ipt]);
    hxjcorrected_data->Add(hCorrected[ipt]);
  }
  hxjA_data->SetTitle("Region A (Data)");
  hxjA_pythia->SetTitle("Region A (Pythia8 #gamma+jet)");
  hxjcorrected_data->SetTitle("Purity-corrected (Data, A-C)");

  // -----------------------------
  // Save
  // -----------------------------
  const char * outRootPath = Form("%s/insitu_xj_comparison.root", insitu_dir);
  TFile * fout = TFile::Open(outRootPath, "RECREATE");
  hxjA_data->Write();
  hxjA_pythia->Write();
  hxjcorrected_data->Write();
  fout->Close();
  cout << "Wrote " << outRootPath << endl;

  // -----------------------------
  // Comparison plot (shape-normalized density, since the three are at different
  // absolute scales - Data counts vs. MC cross-section-weighted counts).
  // -----------------------------
  TH1D * dispA_data       = densityForDisplay(hxjA_data,       "hxjA_data_disp");
  TH1D * dispA_pythia     = densityForDisplay(hxjA_pythia,     "hxjA_pythia_disp");
  TH1D * dispCorrected    = densityForDisplay(hxjcorrected_data,"hxjcorrected_data_disp");
  if (dispA_data->Integral() > 0)    dispA_data->Scale(1./dispA_data->Integral());
  if (dispA_pythia->Integral() > 0)  dispA_pythia->Scale(1./dispA_pythia->Integral());
  if (dispCorrected->Integral() > 0) dispCorrected->Scale(1./dispCorrected->Integral());

  TCanvas * c = new TCanvas("c","",700,700);
  gPad->SetLeftMargin(.15);
  gPad->SetTicks();
  dispA_data->SetLineColor(kBlack);
  dispA_data->SetMarkerColor(kBlack);
  dispA_data->SetMarkerStyle(20);
  dispA_data->SetLineWidth(2);
  dispA_data->GetXaxis()->SetTitle("x_{J#gamma}");
  dispA_data->GetYaxis()->SetTitle("Shape-normalized counts / bin width");
  // Extra headroom (vs. the plain-comparison pages elsewhere in this directory) to fit
  // the sPHENIX label block and the mean-value text below the legend without either
  // colliding with the curves themselves - both are drawn inside the frame per this
  // project's plotting convention, never in an outer pad margin.
  dispA_data->GetYaxis()->SetRangeUser(0, std::max({dispA_data->GetMaximum(),
        dispA_pythia->GetMaximum(), dispCorrected->GetMaximum()})*1.7);
  dispA_data->Draw("p e");

  dispA_pythia->SetLineColor(kAzure+2);
  dispA_pythia->SetMarkerColor(kAzure+2);
  dispA_pythia->SetMarkerStyle(21);
  dispA_pythia->SetLineWidth(2);
  dispA_pythia->Draw("p e same");

  dispCorrected->SetLineColor(kRed);
  dispCorrected->SetMarkerColor(kRed);
  dispCorrected->SetMarkerStyle(22);
  dispCorrected->SetLineWidth(2);
  dispCorrected->Draw("p e same");

  TLegend * l = new TLegend(.5,.68,.85,.85);
  l->SetLineWidth(0);
  l->AddEntry(dispA_data,    "Region A (Data)");
  l->AddEntry(dispA_pythia,  "Region A (Pythia8)");
  l->AddEntry(dispCorrected, "Purity-corrected Data (A-C)");
  l->Draw();

  // Means are taken from the raw (non-density-scaled) histograms, not the display
  // (density) versions above - dividing by the non-uniform bin width reweights the mean
  // calculation by 1/width per bin and biases it, since ana::unfoldXjBins isn't uniform.
  TLatex * texMean = new TLatex();
  texMean->SetNDC();
  texMean->SetTextFont(43);
  texMean->SetTextSize(16);
  texMean->SetTextColor(kBlack);
  texMean->DrawLatex(.5, .63, Form("#LTx_{J#gamma}#GT_{Data} = %.3f #pm %.3f",
        hxjA_data->GetMean(), hxjA_data->GetMeanError()));
  texMean->SetTextColor(kAzure+2);
  texMean->DrawLatex(.5, .585, Form("#LTx_{J#gamma}#GT_{Pythia8} = %.3f #pm %.3f",
        hxjA_pythia->GetMean(), hxjA_pythia->GetMeanError()));
  texMean->SetTextColor(kRed);
  texMean->DrawLatex(.5, .54, Form("#LTx_{J#gamma}#GT_{corr.} = %.3f #pm %.3f",
        hxjcorrected_data->GetMean(), hxjcorrected_data->GetMeanError()));

  drawSPhenixLabel({"p+p Run24 Data"}, {
      Form("Jet R=%.1f", ana::JetRs[ir]),
      Form("|#eta^{#gamma}|<%.1f, |#eta^{jet}|<%.1f", ana::etacut, ana::etacut-ana::JetRs[ir]),
      Form("#Delta#phi>%.0f#pi/%.0f", ana::oppnum, ana::oppden)
    }, .18, .85, 16, gPad->GetWh());

  const char * outPdfPath = Form("%s/insitu_xj_comparison.pdf", insitu_dir);
  c->SaveAs(outPdfPath);
  cout << "Wrote " << outPdfPath << endl;
}
