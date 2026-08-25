#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/unfold_utility.h"
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
using namespace std;

// ana::findPtBin/etc. are implemented in ana.cc, compiled into libgammajet_unfold.so -
// load it explicitly so cling resolves those symbols against the real compiled
// definitions rather than misbinding (see draw_insitu_xj.C).
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Non-purity-corrected companion to draw_insitu_xj.C: reads the same insitutree
// (pho_pt, jet_pt, abcd) written by unfolder.cc, and builds x_{J#gamma} = jet_pt/pho_pt
// histograms for:
//   1. Region A (signal region) in Data - raw, no A-C background subtraction.
//   2. Region A in "Jet12_long" (a single Pythia8 QCD-dijet-triggered MC sample, no
//      truth-level jet-pT cut, at the standard gammajet/trees path - see
//      src/treeuser.h) - also raw, no ABCD subtraction on the MC side either.
//
// Unlike draw_insitu_xj.C there is no purity-corrected histogram here at all - this
// macro is the shape-comparison counterpart of grid_insitu_jet12.C's cross-check, and
// intentionally skips ana::getPurity()/Region C entirely on both sides.
//
// Physics-level comparison uses only ana::ptBinsUsed (15-35 GeV, ana::firstUsedPtBin
// through +ana::nPtBinsUsed) - same restriction as draw_insitu_xj.C and
// grid_insitu_jet12.C.
const int nPtBinsUsed = ana::nPtBinsUsed;
const char * systag = "nominal";
// insitu/ is split into inputs/ (the raw Data/Jet12_long insitu ntuples, written by
// unfolder.h's production pipeline), output/ (this macro's own .root output), and
// pdfs/ (its .pdf output).
const char * insitu_input_dir  = "/home/samson72/sphnx/gammajet_unfold/insitu/inputs";
const char * insitu_output_dir = "/home/samson72/sphnx/gammajet_unfold/insitu/output";
const char * insitu_pdf_dir    = "/home/samson72/sphnx/gammajet_unfold/insitu/pdfs";
// Nominal jet radius index (R=0.4), matching draw_insitu_xj.C - the insitutree files
// read below now hold every jet radius together (one row per radius an event paired
// at), so this needs to be filtered on rather than assumed.
const int ir = 2;

// Cross-section weight for the Jet12 pythia sample - same number as drawer.h's
// scalemap[isphoton=0][12] for sim="pythia" (see grid_insitu_jet12.C).
map<int,double> jet_scale = {{12,3.997e+06}};

// Fills one x_{J} histogram per photon-pT bin (ana::ptBins binning) from a single
// insitutree file, selecting only the requested ABCD region, scaled by `weight` (the
// sample's cross-section weight, 1.0 for Data) times the tree's own "weight" branch
// (the vz/cluster_pt mcWeight from unfolder.cc - always 1.0 for Data, so this is a
// no-op there and only reweights MC).
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

  Float_t pho_pt, jet_pt, mcWeight;
  Int_t abcd, evIr;
  t->SetBranchAddress("pho_pt", &pho_pt);
  t->SetBranchAddress("jet_pt", &jet_pt);
  t->SetBranchAddress("abcd", &abcd);
  t->SetBranchAddress("weight", &mcWeight);
  t->SetBranchAddress("ir", &evIr);

  Long64_t nentries = t->GetEntries();
  for (Long64_t e = 0; e < nentries; e++) {
    t->GetEntry(e);
    if (abcd != abcdSelect) continue;
    if (evIr != ir) continue;
    int ipt = ana::findPtBin(pho_pt);
    if (ipt < 0) continue;
    h[ipt]->Fill(jet_pt/pho_pt, weight*mcWeight);
  }
  f->Close();
  return h;
}

// x_{J} bins are non-uniform (ana::unfoldXjBins) - divide by bin width so the
// comparison plot shows a density, not raw counts with an artificial shelf where the
// bin width changes.
// densityForDisplay now lives in unfold_utility - see src/unfold_utility.h.

void draw_insitu_xj_jet12() {
  gStyle->SetOptStat(0);

  // -----------------------------
  // 1. Region A in Data (raw, no purity correction)
  // -----------------------------
  vector<TH1D*> hA_data = fillXjByPtBin(Form("%s/Data_%s_insitu.root", insitu_input_dir, systag), 0, 1.0, "A_data");

  // -----------------------------
  // 2. Region A in Pythia8 Jet12_long (raw, no purity correction)
  // -----------------------------
  vector<TH1D*> hA_jet12 = fillXjByPtBin(Form("%s/Jet12_long_pythia_%s_insitu.root", insitu_input_dir, systag), 0, jet_scale[12], "A_jet12");

  // -----------------------------
  // Combine the used pT bins (15-35 GeV, ana::ptBinsUsed) into one histogram per
  // version - hA_data/hA_jet12 are indexed by ana::findPtBin's raw (unrestricted)
  // bin index, so start/iterate from ana::firstUsedPtBin.
  // -----------------------------
  TH1D * hxjA_data  = (TH1D*)hA_data[ana::firstUsedPtBin]->Clone("hxjA_data");
  TH1D * hxjA_jet12 = (TH1D*)hA_jet12[ana::firstUsedPtBin]->Clone("hxjA_jet12");
  for (int ipt = ana::firstUsedPtBin+1; ipt < ana::firstUsedPtBin + nPtBinsUsed; ipt++) {
    hxjA_data->Add(hA_data[ipt]);
    hxjA_jet12->Add(hA_jet12[ipt]);
  }
  hxjA_data->SetTitle("Region A (Data, raw)");
  hxjA_jet12->SetTitle("Region A (Pythia8 Jet12, raw)");

  // -----------------------------
  // Save
  // -----------------------------
  const char * outRootPath = Form("%s/insitu_xj_comparison_jet12.root", insitu_output_dir);
  TFile * fout = TFile::Open(outRootPath, "RECREATE");
  hxjA_data->Write();
  hxjA_jet12->Write();
  fout->Close();
  cout << "Wrote " << outRootPath << endl;

  // -----------------------------
  // Comparison plot (shape-normalized density, since the two are at different
  // absolute scales - Data counts vs. MC cross-section-weighted counts).
  // -----------------------------
  TH1D * dispA_data  = unfold_utility::densityForDisplay(hxjA_data,  "hxjA_data_disp");
  TH1D * dispA_jet12 = unfold_utility::densityForDisplay(hxjA_jet12, "hxjA_jet12_disp");
  if (dispA_data->Integral() > 0)  dispA_data->Scale(1./dispA_data->Integral());
  if (dispA_jet12->Integral() > 0) dispA_jet12->Scale(1./dispA_jet12->Integral());

  TCanvas * c = new TCanvas("c","",700,700);
  gPad->SetLeftMargin(.15);
  gPad->SetTicks();
  dispA_data->SetLineColor(kBlack);
  dispA_data->SetMarkerColor(kBlack);
  dispA_data->SetMarkerStyle(20);
  dispA_data->SetLineWidth(2);
  dispA_data->GetXaxis()->SetTitle("x_{J#gamma}");
  dispA_data->GetYaxis()->SetTitle("Shape-normalized counts / bin width");
  dispA_data->GetYaxis()->SetRangeUser(0, std::max(dispA_data->GetMaximum(), dispA_jet12->GetMaximum())*1.4);
  dispA_data->Draw("p e");

  dispA_jet12->SetLineColor(kAzure+2);
  dispA_jet12->SetMarkerColor(kAzure+2);
  dispA_jet12->SetMarkerStyle(21);
  dispA_jet12->SetLineWidth(2);
  dispA_jet12->Draw("p e same");

  TLegend * l = new TLegend(.5,.65,.85,.85);
  l->SetLineWidth(0);
  l->AddEntry(dispA_data,  "Region A (Data, raw)");
  l->AddEntry(dispA_jet12, "Region A (Pythia8 Jet12, raw)");
  l->Draw();

  const char * outPdfPath = Form("%s/insitu_xj_comparison_jet12.pdf", insitu_pdf_dir);
  c->SaveAs(outPdfPath);
  cout << "Wrote " << outPdfPath << endl;
}
