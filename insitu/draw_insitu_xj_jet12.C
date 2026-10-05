#include "../src/ana.h"
#include "../src/unfold_utility.h"
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
using namespace std;

// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Non-purity-corrected companion to draw_insitu_xj.C (counterpart of grid_insitu_jet12.C):
// raw region A xJ in Data and in Jet12_long, used pT bins only.
const int nPtBinsUsed = ana::nPtBinsUsed;
const char * systag = "nominal";
const char * insitu_input_dir  = ana::path("insitu/inputs");
const char * insitu_output_dir = ana::path("insitu/output");
const char * insitu_pdf_dir    = ana::path("insitu/pdfs");
// Nominal R=0.4 (the tree holds every radius).
const int ir = 2;

// Jet12 cross-section weight (drawer.h's scalemap).
map<int,double> jet_scale = {{12,3.997e+06}};

// One xJ histogram per pT bin from one insitu tree for one ABCD region, scaled by weight times
// the tree's weight branch (1 for Data). Empty histograms if the file is missing.
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

void draw_insitu_xj_jet12() {
  gStyle->SetOptStat(0);

  // -----------------------------
  // 1. Region A in Data
  // -----------------------------
  vector<TH1D*> hA_data = fillXjByPtBin(Form("%s/Data_%s_insitu.root", insitu_input_dir, systag), 0, 1.0, "A_data");

  // -----------------------------
  // 2. Region A in Jet12_long
  // -----------------------------
  vector<TH1D*> hA_jet12 = fillXjByPtBin(Form("%s/Jet12_long_pythia_%s_insitu.root", insitu_input_dir, systag), 0, jet_scale[12], "A_jet12");

  // -----------------------------
  // Combine the used pT bins (indexed by findPtBin, so start at firstUsedPtBin)
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
  // Comparison plot (shape-normalized densities)
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
