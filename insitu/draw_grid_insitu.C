#include "../src/ana.h"
#include "../src/insitu_utility.h"
#include "grid_insitu_pages.h"
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
using namespace std;

// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// In-situ scan pages from grid_insitu.C's gammajet-mode output, per jet radius: <x_J> vs pT and
// the per-pT x_J shapes, for region A and for the purity-corrected combination. Drawing only.
// tag: the grid_insitu_<tag>.root suffix (systag, or "shapechi2_<systag>" for the shape method).
// Output: insitu/pdfs/grid_insitu_<tag>.pdf
void draw_grid_insitu(string tag = "nominal") {
  gStyle->SetOptStat(0);
  const char * insitu_output_dir = ana::path("insitu/output");
  const char * insitu_pdf_dir    = ana::path("insitu/pdfs");
  const bool shape = tag.rfind("shapechi2_", 0) == 0;

  string infile = Form("%s/grid_insitu_%s.root", insitu_output_dir, tag.c_str());
  TFile * fin = TFile::Open(infile.c_str());
  if (!fin || fin->IsZombie()) {
    cout << "ERROR: cannot open " << infile << " - run grid_insitu.C first." << endl;
    return;
  }
  string pdfPathStr = Form("%s/grid_insitu_%s.pdf", insitu_pdf_dir, tag.c_str());
  const char * pdfPath = pdfPathStr.c_str();
  TCanvas * c = new TCanvas("c","",700,700);
  c->SaveAs(Form("%s[", pdfPath));

  for (int ir = 0; ir < ana::nJetR; ir++) {
    TDirectory * d = fin->GetDirectory(ana::rnames[ir]);
    if (!d) continue; // radius not scanned (grid_insitu onlyIr)
    TTree * t = (TTree*)d->Get("results");
    float paA, errLowA, errHighA, paCorr, errLowCorr, errHighCorr;
    t->SetBranchAddress("pa_regionA", &paA);
    t->SetBranchAddress("errLow_regionA", &errLowA);
    t->SetBranchAddress("errHigh_regionA", &errHighA);
    t->SetBranchAddress("pa_puritycorrected", &paCorr);
    t->SetBranchAddress("errLow_puritycorrected", &errLowCorr);
    t->SetBranchAddress("errHigh_puritycorrected", &errHighCorr);
    t->GetEntry(0);
    auto graph = [&](const char * name) { return (TGraphErrors*)d->Get(name); };
    auto hist  = [&](const char * prefix, int ipt) { return (TH1D*)d->Get(Form("%s_pt%d", prefix, ipt)); };

    drawJESPage(c, pdfPath, "Region A", ir, graph("gMeanMC"), graph("gMeanData_regionA_raw"),
        graph("gRatio_regionA_raw"), graph("gRatio_regionA_corrected"), paA, errLowA, errHighA, "", shape);
    for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) {
      drawXjPage(c, pdfPath, "Region A", ir, ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1],
          hist("hxjA_pythia", ipt), hist("hxjA_data_raw", ipt), hist("hxjA_data_corr", ipt),
          "Data (raw)", "Data (JES-corrected)");
    }
    drawJESPage(c, pdfPath, "Purity-corrected", ir, graph("gMeanMC"), graph("gMeanData_puritycorrected_raw"),
        graph("gRatio_puritycorrected_raw"), graph("gRatio_puritycorrected_corrected"), paCorr, errLowCorr, errHighCorr, "", shape);
    for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) {
      drawXjPage(c, pdfPath, "Purity-corrected", ir, ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1],
          hist("hxjA_pythia", ipt), hist("hxjcorrected_data_raw", ipt), hist("hxjcorrected_data_bestscale", ipt),
          "Data (raw, purity-corr.)", "Data (JES+purity-corr.)");
    }
  }

  c->SaveAs(Form("%s]", pdfPath));
  cout << "Wrote " << pdfPathStr << endl;
}
