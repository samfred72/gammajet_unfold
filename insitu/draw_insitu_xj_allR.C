#include "../src/ana.h"
#include "../src/insitu_utility.h"
#include "../src/unfold_utility.h"
#include <string>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
using namespace std;

// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// All-radii version of draw_insitu_xj.C: reco xJ for Data region A, Photon5+10+20 region A and
// purity-corrected Data, shape-normalized with a ratio to Pythia8; one page per radius, one
// column per used pT bin. Data jets are uncorrected (the scan at p_a = 1). Means are taken above
// insitu_utility::lowXjFloor (dashed line), as grid_insitu.C does, from the raw histograms.
//
// Usage: root -b -l -q 'draw_insitu_xj_allR.C("nominal")'

const char * insitu_input_dir  = ana::path("insitu/inputs");
const char * insitu_output_dir = ana::path("insitu/output");
const char * insitu_pdf_dir    = ana::path("insitu/pdfs");

// Cross-section weights (drawer.h's scalemap).
map<int,double> photon_scale = {{5,146359.3},{10,6944.675},{20,130.4461}};

typedef vector<vector<TH1D*>> HistGrid; // [ir][ipt]

HistGrid makeGrid(const char * tag) {
  HistGrid h(ana::nJetR, vector<TH1D*>(ana::nPtBins));
  for (int ir = 0; ir < ana::nJetR; ir++) {
    for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
      h[ir][ipt] = new TH1D(Form("hxj_%s_%s_pt%d", tag, ana::rnames[ir], ipt), ";x_{J#gamma};Counts",
          ana::nUnfoldXjBins, ana::unfoldXjBins);
      h[ir][ipt]->Sumw2();
    }
  }
  return h;
}

// One pass fills every radius (the tree's ir branch).
void fillGrid(HistGrid & h, const string & filename, int abcdSelect, double weight) {
  TFile * f = TFile::Open(filename.c_str(), "READ");
  if (!f || f->IsZombie()) {
    cout << "WARNING: could not open " << filename << " - leaving histograms empty." << endl;
    return;
  }
  TTree * t = (TTree*)f->Get("insitutree");
  Float_t pho_pt, jet_pt, mcWeight;
  Int_t abcd, evIr;
  t->SetBranchAddress("pho_pt", &pho_pt);
  t->SetBranchAddress("jet_pt", &jet_pt);
  t->SetBranchAddress("abcd", &abcd);
  t->SetBranchAddress("weight", &mcWeight);
  t->SetBranchAddress("ir", &evIr);
  Long64_t n = t->GetEntries();
  for (Long64_t e = 0; e < n; e++) {
    t->GetEntry(e);
    if (abcd != abcdSelect) continue;
    if (evIr < 0 || evIr >= ana::nJetR) continue;
    int ipt = ana::findPtBin(pho_pt);
    if (ipt < 0) continue;
    h[evIr][ipt]->Fill(jet_pt/pho_pt, weight*mcWeight);
  }
  f->Close();
}

// Zero NaN/Inf bins before drawing.
void sanitize(TH1D * h) {
  for (int b = 0; b <= h->GetNbinsX()+1; b++) {
    if (!std::isfinite(h->GetBinContent(b)) || !std::isfinite(h->GetBinError(b))) {
      h->SetBinContent(b, 0);
      h->SetBinError(b, 0);
    }
  }
}

// Bin-width density, then unit area.
TH1D * shapeDisplay(TH1D * h, const char * name) {
  TH1D * d = unfold_utility::densityForDisplay(h, name);
  if (d->Integral() > 0) d->Scale(1./d->Integral());
  sanitize(d);
  return d;
}

// Mean and error over bins with lower edge >= floor.
void meanAboveFloor(TH1D * h, double floor, double & mean, double & err) {
  double sw = 0, swx = 0;
  for (int b = 1; b <= h->GetNbinsX(); b++) {
    if (h->GetXaxis()->GetBinLowEdge(b) < floor - 1e-9) continue;
    sw  += h->GetBinContent(b);
    swx += h->GetBinContent(b) * h->GetXaxis()->GetBinCenter(b);
  }
  mean = (sw != 0) ? swx/sw : 0;
  double var = 0;
  for (int b = 1; b <= h->GetNbinsX(); b++) {
    if (h->GetXaxis()->GetBinLowEdge(b) < floor - 1e-9) continue;
    double dx = h->GetXaxis()->GetBinCenter(b) - mean;
    var += pow(h->GetBinError(b) * dx, 2);
  }
  err = (sw != 0) ? sqrt(var)/fabs(sw) : 0;
}

void styleHist(TH1D * h, int color, int marker) {
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetLineWidth(2);
}

void draw_insitu_xj_allR(string systag = "nominal") {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  const int nPtUsed = ana::nPtBinsUsed;

  HistGrid hA   = makeGrid("A_data");
  HistGrid hC   = makeGrid("C_data");
  HistGrid hAmc = makeGrid("A_pythia");
  fillGrid(hA, insitu_utility::insituFilename(insitu_input_dir, "Data", "", systag), 0, 1.0);
  fillGrid(hC, insitu_utility::insituFilename(insitu_input_dir, "Data", "", systag), 2, 1.0);
  for (int s : {5, 10, 20}) {
    fillGrid(hAmc, insitu_utility::insituFilename(insitu_input_dir, Form("Photon%d", s), "pythia", systag), 0, photon_scale[s]);
  }

  string pdfPath  = Form("%s/insitu_xj_allR_%s.pdf", insitu_pdf_dir, systag.c_str());
  string rootPath = Form("%s/insitu_xj_allR_%s.root", insitu_output_dir, systag.c_str());
  TFile * fout = TFile::Open(rootPath.c_str(), "RECREATE");

  const int colW = 650, canH = 900;
  TCanvas * c = new TCanvas("c", "", colW*nPtUsed, canH);
  c->SaveAs(Form("%s[", pdfPath.c_str()));

  for (int ir = 0; ir < ana::nJetR; ir++) {
    c->Clear();
    for (int k = 0; k < nPtUsed; k++) {
      int ipt = ana::firstUsedPtBin + k;
      float ptlow = ana::ptBins[ipt], pthigh = ana::ptBins[ipt+1];
      const char * tag = Form("%s_pt%d", ana::rnames[ir], ipt);
      TH1D * A = hA[ir][ipt], * C = hC[ir][ipt], * Amc = hAmc[ir][ipt];

      float pA        = ana::getPurity(ptlow, pthigh, systag, ir);
      float pAErrLow  = ana::getPurityErrorLow(ptlow, pthigh, systag, ir);
      float pAErrHigh = ana::getPurityErrorHigh(ptlow, pthigh, systag, ir);
      float pC        = ana::getPurityC(ptlow, pthigh, systag, ir);
      float pCErrLow  = ana::getPurityCErrorLow(ptlow, pthigh, systag, ir);
      float pCErrHigh = ana::getPurityCErrorHigh(ptlow, pthigh, systag, ir);
      TH1D * corr = unfold_utility::purityCorrect(A, C, pA, pAErrLow, pAErrHigh,
          pC, pCErrLow, pCErrHigh, Form("hxjcorrected_data_%s", tag));

      fout->cd();
      A->Write();
      C->Write();
      Amc->Write();
      if (corr) corr->Write();

      TH1D * dA   = shapeDisplay(A,   Form("%s_disp", A->GetName()));
      TH1D * dAmc = shapeDisplay(Amc, Form("%s_disp", Amc->GetName()));
      TH1D * dCor = corr ? shapeDisplay(corr, Form("%s_disp", corr->GetName())) : nullptr;
      styleHist(dA,   kBlack,   20);
      styleHist(dAmc, kAzure+2, 21);
      if (dCor) styleHist(dCor, kRed, 22);
      double ymax = std::max(dA->GetMaximum(), dAmc->GetMaximum());
      if (dCor) ymax = std::max(ymax, dCor->GetMaximum());

      double floor = insitu_utility::lowXjFloor(ir, ptlow);

      c->cd();
      double x0 = double(k)/nPtUsed, x1 = double(k+1)/nPtUsed;
      TPad * p1 = new TPad(Form("p1_%s", tag), "", x0, .32, x1, 1);
      TPad * p2 = new TPad(Form("p2_%s", tag), "", x0, 0, x1, .32);
      p1->SetLeftMargin(.15);
      p1->SetRightMargin(.04);
      p1->SetTopMargin(.05);
      p1->SetBottomMargin(.02);
      p2->SetLeftMargin(.15);
      p2->SetRightMargin(.04);
      p2->SetTopMargin(.03);
      p2->SetBottomMargin(.3);
      p1->Draw();
      p2->Draw();

      // ---- top: shape-normalized distributions ----
      p1->cd();
      p1->SetTicks();
      dAmc->SetMinimum(0);
      dAmc->SetMaximum(ymax*1.9);
      dAmc->GetYaxis()->SetTitle("Shape-normalized counts / bin width");
      dAmc->GetYaxis()->SetTitleSize(0.05);
      dAmc->GetYaxis()->SetTitleOffset(1.4);
      dAmc->GetYaxis()->SetLabelSize(0.045);
      dAmc->GetXaxis()->SetLabelSize(0);
      dAmc->Draw("p e");
      dA->Draw("p e same");
      if (dCor) dCor->Draw("p e same");
      TLine * lf = new TLine(floor, 0, floor, ymax*1.2);
      lf->SetLineStyle(2);
      lf->SetLineColor(kGray+2);
      lf->Draw("same");

      insitu_utility::drawSPhenixLabel({"p+p Run24 Data"}, {
          Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ptlow, pthigh),
          Form("Jet R=%.1f", ana::JetRs[ir]), "No in-situ JES applied"},
          .2, .88, 16, p1->GetWh()*p1->GetHNDC());

      TLegend * l = new TLegend(.56, .73, .95, .93);
      l->SetBorderSize(0);
      l->SetFillStyle(0);
      l->SetTextFont(43);
      l->SetTextSize(15);
      l->AddEntry(dA,   "Region A (Data)", "p");
      l->AddEntry(dAmc, "Region A (Pythia8)", "p");
      if (dCor) l->AddEntry(dCor, "Purity-corrected Data", "p");
      l->AddEntry(lf, Form("In-situ floor (x_{J#gamma} #geq %.2f)", floor), "l");
      l->Draw();

      double mA, eA, mMC, eMC, mCor = 0, eCor = 0;
      meanAboveFloor(A, floor, mA, eA);
      meanAboveFloor(Amc, floor, mMC, eMC);
      if (corr) meanAboveFloor(corr, floor, mCor, eCor);
      TLatex * tex = new TLatex();
      tex->SetNDC();
      tex->SetTextFont(43);
      tex->SetTextSize(15);
      tex->SetTextColor(kBlack);
      tex->DrawLatex(.56, .64, Form("#LTx_{J#gamma}#GT_{Data} = %.3f #pm %.3f", mA, eA));
      tex->SetTextColor(kAzure+2);
      tex->DrawLatex(.56, .58, Form("#LTx_{J#gamma}#GT_{Pythia8} = %.3f #pm %.3f", mMC, eMC));
      tex->SetTextColor(kRed);
      if (corr) tex->DrawLatex(.56, .52, Form("#LTx_{J#gamma}#GT_{corr.} = %.3f #pm %.3f", mCor, eCor));
      else      tex->DrawLatex(.56, .52, "Purity correction unavailable");
      tex->SetTextColor(kGray+2);
      tex->DrawLatex(.56, .46, "(means above the floor)");

      // ---- bottom: ratio to Pythia8 ----
      p2->cd();
      p2->SetTicks();
      TH1D * rA = (TH1D*)dA->Clone(Form("hratioA_%s", tag));
      rA->Divide(dAmc);
      sanitize(rA);
      rA->SetMinimum(0);
      rA->SetMaximum(2);
      rA->GetYaxis()->SetTitle("Data / Pythia8");
      rA->GetYaxis()->SetNdivisions(505);
      rA->GetYaxis()->SetTitleSize(0.1);
      rA->GetYaxis()->SetTitleOffset(0.65);
      rA->GetYaxis()->SetLabelSize(0.09);
      rA->GetXaxis()->SetTitle("x_{J#gamma}");
      rA->GetXaxis()->SetTitleSize(0.12);
      rA->GetXaxis()->SetTitleOffset(1.0);
      rA->GetXaxis()->SetLabelSize(0.1);
      rA->Draw("p e");
      if (dCor) {
        TH1D * rC = (TH1D*)dCor->Clone(Form("hratioCorr_%s", tag));
        rC->Divide(dAmc);
        sanitize(rC);
        rC->Draw("p e same");
        fout->cd();
        rC->Write();
      }
      TLine * one = new TLine(rA->GetXaxis()->GetXmin(), 1, rA->GetXaxis()->GetXmax(), 1);
      one->SetLineStyle(9);
      one->Draw("same");
      TLine * lf2 = new TLine(floor, 0, floor, 2);
      lf2->SetLineStyle(2);
      lf2->SetLineColor(kGray+2);
      lf2->Draw("same");
      fout->cd();
      rA->Write();
    }
    c->SaveAs(pdfPath.c_str());
  }

  c->SaveAs(Form("%s]", pdfPath.c_str()));
  fout->Close();
  cout << "Wrote " << pdfPath << endl;
  cout << "Wrote " << rootPath << endl;
}
