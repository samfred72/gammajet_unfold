#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/insitu_utility.h"
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
using namespace std;

R__LOAD_LIBRARY(libgammajet_unfold.so);

// Diagnostic for the spike (pa~0.911-0.928) and cliff (pa~0.97) in
// grid_insitu_shapechi2.C's Region-A-only shape-chi2 curve. Rebuilds exactly the same
// per-(ptbin,xjbin) bin-fraction chi2 pieces grid_insitu_shapechi2.C's scan loop computes
// (see that file's referenceShape()/scan-loop comments for the formula/rationale this
// copies), but frozen at three fixed trial pa values instead of scanned, so the actual
// Data-vs-MC xJ shape and the per-bin chi2 pull can be looked at directly instead of only
// the summed curve.
// insitu/ is split into inputs/ (the raw Data/Photon insitu ntuples, written by
// unfolder.h's production pipeline) and pdfs/ (this debug macro's own .pdf output,
// no .root output).
const char * insitu_input_dir = "/home/samson72/sphnx/gammajet_unfold/insitu/inputs";
const char * insitu_pdf_dir   = "/home/samson72/sphnx/gammajet_unfold/insitu/pdfs";
const int nPtBinsUsed = ana::nPtBinsUsed;
// Same explicit low-stat-tail exclusion as grid_insitu_shapechi2.C (see that file's
// header comment / gammajet_unfold/CLAUDE.md's bin-selection ground rule).
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;

map<int,double> photon_scale = {{5,146359.3},{10,6944.675},{20,130.4461}};

// Fixed MC reference xJ shape (bin fraction) and its per-bin error - identical
// computation to grid_insitu_shapechi2.C's referenceShape().
void referenceShape(const vector<pair<string,double>> & samples, int abcdSelect, int ir,
    vector<vector<double>> & refFrac, vector<vector<double>> & refFracErr, const float lowXj[]) {
  vector<vector<double>> sumw(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
  vector<vector<double>> sumw2(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
  for (auto & s : samples) {
    TFile * f = TFile::Open(s.first.c_str(), "READ");
    if (!f || f->IsZombie()) { cout << "WARNING: could not open " << s.first << endl; continue; }
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
      if (ipt < ana::firstUsedPtBin || ipt >= ana::firstUsedPtBin + nPtBinsUsed) continue;
      ipt -= ana::firstUsedPtBin;
      double x = jet_pt/pho_pt;
      if (x < lowXj[ipt]) continue;
      int ixj = ana::findUnfoldXjBin(x);
      if (ixj < 0 || ixj >= ana::nUnfoldXjBins) continue;
      double w = s.second*mcWeight;
      sumw[ipt][ixj]  += w;
      sumw2[ipt][ixj] += w*w;
    }
    f->Close();
  }
  refFrac.assign(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
  refFracErr.assign(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    double N = 0;
    for (double w : sumw[ipt]) N += w;
    if (N <= 0) continue;
    for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
      refFrac[ipt][ixj]    = sumw[ipt][ixj]/N;
      refFracErr[ipt][ixj] = sqrt(sumw2[ipt][ixj])/N;
    }
  }
}

struct BinPull { int ipt, ixj; double fData, fMC, errt, pull, chi2; };

void debug_shapechi2_spike(string systag = "nominal", int ir = 2) {
  gStyle->SetOptStat(0);

  string dataFile = insitu_utility::insituFilename(insitu_input_dir, "Data", "", systag);
  vector<DataEvent> dataA = insitu_utility::cacheDataEvents(dataFile.c_str(), 0, ir);
  cout << "Cached Data Region A events: " << dataA.size() << endl;

  // Low-xJ floor per used pT bin - same cut grid_insitu_shapechi2.C now applies (see
  // src/insitu_utility.h's lowXjFloor comment).
  float lowXj[nPtBinsUsed];
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) lowXj[ipt] = insitu_utility::lowXjFloor(ir, ana::ptBinsUsed[ipt]);

  vector<pair<string,double>> mcSamples = {
    {insitu_utility::insituFilename(insitu_input_dir, "Photon5",  "pythia", systag), photon_scale[5]},
    {insitu_utility::insituFilename(insitu_input_dir, "Photon10", "pythia", systag), photon_scale[10]},
    {insitu_utility::insituFilename(insitu_input_dir, "Photon20", "pythia", systag), photon_scale[20]},
  };
  vector<vector<double>> refFrac, refFracErr;
  referenceShape(mcSamples, 0, ir, refFrac, refFracErr, lowXj);

  vector<float> paValues = {0.92, 0.96, 0.98};
  const char * pdfPath = Form("%s/debug_shapechi2_spike_%s_%s.pdf", insitu_pdf_dir, systag.c_str(), ana::rnames[ir]);
  TCanvas * c = new TCanvas("c", "", 1500, 700);
  c->SaveAs(Form("%s[", pdfPath));

  for (float pa : paValues) {
    // Per-(ptbin,xjbin) raw Data counts at this pa - identical to the inner loop of
    // grid_insitu_shapechi2.C's scan, just at one fixed pa instead of every grid point.
    vector<vector<double>> countA(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
    for (auto & ev : dataA) {
      float x = (ev.jet_pt/pa)/ev.pho_pt;
      if (x < lowXj[ev.ptbin]) continue;
      int ixj = ana::findUnfoldXjBin(x);
      if (ixj < 0 || ixj >= ana::nUnfoldXjBins) continue;
      countA[ev.ptbin][ixj] += 1;
    }

    vector<BinPull> pulls;
    double chisqTotal = 0;
    vector<TH1D*> hData(nPtBinsUsed), hMC(nPtBinsUsed);
    for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
      double NA = 0;
      for (double n : countA[ipt]) NA += n;

      hData[ipt] = new TH1D(Form("hData_pa%.2f_pt%d", pa, ipt), ";x_{J#gamma};Bin fraction",
          ana::nUnfoldXjBins, ana::unfoldXjBins);
      hMC[ipt]   = new TH1D(Form("hMC_pa%.2f_pt%d", pa, ipt), ";x_{J#gamma};Bin fraction",
          ana::nUnfoldXjBins, ana::unfoldXjBins);

      for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
        double fData = NA > 0 ? countA[ipt][ixj]/NA : 0;
        double errData = NA > 0 ? sqrt(countA[ipt][ixj])/NA : 0;
        hData[ipt]->SetBinContent(ixj+1, fData);
        hData[ipt]->SetBinError(ixj+1, errData);
        hMC[ipt]->SetBinContent(ixj+1, refFrac[ipt][ixj]);
        hMC[ipt]->SetBinError(ixj+1, refFracErr[ipt][ixj]);

        if (NA <= 0 || ixj >= nXjBinsForChi2) continue;
        double errt = sqrt(errData*errData + refFracErr[ipt][ixj]*refFracErr[ipt][ixj]);
        errt = std::max(errt, 1.0/NA); // same errt floor as grid_insitu_shapechi2.C
        double diff = fData - refFrac[ipt][ixj];
        double chi2bin = diff*diff/(errt*errt);
        chisqTotal += chi2bin;
        pulls.push_back({ipt, ixj, fData, refFrac[ipt][ixj], errt, diff/errt, chi2bin});
      }
    }

    sort(pulls.begin(), pulls.end(), [](const BinPull & a, const BinPull & b) { return a.chi2 > b.chi2; });
    cout << "\n=== pa=" << pa << "  chi2A(total, first " << nXjBinsForChi2 << " xJ bins)="
         << chisqTotal << " ===" << endl;
    cout << "Top contributing (ptbin,xjbin) - ptbin edges " << ana::ptBinsUsed[0] << ".."
         << ana::ptBinsUsed[nPtBinsUsed] << " GeV, xJ edges below:" << endl;
    for (int k = 0; k < (int)pulls.size() && k < 6; k++) {
      const BinPull & p = pulls[k];
      cout << "  ptbin=" << p.ipt << " [" << ana::ptBinsUsed[p.ipt] << "," << ana::ptBinsUsed[p.ipt+1]
           << ") xjbin=" << p.ixj << " [" << ana::unfoldXjBins[p.ixj] << "," << ana::unfoldXjBins[p.ixj+1]
           << ")  fData=" << p.fData << " fMC=" << p.fMC << " errt=" << p.errt
           << " pull=" << p.pull << " chi2=" << p.chi2 << endl;
    }

    // ---- Page: nPtBinsUsed columns, each with a fraction-comparison pad on top and a
    // signed-pull pad on the bottom (pull = (fData-fMC)/errt, the actual per-bin
    // chi2-sum term's sign+magnitude) - the vertical dashed line marks
    // ana::unfoldXjBins[nXjBinsForChi2], the edge beyond which bins are dropped from
    // the chi2 sum (the low-stat tail, same 3 bins every chi2 in this project drops).
    c->Clear();
    TLatex pageTitle;
    pageTitle.SetNDC();
    pageTitle.SetTextSize(.03);
    pageTitle.SetTextFont(62);
    pageTitle.DrawLatex(.4, .96, Form("Region A shape check, p_{a}=%.2f, systag=%s", pa, systag.c_str()));

    double xcut = ana::unfoldXjBins[nXjBinsForChi2];
    for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
      // Each column's pads must attach to the canvas itself, not to whichever pad
      // happened to be current (pBot, from the previous column's cd()) - without this
      // cd() back to c, TPad::Draw() nests the new pads inside the last column's pBot
      // instead of placing them side by side.
      c->cd();
      double x0 = ipt/(double)nPtBinsUsed, x1 = (ipt+1)/(double)nPtBinsUsed;
      TPad * pTop = new TPad(Form("pTop%d",ipt), "", x0, .32, x1, .92);
      TPad * pBot = new TPad(Form("pBot%d",ipt), "", x0, .05, x1, .32);
      pTop->SetLeftMargin(.18); pTop->SetBottomMargin(.02); pTop->SetRightMargin(.03);
      pBot->SetLeftMargin(.18); pBot->SetTopMargin(.02); pBot->SetBottomMargin(.35); pBot->SetRightMargin(.03);
      pTop->Draw(); pBot->Draw();

      pTop->cd();
      double ymax = std::max(hMC[ipt]->GetMaximum(), hData[ipt]->GetMaximum())*1.4;
      if (ymax <= 0) ymax = 1;
      hMC[ipt]->SetLineColor(kMagenta+1);
      hMC[ipt]->SetLineWidth(2);
      hMC[ipt]->GetYaxis()->SetRangeUser(0, ymax);
      hMC[ipt]->GetXaxis()->SetLabelSize(0);
      hMC[ipt]->GetYaxis()->SetTitleSize(.06);
      hMC[ipt]->GetYaxis()->SetLabelSize(.05);
      hMC[ipt]->Draw("hist");
      hData[ipt]->SetLineColor(kBlue);
      hData[ipt]->SetMarkerColor(kBlue);
      hData[ipt]->SetMarkerStyle(20);
      hData[ipt]->Draw("p e same");
      TLine * vcut1 = new TLine(xcut, 0, xcut, ymax);
      vcut1->SetLineStyle(2);
      vcut1->SetLineColor(kGray+2);
      vcut1->Draw();
      if (ipt == 0) {
        TLegend * l = new TLegend(.35, .7, .95, .9);
        l->SetLineWidth(0);
        l->SetTextSize(.05);
        l->AddEntry(hMC[ipt], "Pythia8 reference", "l");
        l->AddEntry(hData[ipt], "Data (Region A)", "p");
        l->Draw();
      }
      TLatex ptlabel;
      ptlabel.SetNDC();
      ptlabel.SetTextSize(.06);
      ptlabel.DrawLatex(.22, .3, Form("%.0f<p_{T}^{#gamma}<%.0f GeV", ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1]));

      pBot->cd();
      TH1D * hPull = new TH1D(Form("hPull_pa%.2f_pt%d", pa, ipt), ";x_{J#gamma};pull", ana::nUnfoldXjBins, ana::unfoldXjBins);
      for (auto & p : pulls) {
        if (p.ipt != ipt) continue;
        hPull->SetBinContent(p.ixj+1, p.pull);
      }
      double pullmax = std::max(3.0, hPull->GetMaximum()*1.2);
      double pullmin = std::min(-3.0, hPull->GetMinimum()*1.2);
      hPull->SetLineColor(kBlack);
      hPull->SetFillColor(kAzure-4);
      hPull->GetYaxis()->SetRangeUser(pullmin, pullmax);
      hPull->GetYaxis()->SetTitleSize(.13);
      hPull->GetYaxis()->SetLabelSize(.11);
      hPull->GetYaxis()->SetTitleOffset(.5);
      hPull->GetXaxis()->SetTitleSize(.13);
      hPull->GetXaxis()->SetLabelSize(.11);
      hPull->GetXaxis()->SetTitleOffset(1.1);
      hPull->Draw("hist");
      TLine * vcut2 = new TLine(xcut, pullmin, xcut, pullmax);
      vcut2->SetLineStyle(2);
      vcut2->SetLineColor(kGray+2);
      vcut2->Draw();
      TLine * zero = new TLine(ana::unfoldXjBins[0], 0, ana::unfoldXjBins[ana::nUnfoldXjBins], 0);
      zero->Draw();
    }

    c->SaveAs(pdfPath);
  }

  c->SaveAs(Form("%s]", pdfPath));
  cout << "\nWrote " << pdfPath << endl;
}
