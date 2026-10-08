#include "../../src/ana.h"
#include "../../src/insitu_utility.h"
// Explicit load - see drawing/draw_final_result.C. Run interpreted, never with ACLiC "+".
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Page 3 of draw_gammajet_extrap.C drawn large, all modes normalized by the Default histogram's integral: the MC reference x_J (R=0.4, region A, x_J floor,
// cross section x in-situ weight) for each low-pT JER width extrapolation, photon-pT bins side by side on one wide page,
// with the ratio to the template ("Default") below; the no-extra-smearing reference is not drawn here. Inputs and mode definitions as in draw_gammajet_extrap.C.

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
}

void draw_gammajet_extrap_xj() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();

  float lowXj[ana::nPtBinsUsed];
  for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) lowXj[ipt] = insitu_utility::lowXjFloor(ir, ana::ptBinsUsed[ipt]);

  vector<TH1D*> hXj[nMode];
  float mean[nMode][ana::nPtBinsUsed], err[nMode][ana::nPtBinsUsed];
  for (int m = 0; m < nMode; m++) {
    hXj[m] = insitu_utility::buildMCXjByPtBin(samples(modes[m]), 0, ir, Form("hxj_%s", modes[m]), lowXj);
    insitu_utility::referenceMeans(samples(modes[m]), 0, ir, mean[m], err[m], lowXj);
  }
  // Every mode scaled by the Default (template) histogram's normalization, so Default is 1/N dN/dx and
  // the others keep their yield relative to it (same cross-section weights for all modes).
  for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) {
    double norm = hXj[1][ipt]->Integral();
    for (int m = 0; m < nMode; m++) {
      if (norm > 0) hXj[m][ipt]->Scale(1.0/norm, "width");
      sanitize(hXj[m][ipt]);
    }
  }

  string outdir = ana::path("claude_checks/smear_test/pdfs");
  gSystem->mkdir(outdir.c_str(), true);
  string pdf = outdir + "/gammajet_extrap_xj_R04.pdf";
  // One wide page: one column (1000 x 1000 each) per photon-pT bin.
  const int nCol = ana::nPtBinsUsed;
  TCanvas * c = new TCanvas("c", "", 1000*nCol, 1000);
  for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) {
    const double x0 = double(ipt)/nCol, x1 = double(ipt + 1)/nCol;
    c->cd();
    TPad * top = new TPad(Form("top%d", ipt), "", x0, 0.3, x1, 1);
    top->SetBottomMargin(0.02); top->SetLeftMargin(0.1); top->SetRightMargin(0.03); top->SetTopMargin(0.05);
    top->Draw();
    TPad * bot = new TPad(Form("bot%d", ipt), "", x0, 0, x1, 0.3);
    bot->SetTopMargin(0.03); bot->SetBottomMargin(0.32); bot->SetLeftMargin(0.1); bot->SetRightMargin(0.03);
    bot->Draw();

    top->cd();
    TH1D * href = hXj[1][ipt];
    double ymax = 0;
    for (int m = 1; m < nMode; m++) ymax = std::max(ymax, hXj[m][ipt]->GetMaximum());
    href->SetMaximum(1.45*ymax); href->SetMinimum(0);
    href->GetYaxis()->SetTitle("1/N_{Default} dN/dx_{J#gamma}");
    href->GetYaxis()->SetTitleSize(0.05); href->GetYaxis()->SetLabelSize(0.045); href->GetYaxis()->SetTitleOffset(0.9);
    href->GetXaxis()->SetLabelSize(0);
    TLegend * l = new TLegend(0.42, 0.62, 0.96, 0.92);
    l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.036);
    for (int m : {1, 2, 3}) {
      TH1D * h = hXj[m][ipt];
      h->SetLineColor(cols[m]); h->SetLineWidth(1); h->SetLineStyle(m == 1 ? 1 : (m == 0 ? 3 : 2));
      h->Draw(m == 1 ? "hist" : "hist same");
      l->AddEntry(h, Form("%s: #LTx_{J#gamma}#GT = %.4f", labels[m], mean[m][ipt]), "l");
    }
    l->Draw();
    insitu_utility::drawSPhenixLabel({"Pythia8 #gamma+jet MC"},
        {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", ana::ptBinsUsed[ipt], ana::ptBinsUsed[ipt+1]),
         "Jet R=0.4, Region A, reco level", Form("x_{J#gamma} floor %.2f", lowXj[ipt])}, .13, .88, 18, top->GetWh());

    bot->cd();
    TH1D * frame = (TH1D*)href->Clone(Form("frame%d", ipt)); frame->Reset();
    frame->GetYaxis()->SetTitle("ratio to Default"); frame->SetMinimum(0.8); frame->SetMaximum(1.2);
    frame->GetXaxis()->SetTitle("x_{J#gamma}");
    frame->GetXaxis()->SetLabelSize(0.1); frame->GetXaxis()->SetTitleSize(0.12); frame->GetXaxis()->SetTitleOffset(1.0);
    frame->GetYaxis()->SetLabelSize(0.09); frame->GetYaxis()->SetTitleSize(0.1); frame->GetYaxis()->SetTitleOffset(0.42);
    frame->GetYaxis()->SetNdivisions(505);
    frame->Draw("axis");
    TLine * one = new TLine(frame->GetXaxis()->GetXmin(), 1, frame->GetXaxis()->GetXmax(), 1);
    one->SetLineStyle(2); one->Draw();
    for (int m : {2, 3}) {
      TH1D * r = (TH1D*)hXj[m][ipt]->Clone(Form("r%d_%d", m, ipt));
      r->Divide(href);
      sanitize(r);
      r->SetLineColor(cols[m]); r->SetMarkerColor(cols[m]); r->SetMarkerStyle(20); r->SetMarkerSize(0.4); r->SetLineStyle(1); r->SetLineWidth(1);
      r->Draw("same e");
    }
  }
  c->SaveAs(pdf.c_str());
  printf("Wrote %s\n", pdf.c_str());
}
