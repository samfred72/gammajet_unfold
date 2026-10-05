#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/unfold_utility.h"
// The original gammajet project builds its OWN, differently-laid-out drawer/ana classes
// into /home/samson72/root/lib/libgammajet.so, sitting on the same library search path
// as this project's libgammajet_unfold.so. Without forcing which one loads first, ROOT's
// implicit symbol autoload can bind drawer/ana calls to the wrong (mismatched-layout)
// library and segfault - explicit load removes the ambiguity.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Diagnostic-only macro: compares the reco-level (pre-unfolding) purity-corrected xJ
// spectrum from two methods -
//   1) raw region A, uncorrected
//   2) the OLD, superseded production method (purityCorrectOld below), which assumes
//      region C is 100% background
//   3) the exact two-purity method (unfold_utility::purityCorrect), which additionally
//      uses region C's own purity P_C (from puritymaker.C's same leakage-corrected
//      bootstrap, via the ana::getPurityC getter) and solves the 2x2 linear system
//        A_i = S^A s_i + B^A bkg_i,  C_i = S^C s_i + B^C bkg_i
//      for the shared signal/background shapes s_i, bkg_i, instead of assuming region C's
//      raw shape IS the background shape. See src/unfold_utility.h for the derivation;
//      this reduces to method (2) exactly when P_C -> 0.
// No unfolding here - this only tests whether the correction method itself looks
// reasonable at reco level. This macro is what motivated switching production code
// (drawing/draw_purity_corrected.C and the insitu/ macros) over to the two-purity
// method - kept around as a standing sanity check / regression comparison.
const int ir = 2; // nominal jet radius index (R=0.4), matches draw_purity_corrected.C
const int nPtBinsUsed = ana::nPtBinsUsed;
// densityForDisplay now lives in unfold_utility - see src/unfold_utility.h.

// The OLD, superseded single-purity method - no longer used anywhere in production
// (draw_purity_corrected.C and the rest of the pipeline call unfold_utility::purityCorrect
// now), kept here only as the comparison baseline this macro exists to check against.
// Signal(xJ) = A(xJ) - (1-P_A)*(N_A/N_C)*C(xJ), i.e. assumes region C is pure background.
TH1D * purityCorrectOld(TH1D * A, TH1D * C, float p, float pErrLow, float pErrHigh, const char * name) {
  float NA = A->Integral();
  float NC = C->Integral();
  if (NC <= 0) {
    cout << "WARNING: " << name << " has zero region-C statistics - cannot cross-normalize, skipping." << endl;
    return nullptr;
  }
  float scale = (1-p)*(NA/NC);
  TH1D * h = (TH1D*)A->Clone(name);
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

// (1-P_A)*(N_A/N_C)*C(xJ): the piece subtracted from region A by the old method - same
// as bkgFromRegionC in drawing/draw_purity_corrected.C. TH1D::Scale() carries C's own
// bin errors through proportionally, so no separate error handling is needed here.
TH1D * bkgOld(TH1D * C, float p, float NA, float NC, const char * name) {
  float scale = (1-p)*(NA/NC);
  TH1D * h = (TH1D*)C->Clone(name);
  h->Scale(scale);
  return h;
}

// The exact two-purity method itself (signal and background) now lives in
// unfold_utility::purityCorrect/purityCorrectBkg, shared with drawing/draw_purity_corrected.C
// and the insitu/ macros - see src/unfold_utility.h for the derivation. Only the OLD,
// deliberately-superseded single-purity formula (purityCorrectOld/bkgOld) stays local
// here, since this macro's whole point is comparing the two.

void draw_purity_method_compare(string systag = "nominal") {
  gStyle->SetOptStat(0);

  drawer d("pythia", systag);
  string rootPath = Form("%s/hists/purity_method_compare_%s.root", ana::dir(), systag.c_str());
  string pdfPath  = Form("%s/pdfs/purity_method_compare_%s.pdf", ana::dir(), systag.c_str());

  TFile * fout = TFile::Open(rootPath.c_str(), "RECREATE");
  TCanvas * c = new TCanvas("c","",700,700);
  c->SaveAs(Form("%s[", pdfPath.c_str()));

  TH1D * flatA = d.get(Form("hrecoxj%i_0",ir), 0);
  TH1D * flatC = d.get(Form("hrecoxj%i_2",ir), 0);

  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    float ptlow  = ana::ptBins[ipt];
    float pthigh = ana::ptBins[ipt+1];
    float pA        = ana::getPurity(ptlow, pthigh, systag);
    float pAErrLow  = ana::getPurityErrorLow(ptlow, pthigh, systag);
    float pAErrHigh = ana::getPurityErrorHigh(ptlow, pthigh, systag);
    float pC        = ana::getPurityC(ptlow, pthigh, systag);
    float pCErrLow  = ana::getPurityCErrorLow(ptlow, pthigh, systag);
    float pCErrHigh = ana::getPurityCErrorHigh(ptlow, pthigh, systag);

    TH1D * A = unfold_utility::unflattenXj(flatA, ipt, Form("hxjA_pt%d", ipt));
    TH1D * C = unfold_utility::unflattenXj(flatC, ipt, Form("hxjC_pt%d", ipt));

    TH1D * hOld = purityCorrectOld(A, C, pA, pAErrLow, pAErrHigh, Form("hxjOld_pt%d", ipt));
    TH1D * hNew = unfold_utility::purityCorrect(A, C, pA, pAErrLow, pAErrHigh, pC, pCErrLow, pCErrHigh, Form("hxjNew_pt%d", ipt));

    float NA = A->Integral();
    float NC = C->Integral();
    TH1D * hBkgOld = bkgOld(C, pA, NA, NC, Form("hxjBkgOld_pt%d", ipt));
    TH1D * hBkgNew = unfold_utility::purityCorrectBkg(A, C, pA, pC, Form("hxjBkgNew_pt%d", ipt));

    // Page 1: region A (raw) with the two competing background estimates overlaid -
    // same layout as drawing/draw_purity_corrected.C's page 1 (Region A / (1-P)*N_A/N_C*C),
    // just with both the old and new background curves shown together for comparison.
    c->Clear();
    c->cd();
    gPad->SetTicks();
    gPad->SetLeftMargin(.15);
    TH1D * AdispBkgPage = unfold_utility::densityForDisplay(A, Form("hxjA_pt%d_bkgpage_disp", ipt));
    AdispBkgPage->SetLineColor(kBlack);
    AdispBkgPage->SetLineWidth(2);
    AdispBkgPage->GetXaxis()->SetTitle("x_{J#gamma}");
    AdispBkgPage->Draw("hist");
    TLegend * lbkg = new TLegend(.55,.35,.85,.48);
    lbkg->SetLineWidth(0);
    lbkg->SetTextSize(0.022);
    lbkg->AddEntry(AdispBkgPage, "Region A (raw)");
    double ymaxBkg = AdispBkgPage->GetMaximum();
    TH1D * hBkgOlddisp = nullptr;
    TH1D * hBkgNewdisp = nullptr;
    if (hBkgOld) {
      hBkgOlddisp = unfold_utility::densityForDisplay(hBkgOld, Form("hxjBkgOld_pt%d_disp", ipt));
      hBkgOlddisp->SetLineColor(kAzure+2);
      hBkgOlddisp->SetLineWidth(2);
      ymaxBkg = std::max(ymaxBkg, hBkgOlddisp->GetMaximum());
      lbkg->AddEntry(hBkgOlddisp, "Bkg, old: (1-P_{A})#times#frac{N_{A}}{N_{C}}#times C");
    }
    if (hBkgNew) {
      hBkgNewdisp = unfold_utility::densityForDisplay(hBkgNew, Form("hxjBkgNew_pt%d_disp", ipt));
      hBkgNewdisp->SetLineColor(kOrange+7);
      hBkgNewdisp->SetLineWidth(2);
      ymaxBkg = std::max(ymaxBkg, hBkgNewdisp->GetMaximum());
      lbkg->AddEntry(hBkgNewdisp, "Bkg, new: solved from P_{A} & P_{C}");
    }
    AdispBkgPage->GetYaxis()->SetRangeUser(0, ymaxBkg*1.4);
    if (hBkgOlddisp) hBkgOlddisp->Draw("hist e same");
    if (hBkgNewdisp) hBkgNewdisp->Draw("hist e same");
    lbkg->Draw();
    d.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ptlow,pthigh),
        Form("P_{A} = %.3f +%.3f/-%.3f", pA, pAErrHigh, pAErrLow),
        Form("P_{C} = %.3f +%.3f/-%.3f", pC, pCErrHigh, pCErrLow),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir])}, .5, .85, 15, 700);
    c->SaveAs(pdfPath.c_str());
    fout->cd();
    if (hBkgOld) hBkgOld->Write();
    if (hBkgNew) hBkgNew->Write();

    // Page 2: the raw/old-corrected/new-corrected signal comparison, with a ratio panel.
    c->Clear();
    c->cd();
    TPad * p1 = new TPad(Form("pcmp1_%d",ipt),"",0,.35,1,1);
    TPad * p2 = new TPad(Form("pcmp2_%d",ipt),"",0,0,1,.35);
    p1->Draw();
    p2->Draw();

    p1->cd();
    p1->SetBottomMargin(0.02);
    p1->SetLeftMargin(.15);
    gPad->SetTicks(1,1);

    TH1D * Adisp = unfold_utility::densityForDisplay(A, Form("hxjA_pt%d_disp", ipt));
    Adisp->SetLineColor(kBlack);
    Adisp->SetMarkerColor(kBlack);
    Adisp->SetMarkerStyle(20);
    Adisp->SetLineWidth(2);
    Adisp->GetXaxis()->SetLabelSize(0);
    Adisp->GetXaxis()->SetTitle("");
    Adisp->GetYaxis()->SetTitle("Counts / bin width");

    TLegend * l = new TLegend(.55,.32,.85,.45);
    l->SetLineWidth(0);
    l->SetTextSize(0.024);
    l->AddEntry(Adisp, "Region A (raw)");

    double ymax = Adisp->GetMaximum();
    TH1D * Olddisp = nullptr;
    TH1D * Newdisp = nullptr;
    if (hOld) {
      Olddisp = unfold_utility::densityForDisplay(hOld, Form("hxjOld_pt%d_disp", ipt));
      Olddisp->SetLineColor(kAzure+2);
      Olddisp->SetMarkerColor(kAzure+2);
      Olddisp->SetMarkerStyle(21);
      Olddisp->SetLineWidth(2);
      ymax = std::max(ymax, Olddisp->GetMaximum());
      l->AddEntry(Olddisp, "Old method (P_{A} only)");
    }
    if (hNew) {
      Newdisp = unfold_utility::densityForDisplay(hNew, Form("hxjNew_pt%d_disp", ipt));
      Newdisp->SetLineColor(kOrange+7);
      Newdisp->SetMarkerColor(kOrange+7);
      Newdisp->SetMarkerStyle(22);
      Newdisp->SetLineWidth(2);
      ymax = std::max(ymax, Newdisp->GetMaximum());
      l->AddEntry(Newdisp, "New method (P_{A} & P_{C})");
    }
    Adisp->GetYaxis()->SetRangeUser(0, ymax*1.4);
    Adisp->Draw("p e");
    if (Olddisp) Olddisp->Draw("p e same");
    if (Newdisp) Newdisp->Draw("p e same");
    l->Draw();

    d.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ptlow,pthigh),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir]),
        Form("P_{A} = %.3f +%.3f/-%.3f", pA, pAErrHigh, pAErrLow),
        Form("P_{C} = %.3f +%.3f/-%.3f", pC, pCErrHigh, pCErrLow)}, .5, .85, 15, gPad->GetWh()*0.8);

    // Means of the (shape, not signal-count) background-subtracted xJ distributions -
    // same "Mean raw/corr" convention as drawing/draw_purity_corrected.C, extended to
    // both correction methods. Drawn below the legend so it doesn't collide with it.
    d.drawText(Form("Mean raw: %.2f #pm %.2f", Adisp->GetMean(), Adisp->GetMeanError()), .6, .28);
    if (Olddisp) d.drawText(Form("Mean old: %.2f #pm %.2f", Olddisp->GetMean(), Olddisp->GetMeanError()), .6, .23);
    if (Newdisp) d.drawText(Form("Mean new: %.2f #pm %.2f", Newdisp->GetMean(), Newdisp->GetMeanError()), .6, .18);

    p2->cd();
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.3);
    p2->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    if (Olddisp && Newdisp) {
      TH1D * hratio = (TH1D*)Newdisp->Clone(Form("hxjratio_pt%d", ipt));
      hratio->Divide(Olddisp);
      hratio->SetLineColor(kBlack);
      hratio->SetMarkerColor(kBlack);
      hratio->SetMarkerStyle(20);
      hratio->GetYaxis()->SetRangeUser(0.5,1.5);
      hratio->GetYaxis()->SetTitle("New / Old");
      hratio->GetYaxis()->SetTitleSize(0.09);
      hratio->GetYaxis()->SetLabelSize(0.08);
      hratio->GetXaxis()->SetTitle("x_{J#gamma}");
      hratio->GetXaxis()->SetTitleSize(0.09);
      hratio->GetXaxis()->SetLabelSize(0.08);
      hratio->Draw("p e");
      TLine * line = new TLine(ana::unfoldXjBins[0],1,ana::unfoldXjBins[ana::nUnfoldXjBins],1);
      line->SetLineStyle(9);
      line->Draw("same");
      fout->cd();
      hratio->Write();
    }
    c->SaveAs(pdfPath.c_str());

    fout->cd();
    A->Write();
    C->Write();
    if (hOld) hOld->Write();
    if (hNew) hNew->Write();
  }

  c->SaveAs(Form("%s]", pdfPath.c_str()));
  fout->Close();
  cout << "Wrote " << pdfPath << endl;
}
