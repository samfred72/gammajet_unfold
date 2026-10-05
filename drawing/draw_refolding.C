#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Refolding closure: fold the nominal unfolded result back through the response
// (ApplyToTruth) and compare with the purity-corrected measurement, on the absolute scale.
// ApplyToTruth's output has no errors, so chi2 uses the measured error only. Finite iterations
// and prior mismatch mean closure is not exact; residuals comparable to the statistical error
// are the red flag.
//
// ApplyToTruth does not add fakes. Hfakes() is on the MC training scale, so the Data-scale fakes
// are (Hfakes/Hmeasured) x Data, as RooUnfoldBayes uses them.
//
// Also a pure-MC check: Pythia truth folded forward (+ training fakes) vs Pythia reco.

const int ir = 2; // nominal jet radius index (R=0.4)
const int nPtBinsUsed = ana::nPtBinsUsed; // reported bins start at ana::firstUsedPtBin
const int niterate = 2; // nominal iteration count
const vector<int> iterationsToScan = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; // for the bonus niter-dependence page

// The last 3 xJ bins per pT bin are low-count noise: excluded from chi2 only (still drawn).
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;

// chi2/NDF of refolded vs measured, measured error only.
double computeChi2NDF(TH1D * hRefolded, TH1D * hMeasured) {
  double chi2 = 0;
  int ndf = 0;
  for (int b = 1; b <= nXjBinsForChi2; b++) {
    double vR = hRefolded->GetBinContent(b);
    double vM = hMeasured->GetBinContent(b);
    double eM = hMeasured->GetBinError(b);
    if (eM <= 0) continue;
    chi2 += pow(vR-vM,2)/(eM*eM);
    ndf++;
  }
  return (ndf > 0) ? chi2/ndf : 0;
}

void draw_refolding(string systag = "nominal") {
  gStyle->SetOptStat(0);

  drawer d("pythia", systag);
  string pdfPath  = Form("%s/pdfs/draw_refolding_%s.pdf", ana::dir(), systag.c_str());
  string rootPath = Form("%s/hists/refolding_%s.root", ana::dir(), systag.c_str());

  // Response built explicitly so the same object serves ApplyToTruth.
  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i",ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i",ir), 1);
  RooUnfoldResponse * response = new RooUnfoldResponse(respRecoTemplate, respTruthTemplate, respMatrix2D);

  // Training-scale fakes.
  TH1D * flatFakesTraining = (TH1D*)response->Hfakes()->Clone("hFakesTrainingFlat");

  // Pure-MC check: truth folded forward plus training fakes.
  TH1D * flatRefoldedTruth = (TH1D*)response->ApplyToTruth(respTruthTemplate, "hRefoldedTruth");
  flatRefoldedTruth->Add(flatFakesTraining);

  TH1D * flatA = d.get(Form("hrecoxj%i_0",ir), 0);
  TH1D * flatC = d.get(Form("hrecoxj%i_2",ir), 0);
  TH1D * flatMeasured = unfold_utility::buildFullyCorrected(flatA, flatC, "data", systag);

  // Data-scale fakes: training fakes fraction x Data.
  TH1D * flatFakes = (TH1D*)flatMeasured->Clone("hFakesInDataFlat");
  flatFakes->Reset("ICES");
  for (int b = 1; b <= flatFakes->GetNbinsX(); b++) {
    double mesTraining = respRecoTemplate->GetBinContent(b);
    if (mesTraining <= 0) continue;
    double fakeFrac = flatFakesTraining->GetBinContent(b) / mesTraining;
    flatFakes->SetBinContent(b, fakeFrac * flatMeasured->GetBinContent(b));
  }

  // Refolded prediction including fakes.
  TH1D * flatUnfolded = unfold_utility::unfoldOnce(response, flatMeasured, niterate, "hUnfoldedNominal");
  TH1D * flatRefolded = (TH1D*)response->ApplyToTruth(flatUnfolded, "hRefoldedNominal");
  flatRefolded->Add(flatFakes);

  TFile * fout = TFile::Open(rootPath.c_str(), "RECREATE");
  TCanvas * c = new TCanvas("c","",700,900);
  c->SaveAs(Form("%s[", pdfPath.c_str()));

  cout << "Refolding closure (niter=" << niterate << "): pT bin, chi2/NDF (first "
       << nXjBinsForChi2 << " of " << ana::nUnfoldXjBins << " xJ bins), fakes integral" << endl;

  // Pages 1..nPtBinsUsed: measured vs refolded (absolute densities) with a ratio panel.
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hMeasured = unfold_utility::unflattenXj(flatMeasured, ipt, Form("hMeasured_pt%d", ipt));
    TH1D * hRefolded = unfold_utility::unflattenXj(flatRefolded, ipt, Form("hRefolded_pt%d", ipt));
    TH1D * hFakes    = unfold_utility::unflattenXj(flatFakes, ipt, Form("hFakes_pt%d", ipt));
    double chi2ndf = computeChi2NDF(hRefolded, hMeasured);
    cout << "  pt" << ipt << " (" << ana::ptBins[ipt] << "-" << ana::ptBins[ipt+1] << " GeV): chi2/NDF = " << chi2ndf
         << ", fakes integral = " << hFakes->Integral() << endl;

    TH1D * hMeasuredDisp = unfold_utility::densityForDisplay(hMeasured, Form("hMeasuredDisp_pt%d", ipt));
    TH1D * hRefoldedDisp = unfold_utility::densityForDisplay(hRefolded, Form("hRefoldedDisp_pt%d", ipt));

    c->Clear();
    c->cd();
    TPad * p1 = new TPad(Form("p1_%d",ipt),"",0,.35,1,1);
    TPad * p2 = new TPad(Form("p2_%d",ipt),"",0,0,1,.35);
    p1->Draw();
    p2->Draw();

    p1->cd();
    p1->SetBottomMargin(0.02);
    p1->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    hMeasuredDisp->SetLineColor(kBlack);
    hMeasuredDisp->SetMarkerColor(kBlack);
    hMeasuredDisp->SetMarkerStyle(20);
    hMeasuredDisp->SetLineWidth(2);
    hMeasuredDisp->GetXaxis()->SetLabelSize(0);
    hMeasuredDisp->GetXaxis()->SetTitle("");
    hMeasuredDisp->GetYaxis()->SetRangeUser(0, std::max(hMeasuredDisp->GetMaximum(), hRefoldedDisp->GetMaximum())*1.4);
    hMeasuredDisp->Draw("p e");
    // Line, not points: the refolded prediction has no per-bin error.
    hRefoldedDisp->SetLineColor(kRed);
    hRefoldedDisp->SetLineWidth(3);
    hRefoldedDisp->Draw("hist same");
    TLegend * l = new TLegend(.5,.65,.85,.80);
    l->SetLineWidth(0);
    l->SetTextSize(0.032);
    l->AddEntry(hMeasuredDisp, "Measured (purity-corrected)", "lp");
    l->AddEntry(hRefoldedDisp, Form("Refolded + fakes (%d iter.)", niterate), "l");
    l->Draw();
    d.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV, #chi^{2}/NDF = %.2f", ana::JetRs[ir], ana::jet_calib_pt_cut[ir], chi2ndf)}, .18, .85, 14, gPad->GetWh()*0.8);

    p2->cd();
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.3);
    p2->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    TH1D * hratio = (TH1D*)hRefoldedDisp->Clone(Form("hratio_pt%d", ipt));
    hratio->Divide(hMeasuredDisp);
    hratio->SetLineColor(kBlack);
    hratio->SetMarkerColor(kBlack);
    hratio->SetMarkerStyle(20);
    hratio->GetYaxis()->SetRangeUser(0.5,1.5);
    hratio->GetYaxis()->SetTitle("Refolded / Measured");
    hratio->GetYaxis()->SetTitleSize(0.09);
    hratio->GetYaxis()->SetTitleOffset(0.7);
    hratio->GetYaxis()->SetLabelSize(0.08);
    hratio->GetXaxis()->SetTitle("x_{J#gamma}");
    hratio->GetXaxis()->SetTitleSize(0.09);
    hratio->GetXaxis()->SetLabelSize(0.08);
    hratio->Draw("p e");
    TLine * line = new TLine(ana::unfoldXjBins[0],1,ana::unfoldXjBins[ana::nUnfoldXjBins],1);
    line->SetLineStyle(9);
    line->Draw("same");
    c->SaveAs(pdfPath.c_str());

    fout->cd();
    hMeasured->Write();
    hRefolded->Write();
    hFakes->Write();
    hratio->Write();
    delete hMeasured; delete hRefolded; delete hFakes;
  }

  // Pages (nPtBinsUsed+1)..: pure-MC check (no niter dependence).
  cout << "Pure-MC truth-refolding closure: pT bin, chi2/NDF (first "
       << nXjBinsForChi2 << " of " << ana::nUnfoldXjBins << " xJ bins)" << endl;
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hPythiaReco      = unfold_utility::unflattenXj(respRecoTemplate, ipt, Form("hPythiaReco_pt%d", ipt));
    TH1D * hRefoldedTruthPt = unfold_utility::unflattenXj(flatRefoldedTruth, ipt, Form("hRefoldedTruth_pt%d", ipt));
    double chi2ndf = computeChi2NDF(hRefoldedTruthPt, hPythiaReco);
    cout << "  pt" << ipt << " (" << ana::ptBins[ipt] << "-" << ana::ptBins[ipt+1] << " GeV): chi2/NDF = " << chi2ndf << endl;

    TH1D * hPythiaRecoDisp      = unfold_utility::densityForDisplay(hPythiaReco, Form("hPythiaRecoDisp_pt%d", ipt));
    TH1D * hRefoldedTruthDisp   = unfold_utility::densityForDisplay(hRefoldedTruthPt, Form("hRefoldedTruthDisp_pt%d", ipt));

    c->Clear();
    c->cd();
    TPad * p1 = new TPad(Form("pt1_%d",ipt),"",0,.35,1,1);
    TPad * p2 = new TPad(Form("pt2_%d",ipt),"",0,0,1,.35);
    p1->Draw();
    p2->Draw();

    p1->cd();
    p1->SetBottomMargin(0.02);
    p1->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    hPythiaRecoDisp->SetLineColor(kBlack);
    hPythiaRecoDisp->SetMarkerColor(kBlack);
    hPythiaRecoDisp->SetMarkerStyle(20);
    hPythiaRecoDisp->SetLineWidth(2);
    hPythiaRecoDisp->GetXaxis()->SetLabelSize(0);
    hPythiaRecoDisp->GetXaxis()->SetTitle("");
    hPythiaRecoDisp->GetYaxis()->SetRangeUser(0, std::max(hPythiaRecoDisp->GetMaximum(), hRefoldedTruthDisp->GetMaximum())*1.4);
    hPythiaRecoDisp->Draw("p e");
    hRefoldedTruthDisp->SetLineColor(kRed);
    hRefoldedTruthDisp->SetLineWidth(3);
    hRefoldedTruthDisp->Draw("hist same");
    TLegend * l = new TLegend(.5,.65,.85,.80);
    l->SetLineWidth(0);
    l->SetTextSize(0.032);
    l->AddEntry(hPythiaRecoDisp, "Pythia8 reco (training)", "lp");
    l->AddEntry(hRefoldedTruthDisp, "Refolded truth + fakes", "l");
    l->Draw();
    d.drawAll({"Pythia8 #gamma+jet MC"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV, #chi^{2}/NDF = %.2f", ana::JetRs[ir], ana::jet_calib_pt_cut[ir], chi2ndf)}, .18, .85, 14, gPad->GetWh()*0.8);

    p2->cd();
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.3);
    p2->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    TH1D * hratio = (TH1D*)hRefoldedTruthDisp->Clone(Form("hratio_truth_pt%d", ipt));
    hratio->Divide(hPythiaRecoDisp);
    hratio->SetLineColor(kBlack);
    hratio->SetMarkerColor(kBlack);
    hratio->SetMarkerStyle(20);
    hratio->GetYaxis()->SetRangeUser(0.5,1.5);
    hratio->GetYaxis()->SetTitle("Refolded truth / Pythia reco");
    hratio->GetYaxis()->SetTitleSize(0.09);
    hratio->GetYaxis()->SetTitleOffset(0.7);
    hratio->GetYaxis()->SetLabelSize(0.08);
    hratio->GetXaxis()->SetTitle("x_{J#gamma}");
    hratio->GetXaxis()->SetTitleSize(0.09);
    hratio->GetXaxis()->SetLabelSize(0.08);
    hratio->Draw("p e");
    TLine * line = new TLine(ana::unfoldXjBins[0],1,ana::unfoldXjBins[ana::nUnfoldXjBins],1);
    line->SetLineStyle(9);
    line->Draw("same");
    c->SaveAs(pdfPath.c_str());

    fout->cd();
    hPythiaReco->Write();
    hRefoldedTruthPt->Write();
    hratio->Write();
    delete hPythiaReco; delete hRefoldedTruthPt;
  }

  // Last page: refolding chi2/NDF vs iteration count.
  cout << "Refolding niter-dependence scan..." << endl;
  TGraph * gChi2 = new TGraph((int)iterationsToScan.size());
  for (unsigned k = 0; k < iterationsToScan.size(); k++) {
    int iter = iterationsToScan[k];
    TH1D * hUnfoldedIter = unfold_utility::unfoldOnce(response, flatMeasured, iter, Form("hUnfoldedIter_%d", iter));
    TH1D * hRefoldedIter = (TH1D*)response->ApplyToTruth(hUnfoldedIter, Form("hRefoldedIter_%d", iter));
    hRefoldedIter->Add(flatFakes);
    double chi2Sum = 0;
    int nCounted = 0;
    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      TH1D * hM = unfold_utility::unflattenXj(flatMeasured, ipt, Form("hMeasuredScan_%d_%d", iter, ipt));
      TH1D * hR = unfold_utility::unflattenXj(hRefoldedIter, ipt, Form("hRefoldedScan_%d_%d", iter, ipt));
      for (int b = 1; b <= nXjBinsForChi2; b++) {
        double vR = hR->GetBinContent(b);
        double vM = hM->GetBinContent(b);
        double eM = hM->GetBinError(b);
        if (eM <= 0) continue;
        chi2Sum += pow(vR-vM,2)/(eM*eM);
        nCounted++;
      }
      delete hM; delete hR;
    }
    double chi2ndf = nCounted > 0 ? chi2Sum/nCounted : 0;
    gChi2->SetPoint(k, iter, chi2ndf);
    cout << "  niter=" << iter << ": refolding chi2/NDF = " << chi2ndf << endl;
    delete hUnfoldedIter;
    delete hRefoldedIter;
  }

  c->Clear();
  c->cd();
  gPad->SetTicks(1,1);
  gPad->SetLeftMargin(.15);
  gChi2->SetMarkerStyle(20);
  gChi2->SetMarkerColor(kBlack);
  gChi2->SetLineColor(kBlack);
  gChi2->SetLineWidth(2);
  gChi2->GetXaxis()->SetTitle("Bayesian unfolding iterations");
  gChi2->GetYaxis()->SetTitle("#chi^{2}/NDF (Refolded+Fakes vs Measured)");
  double ymax = 0;
  for (int k = 0; k < gChi2->GetN(); k++) { double x,y; gChi2->GetPoint(k,x,y); ymax = std::max(ymax,y); }
  gChi2->SetMinimum(0);
  gChi2->SetMaximum(ymax*1.3);
  gChi2->Draw("APL");
  TLine * lnom = new TLine(niterate, 0, niterate, ymax*1.3);
  lnom->SetLineStyle(9);
  lnom->SetLineColor(kRed);
  lnom->Draw("same");
  d.drawAll({"p+p Run24 Data"},{Form("Jet R=%.1f",ana::JetRs[ir]),
      Form("Nominal: %d iterations (dashed line)",niterate)}, .5, .85, 16, 700);
  c->SaveAs(pdfPath.c_str());
  fout->cd();
  gChi2->Write("gRefoldingChi2");

  c->SaveAs(Form("%s]", pdfPath.c_str()));
  fout->Close();
  cout << "Done. Wrote " << pdfPath << " and " << rootPath << endl;
}
