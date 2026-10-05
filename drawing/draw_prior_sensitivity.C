#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Prior sensitivity, the data-driven reweighting of PPG08 (Sec. 4.4/6.3) and ATLAS
// (arXiv:1809.07280, Sec. 5.3), without ATLAS's smoothing fit:
//   1. Unfold Data through the nominal response (the analysis result).
//   2. w = (shape-normalized unfolded Data) / (shape-normalized Pythia8 truth), per bin.
//   3. Reweight the response's truth axis by w^0.5, w, w^1.5 and re-project onto reco; training
//      fakes are added back unweighted (no truth to reweight by).
//   4. Re-unfold the same Data through each (the prior is the response's truth marginal).
//   5. Only w vs nominal is the systematic (one symmetric source in draw_systematics.C); our
//      nominal is un-reweighted, so sqrt(w) and w^1.5 are diagnostics only.

// Jet radius index, from jetRadiusIndex (default R=0.4).
int ir = 2;
const int nPtBinsUsed = ana::nPtBinsUsed; // reported bins start at ana::firstUsedPtBin
const int niterate = 2; // nominal iteration count
const vector<int> iterationsToScan = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; // for the bonus niter-dependence page
const vector<double> exponents = {0.5, 1.0, 1.5}; // ATLAS's sqrt(w), w, w^1.5 - see file header
const vector<string> variantNames = {"sqrtw", "w", "w1p5"};
const vector<string> variantLabels = {"#sqrt{w} prior", "w prior", "w^{1.5} prior"};
const vector<int> variantColors = {kAzure+2, kGreen+2, kRed};
const int wIndex = 1; // the reported variant (w)

// The last 3 xJ bins per pT bin are low-count noise: excluded from chi2 only (still drawn).
const int nXjBinsForChi2 = ana::nUnfoldXjBins - 3;

// w per bin (step 2); 1 where either side is empty.
TH1D * buildDataInformedWeights(TH1D * unfoldedNominal, TH1D * priorTruth) {
  TH1D * w = (TH1D*)priorTruth->Clone("hPriorWeight");
  for (int b = 0; b <= w->GetNbinsX()+1; b++) w->SetBinContent(b, 1.0);
  for (int ipt = 0; ipt < ana::nPtBins; ipt++) {
    TH1D * unfPt   = unfold_utility::unflattenXj(unfoldedNominal, ipt, Form("hUnfPt_tmp_%d", ipt));
    TH1D * truthPt = unfold_utility::unflattenXj(priorTruth, ipt, Form("hTruthPt_tmp_%d", ipt));
    double unfInt   = unfPt->Integral();
    double truthInt = truthPt->Integral();
    for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
      int flatbin = ipt*(ana::nUnfoldXjBins+2) + ixj + 1;
      double u = (unfInt   > 0) ? unfPt->GetBinContent(ixj+1)/unfInt     : 0;
      double t = (truthInt > 0) ? truthPt->GetBinContent(ixj+1)/truthInt : 0;
      if (t > 0) w->SetBinContent(flatbin+1, u/t);
    }
    delete unfPt; delete truthPt;
  }
  return w;
}

// One reweighted-prior response at weight^power (step 3).
void buildReweightedResponse(TH2D * nominalMatrix, TH1D * nominalTruth, TH1D * nominalReco,
    TH1D * fakesTraining, TH1D * baseWeight, double power, const char * tag,
    TH2D *& altMatrixOut, TH1D *& altTruthOut, TH1D *& altRecoOut) {
  TH1D * w = (TH1D*)baseWeight->Clone(Form("hWeight_%s", tag));
  for (int b = 0; b <= w->GetNbinsX()+1; b++)
    w->SetBinContent(b, pow(std::max(w->GetBinContent(b), 0.), power));

  altTruthOut = (TH1D*)nominalTruth->Clone(Form("hAltTruth_%s", tag));
  for (int b = 1; b <= altTruthOut->GetNbinsX(); b++)
    altTruthOut->SetBinContent(b, nominalTruth->GetBinContent(b) * w->GetBinContent(b));

  altMatrixOut = (TH2D*)nominalMatrix->Clone(Form("hAltMatrix_%s", tag));
  for (int by = 1; by <= altMatrixOut->GetNbinsY(); by++) {
    double wgt = w->GetBinContent(by);
    for (int bx = 1; bx <= altMatrixOut->GetNbinsX(); bx++)
      altMatrixOut->SetBinContent(bx, by, nominalMatrix->GetBinContent(bx,by) * wgt);
  }

  TH1D * altMatched = (TH1D*)nominalReco->Clone(Form("hAltMatched_%s", tag));
  for (int bx = 1; bx <= altMatrixOut->GetNbinsX(); bx++) {
    double sum = 0;
    for (int by = 1; by <= altMatrixOut->GetNbinsY(); by++) sum += altMatrixOut->GetBinContent(bx,by);
    altMatched->SetBinContent(bx, sum);
  }
  altRecoOut = (TH1D*)altMatched->Clone(Form("hAltReco_%s", tag));
  altRecoOut->Add(fakesTraining); // unweighted: fakes have no truth partner
  delete altMatched;
  delete w;
}

// chi2/NDF of a variant vs nominal, with the nominal's statistical error as yardstick.
double computeChi2NDF(TH1D * hVariant, TH1D * hNominal) {
  double chi2 = 0;
  int ndf = 0;
  for (int b = 1; b <= nXjBinsForChi2; b++) {
    double vV = hVariant->GetBinContent(b);
    double vN = hNominal->GetBinContent(b);
    double eN = hNominal->GetBinError(b);
    if (eN <= 0) continue;
    chi2 += pow(vV-vN,2)/(eN*eN);
    ndf++;
  }
  return (ndf > 0) ? chi2/ndf : 0;
}

void draw_prior_sensitivity(string systag = "nominal", int jetRadiusIndex = 2) {
  gStyle->SetOptStat(0);
  ir = jetRadiusIndex;

  drawer d("pythia", systag);
  // R=0.4 keeps the unsuffixed filenames; other radii get their own.
  string pdfPath  = (ir == 2) ? Form("%s/pdfs/draw_prior_sensitivity_%s.pdf", ana::dir(), systag.c_str())
                              : Form("%s/pdfs/draw_prior_sensitivity_%s_%s.pdf", ana::dir(), systag.c_str(), ana::rnames[ir]);
  string rootPath = (ir == 2) ? Form("%s/hists/prior_sensitivity_%s.root", ana::dir(), systag.c_str())
                              : Form("%s/hists/prior_sensitivity_%s_%s.root", ana::dir(), systag.c_str(), ana::rnames[ir]);

  // Response and purity-corrected Data, as in draw_purity_corrected.C.
  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i",ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i",ir), 1);
  RooUnfoldResponse * response = new RooUnfoldResponse(respRecoTemplate, respTruthTemplate, respMatrix2D);
  TH1D * flatFakesTraining = (TH1D*)response->Hfakes()->Clone("hFakesTrainingFlat");

  TH1D * flatA = d.get(Form("hrecoxj%i_0",ir), 0);
  TH1D * flatC = d.get(Form("hrecoxj%i_2",ir), 0);
  TH1D * flatMeasured = unfold_utility::buildFullyCorrected(flatA, flatC, "data", systag, ir);

  // Step 1: nominal.
  TH1D * flatUnfoldedNominal = unfold_utility::unfoldOnce(response, flatMeasured, niterate, "hUnfoldedNominal");

  // Step 2: base weight.
  TH1D * baseWeight = buildDataInformedWeights(flatUnfoldedNominal, respTruthTemplate);

  // Steps 3-4: reweighted responses, same flatMeasured.
  vector<TH1D*> flatUnfoldedVariant(exponents.size());
  for (unsigned iv = 0; iv < exponents.size(); iv++) {
    TH2D * altMatrix; TH1D * altTruth; TH1D * altReco;
    buildReweightedResponse(respMatrix2D, respTruthTemplate, respRecoTemplate, flatFakesTraining,
        baseWeight, exponents[iv], variantNames[iv].c_str(), altMatrix, altTruth, altReco);
    RooUnfoldResponse * altResponse = new RooUnfoldResponse(altReco, altTruth, altMatrix);
    flatUnfoldedVariant[iv] = unfold_utility::unfoldOnce(altResponse, flatMeasured, niterate,
        Form("hUnfolded_%s", variantNames[iv].c_str()));
  }

  TFile * fout = TFile::Open(rootPath.c_str(), "RECREATE");
  TCanvas * c = new TCanvas("c","",700,900);
  c->SaveAs(Form("%s[", pdfPath.c_str()));

  cout << "Prior sensitivity (niter=" << niterate << "): pT bin, chi2/NDF per variant (first "
       << nXjBinsForChi2 << " of " << ana::nUnfoldXjBins << " xJ bins)" << endl;

  // Pages 1..nPtBinsUsed: nominal vs the three variants (absolute densities) with ratios.
  // Only the w variant (wIndex) enters fracDiff.
  vector<vector<double>> fracDiff(nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0));
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hNominal = unfold_utility::unflattenXj(flatUnfoldedNominal, ipt, Form("hNominal_pt%d", ipt));
    TH1D * hNominalDisp = unfold_utility::densityForDisplay(hNominal, Form("hNominalDisp_pt%d", ipt));

    vector<TH1D*> hVariant(exponents.size()), hVariantDisp(exponents.size());
    double ymax = hNominalDisp->GetMaximum();
    int idisplay = ipt - ana::firstUsedPtBin;
    cout << "  pt" << ipt << " (" << ana::ptBins[ipt] << "-" << ana::ptBins[ipt+1] << " GeV):";
    for (unsigned iv = 0; iv < exponents.size(); iv++) {
      hVariant[iv] = unfold_utility::unflattenXj(flatUnfoldedVariant[iv], ipt, Form("h%s_pt%d", variantNames[iv].c_str(), ipt));
      hVariantDisp[iv] = unfold_utility::densityForDisplay(hVariant[iv], Form("h%s_pt%d_disp", variantNames[iv].c_str(), ipt));
      ymax = std::max(ymax, hVariantDisp[iv]->GetMaximum());
      double chi2ndf = computeChi2NDF(hVariant[iv], hNominal);
      cout << " " << variantNames[iv] << " chi2/NDF=" << chi2ndf;

      if ((int)iv != wIndex) continue; // sqrt(w)/w^1.5: diagnostic only, see file header
      for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
        double vn = hNominal->GetBinContent(ixj+1);
        double vv = hVariant[iv]->GetBinContent(ixj+1);
        if (vn <= 0) continue;
        fracDiff[idisplay][ixj] = (vv-vn)/vn;
      }
    }
    cout << endl;

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
    hNominalDisp->SetLineColor(kBlack);
    hNominalDisp->SetMarkerColor(kBlack);
    hNominalDisp->SetMarkerStyle(20);
    hNominalDisp->SetLineWidth(2);
    hNominalDisp->GetXaxis()->SetLabelSize(0);
    hNominalDisp->GetXaxis()->SetTitle("");
    hNominalDisp->GetYaxis()->SetRangeUser(0, ymax*1.4);
    hNominalDisp->Draw("p e");
    TLegend * l = new TLegend(.45,.58,.85,.80);
    l->SetLineWidth(0);
    l->SetTextSize(0.03);
    l->AddEntry(hNominalDisp, "Nominal (un-reweighted prior)", "lp");
    for (unsigned iv = 0; iv < exponents.size(); iv++) {
      hVariantDisp[iv]->SetLineColor(variantColors[iv]);
      hVariantDisp[iv]->SetLineWidth(2);
      hVariantDisp[iv]->Draw("hist same");
      l->AddEntry(hVariantDisp[iv], variantLabels[iv].c_str(), "l");
    }
    l->Draw();
    d.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
        Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir])}, .18, .85, 14, gPad->GetWh()*0.8);

    p2->cd();
    p2->SetTopMargin(0.02);
    p2->SetBottomMargin(0.3);
    p2->SetLeftMargin(.15);
    gPad->SetTicks(1,1);
    TH1D * hratioFrame = (TH1D*)hNominalDisp->Clone(Form("hratioFrame_pt%d", ipt));
    hratioFrame->Reset("ICES");
    hratioFrame->SetLineColor(kWhite);
    hratioFrame->SetMarkerColor(kWhite);
    hratioFrame->GetYaxis()->SetRangeUser(0.5,1.5);
    hratioFrame->GetYaxis()->SetTitle("Variant / Nominal");
    hratioFrame->GetYaxis()->SetTitleSize(0.09);
    hratioFrame->GetYaxis()->SetTitleOffset(0.7);
    hratioFrame->GetYaxis()->SetLabelSize(0.08);
    hratioFrame->GetXaxis()->SetTitle("x_{J#gamma}");
    hratioFrame->GetXaxis()->SetTitleSize(0.09);
    hratioFrame->GetXaxis()->SetLabelSize(0.08);
    hratioFrame->Draw("p");
    for (unsigned iv = 0; iv < exponents.size(); iv++) {
      TH1D * hratio = (TH1D*)hVariantDisp[iv]->Clone(Form("hratio_%s_pt%d", variantNames[iv].c_str(), ipt));
      hratio->Divide(hNominalDisp);
      hratio->SetLineColor(variantColors[iv]);
      hratio->SetMarkerColor(variantColors[iv]);
      hratio->SetMarkerStyle(20);
      hratio->Draw("p e same");
      fout->cd();
      hratio->Write();
    }
    TLine * line = new TLine(ana::unfoldXjBins[0],1,ana::unfoldXjBins[ana::nUnfoldXjBins],1);
    line->SetLineStyle(9);
    line->Draw("same");
    c->SaveAs(pdfPath.c_str());

    fout->cd();
    hNominal->Write();
    for (unsigned iv = 0; iv < exponents.size(); iv++) hVariant[iv]->Write();
    delete hNominal;
    for (unsigned iv = 0; iv < exponents.size(); iv++) delete hVariant[iv];
  }

  // Summary: the w variant's fractional shift (the systematic), drawn with its negation as a
  // symmetric source.
  c->Clear();
  c->cd();
  vector<TPad*> pads(nPtBinsUsed);
  for (int idisplay = 0; idisplay < nPtBinsUsed; idisplay++) {
    double y1 = 1.0 - (idisplay+1)*(1.0/nPtBinsUsed);
    double y2 = 1.0 - idisplay*(1.0/nPtBinsUsed);
    pads[idisplay] = new TPad(Form("psum_%d",idisplay),"",0,y1,1,y2);
    pads[idisplay]->Draw();
  }
  for (int idisplay = 0; idisplay < nPtBinsUsed; idisplay++) {
    int ipt = idisplay + ana::firstUsedPtBin;
    pads[idisplay]->cd();
    pads[idisplay]->SetLeftMargin(.15);
    pads[idisplay]->SetBottomMargin(idisplay == nPtBinsUsed-1 ? 0.2 : 0.02);
    pads[idisplay]->SetTopMargin(0.05);
    gPad->SetTicks(1,1);

    // Binning only; contents overwritten below.
    TH1D * hNominalPt = unfold_utility::unflattenXj(flatUnfoldedNominal, ipt, Form("hFracFrame_pt%d",ipt));
    TH1D * hFrac    = (TH1D*)hNominalPt->Clone(Form("hFracDiff_pt%d",ipt));
    TH1D * hFracNeg = (TH1D*)hNominalPt->Clone(Form("hFracDiffNeg_pt%d",ipt));
    delete hNominalPt;
    double ymax = 0.05;
    for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
      hFrac->SetBinContent(ixj+1, fracDiff[idisplay][ixj]);
      hFrac->SetBinError(ixj+1, 0);
      hFracNeg->SetBinContent(ixj+1, -fracDiff[idisplay][ixj]);
      hFracNeg->SetBinError(ixj+1, 0);
      ymax = std::max(ymax, std::abs(fracDiff[idisplay][ixj]));
    }
    hFrac->GetYaxis()->SetRangeUser(-ymax*1.3, ymax*1.3);
    hFrac->GetYaxis()->SetTitle("(w - Nom.)/Nom.");
    hFrac->GetYaxis()->SetTitleSize(0.08);
    hFrac->GetYaxis()->SetTitleOffset(0.8);
    hFrac->GetYaxis()->SetLabelSize(0.07);
    hFrac->GetXaxis()->SetTitle(idisplay == nPtBinsUsed-1 ? "x_{J#gamma}" : "");
    hFrac->GetXaxis()->SetLabelSize(idisplay == nPtBinsUsed-1 ? 0.07 : 0);
    hFrac->GetXaxis()->SetTitleSize(0.08);
    hFrac->SetLineColor(kBlack);
    hFrac->SetLineWidth(2);
    hFrac->Draw("hist");
    hFracNeg->SetLineColor(kGray+1);
    hFracNeg->SetLineStyle(2);
    hFracNeg->SetLineWidth(2);
    hFracNeg->Draw("hist same");
    TLine * zero = new TLine(ana::unfoldXjBins[0],0,ana::unfoldXjBins[ana::nUnfoldXjBins],0);
    zero->SetLineStyle(9);
    zero->Draw("same");
    TLatex * t = new TLatex(.18,.85,Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]));
    t->SetNDC();
    t->SetTextFont(43);
    t->SetTextSize(16);
    t->Draw();
    fout->cd();
    hFrac->Write(Form("hPriorSensFracDiff_pt%d",ipt));
  }
  c->SaveAs(pdfPath.c_str());

  // Bonus page: mean chi2/NDF of the w variant vs nominal at each niter (should fall with niter).
  cout << "Prior-sensitivity niter-dependence scan (w-power variant)..." << endl;
  TGraph * gChi2 = new TGraph((int)iterationsToScan.size());
  for (unsigned k = 0; k < iterationsToScan.size(); k++) {
    int iter = iterationsToScan[k];
    TH1D * hNomIter = unfold_utility::unfoldOnce(response, flatMeasured, iter, Form("hNomIter_%d", iter));
    TH2D * altMatrix; TH1D * altTruth; TH1D * altReco;
    buildReweightedResponse(respMatrix2D, respTruthTemplate, respRecoTemplate, flatFakesTraining,
        baseWeight, 1.0, Form("wscan_%d", iter), altMatrix, altTruth, altReco);
    RooUnfoldResponse * altResponse = new RooUnfoldResponse(altReco, altTruth, altMatrix);
    TH1D * hVarIter = unfold_utility::unfoldOnce(altResponse, flatMeasured, iter, Form("hVarIter_%d", iter));

    double chi2Sum = 0;
    int nCounted = 0;
    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      TH1D * hN = unfold_utility::unflattenXj(hNomIter, ipt, Form("hNomScan_%d_%d", iter, ipt));
      TH1D * hV = unfold_utility::unflattenXj(hVarIter, ipt, Form("hVarScan_%d_%d", iter, ipt));
      for (int b = 1; b <= nXjBinsForChi2; b++) {
        double eN = hN->GetBinError(b);
        if (eN <= 0) continue;
        chi2Sum += pow(hV->GetBinContent(b)-hN->GetBinContent(b),2)/(eN*eN);
        nCounted++;
      }
      delete hN; delete hV;
    }
    double chi2ndf = nCounted > 0 ? chi2Sum/nCounted : 0;
    gChi2->SetPoint(k, iter, chi2ndf);
    cout << "  niter=" << iter << ": prior-sensitivity chi2/NDF = " << chi2ndf << endl;
    delete hNomIter; delete hVarIter;
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
  gChi2->GetYaxis()->SetTitle("Prior-sensitivity #chi^{2}/NDF (w vs nominal)");
  double ymaxg = 0;
  for (int k = 0; k < gChi2->GetN(); k++) { double x,y; gChi2->GetPoint(k,x,y); ymaxg = std::max(ymaxg,y); }
  gChi2->SetMinimum(0);
  gChi2->SetMaximum(ymaxg*1.3);
  gChi2->Draw("APL");
  TLine * lnom = new TLine(niterate, 0, niterate, ymaxg*1.3);
  lnom->SetLineStyle(9);
  lnom->SetLineColor(kRed);
  lnom->Draw("same");
  d.drawAll({"p+p Run24 Data"},{Form("Jet R=%.1f",ana::JetRs[ir]),
      Form("Nominal: %d iterations (dashed line)",niterate)}, .5, .85, 16, 700);
  c->SaveAs(pdfPath.c_str());
  fout->cd();
  gChi2->Write("gPriorSensChi2");

  c->SaveAs(Form("%s]", pdfPath.c_str()));
  fout->Close();
  cout << "Done. Wrote " << pdfPath << " and " << rootPath << endl;
}
