#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/systematics_config.h"
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Systematic-uncertainty pages from macros/compute_systematics.C's systematics.root: per source and
// pT bin, the variation vs nominal with their ratio; then every source's fractional difference with
// the combined Purity, Unfolding and total uncertainties, one stacked panel per pT bin (plus the
// 20-25 GeV panel alone, for talks). Drawing only.
// Output: pdfs/draw_systematics[_<rname>].pdf, pdfs/syst_total_<rname>_pt20_25.pdf

void draw_systematics(int ir = 2) {
  gStyle->SetOptStat(0);
  const int nPtBinsUsed = ana::nPtBinsUsed;
  // R=0.4 keeps the unsuffixed filenames; other radii get _<rname>.
  string pdfPathStr  = (ir == 2) ? ana::path("pdfs/draw_systematics.pdf")
                                 : Form("%s/pdfs/draw_systematics_%s.pdf", ana::dir(), ana::rnames[ir]);
  string rootPathStr = (ir == 2) ? ana::path("hists/systematics.root")
                                 : Form("%s/hists/systematics_%s.root", ana::dir(), ana::rnames[ir]);
  const char * pdfPath = pdfPathStr.c_str();
  TFile * fin = TFile::Open(rootPathStr.c_str());
  if (!fin || fin->IsZombie()) {
    cout << "ERROR: cannot open " << rootPathStr << " - run macros/compute_systematics.C(" << ir << ") first." << endl;
    return;
  }
  auto get = [&](const char * name) { return (TH1D*)fin->Get(name); };

  // Generic drawer, used only for labels.
  drawer dLabel;
  TCanvas * c = new TCanvas("c","",700,900);
  c->SaveAs(Form("%s[", pdfPath));

  vector<TH1D*> nominalDisp(ana::nPtBins);
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++)
    nominalDisp[ipt] = get(Form("hNominalDisp_pt%d", ipt));

  // Per source and pT bin: variation vs nominal, and their ratio.
  map<string, vector<TH1D*>> fracDiff;
  for (const string & systag : systematics) {
    fracDiff[systag] = vector<TH1D*>(ana::nPtBins);
    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      TH1D * hVarDisp = get(Form("h%s_pt%d_disp", systag.c_str(), ipt));
      TH1D * hratio   = get(Form("hratio_%s_pt%d", systag.c_str(), ipt));
      fracDiff[systag][ipt] = get(Form("hfracdiff_%s_pt%d", systag.c_str(), ipt));
      if (!hVarDisp || !hratio || !fracDiff[systag][ipt]) continue; // skipped in compute_systematics.C

      c->Clear();
      c->cd();
      TPad * p1 = new TPad(Form("p1_%s_%d",systag.c_str(),ipt),"",0,.35,1,1);
      TPad * p2 = new TPad(Form("p2_%s_%d",systag.c_str(),ipt),"",0,0,1,.35);
      p1->Draw();
      p2->Draw();

      p1->cd();
      p1->SetBottomMargin(0.02);
      p1->SetLeftMargin(.15);
      gPad->SetTicks(1,1);
      nominalDisp[ipt]->SetLineColor(kBlack);
      nominalDisp[ipt]->SetMarkerColor(kBlack);
      nominalDisp[ipt]->SetMarkerStyle(20);
      nominalDisp[ipt]->SetLineWidth(2);
      nominalDisp[ipt]->GetXaxis()->SetLabelSize(0);
      nominalDisp[ipt]->GetXaxis()->SetTitle("");
      nominalDisp[ipt]->GetYaxis()->SetRangeUser(0, 0.25);
      nominalDisp[ipt]->Draw("p e");
      hVarDisp->SetLineColor(systColors.at(systag));
      hVarDisp->SetMarkerColor(systColors.at(systag));
      hVarDisp->SetMarkerStyle(21);
      hVarDisp->SetLineWidth(2);
      hVarDisp->Draw("p e same");
      TLegend * l = new TLegend(.55,.55,.85,.70);
      l->SetLineWidth(0);
      l->SetTextSize(0.035);
      l->AddEntry(nominalDisp[ipt], "Nominal");
      l->AddEntry(hVarDisp, systag.c_str());
      l->Draw();
      dLabel.drawAll({ana::dataLabel()},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
          },
          .5, .85, 16, gPad->GetWh()*0.8);

      p2->cd();
      p2->SetTopMargin(0.02);
      p2->SetBottomMargin(0.3);
      p2->SetLeftMargin(.15);
      gPad->SetTicks(1,1);
      hratio->SetLineColor(kBlack);
      hratio->SetMarkerColor(kBlack);
      hratio->SetMarkerStyle(20);
      hratio->GetYaxis()->SetRangeUser(0.5,1.5);
      hratio->GetYaxis()->SetTitle("Variation / Nominal");
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
      c->SaveAs(pdfPath);
    }
  }

  vector<TH1D*> purityUnc(ana::nPtBins), unfoldingUncUp(ana::nPtBins), unfoldingUncDown(ana::nPtBins);
  vector<TH1D*> totalUncUp(ana::nPtBins), totalUncDown(ana::nPtBins);
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    purityUnc[ipt]        = get(Form("hpurity_pt%d", ipt));
    unfoldingUncUp[ipt]   = get(Form("hunfolding_up_pt%d", ipt));
    unfoldingUncDown[ipt] = get(Form("hunfolding_down_pt%d", ipt));
    totalUncUp[ipt]       = get(Form("hquadsum_up_pt%d", ipt));
    totalUncDown[ipt]     = get(Form("hquadsum_down_pt%d", ipt));
  }

  // Final page: every source's fractional difference plus the totals, one panel per pT bin.
  // No error bars: these are ratios of already-unfolded results.
  c->Clear();
  c->cd();
  vector<TPad*> pads(nPtBinsUsed);
  for (int ipt = 0; ipt < nPtBinsUsed; ipt++) {
    double y1 = 1.0 - (ipt+1)*(1.0/nPtBinsUsed);
    double y2 = 1.0 - ipt*(1.0/nPtBinsUsed);
    pads[ipt] = new TPad(Form("psum_%d",ipt),"",0,y1,1,y2);
    pads[ipt]->Draw();
  }
  // One summary panel for pT bin ipt in the current pad; standalone = its own page.
  auto drawSummaryPanel = [&](int ipt, bool standalone) {
    // idisplay: position among the stacked pads; ipt: the ana::ptBins index.
    int idisplay = ipt - ana::firstUsedPtBin;
    bool xaxis = standalone || idisplay == nPtBinsUsed-1;
    // unique clone names for the standalone page
    const char * sfx = standalone ? "_single" : "";
    gPad->SetTicks(1,1);

    double ymax = 0.05; // headroom floor
    for (const string & systag : systematics)
      ymax = std::max(ymax, fracDiff[systag][ipt]->GetMaximum());
    for (const string & systag : systematics)
      ymax = std::max(ymax, -fracDiff[systag][ipt]->GetMinimum());
    ymax = std::max(ymax, purityUnc[ipt]->GetMaximum());
    ymax = std::max(ymax, unfoldingUncUp[ipt]->GetMaximum());
    ymax = std::max(ymax, unfoldingUncDown[ipt]->GetMaximum());
    ymax = std::max(ymax, totalUncUp[ipt]->GetMaximum());
    ymax = std::max(ymax, totalUncDown[ipt]->GetMaximum());

    TH1D * frame = (TH1D*)fracDiff[systematics[0]][ipt]->Clone(Form("hsumframe_pt%d%s",ipt,sfx));
    frame->Reset("ICES");
    frame->SetLineColor(kWhite);
    frame->GetYaxis()->SetRangeUser(-1.3, 1.3);
    frame->GetYaxis()->SetTitle("(Var.-Nom.)/Nom.");
    // stacked panels are ~1/3 of the canvas, so larger text fractions
    frame->GetYaxis()->SetTitleSize(standalone ? 0.05 : 0.08);
    frame->GetYaxis()->SetTitleOffset(standalone ? 1.3 : 0.8);
    frame->GetYaxis()->SetLabelSize(standalone ? 0.045 : 0.07);
    frame->GetXaxis()->SetTitle(xaxis ? "x_{J#gamma}" : "");
    frame->GetXaxis()->SetLabelSize(xaxis ? (standalone ? 0.045 : 0.07) : 0);
    frame->GetXaxis()->SetTitleSize(standalone ? 0.05 : 0.08);
    frame->Draw("p");

    TLine * zero = new TLine(ana::unfoldXjBins[0],0,ana::unfoldXjBins[ana::nUnfoldXjBins],0);
    zero->SetLineStyle(9);
    zero->Draw("same");

    // Stacked page: sPHENIX label in the first panel, legend only in the second.
    TLegend * ls = standalone ? new TLegend(.17,.72,.62,.93) : new TLegend(.17,.60,.62,.93);
    ls->SetLineWidth(0);
    ls->SetFillStyle(0);
    ls->SetTextSize(standalone ? 0.035 : 0.05);
    ls->SetNColumns(2);
    for (const DisplayGroup & g : displayGroups) {
      bool firstMember = true;
      for (const string & systag : g.members) {
        // Display clone without errors; fracDiff keeps its errors in the file.
        TH1D * hf = (TH1D*)fracDiff[systag][ipt]->Clone(Form("hfracdiff_%s_pt%d_disp%s", systag.c_str(), ipt, sfx));
        for (int b = 1; b <= hf->GetNbinsX(); b++) hf->SetBinError(b, 0);
        hf->SetLineColor(g.color);
        hf->SetLineWidth(1);
        hf->Draw("hist same");
        if (firstMember) { ls->AddEntry(hf, g.label.c_str(), "l"); firstMember = false; }

        if (g.symmetric) {
          TH1D * hfNeg = (TH1D*)hf->Clone(Form("hfracdiff_%s_pt%d_disp_neg%s", systag.c_str(), ipt, sfx));
          hfNeg->Scale(-1);
          hfNeg->Draw("hist same");
        }
      }
    }
    TH1D * hPurityDisp = (TH1D*)purityUnc[ipt]->Clone(Form("hpurity_pt%d_disp%s", ipt, sfx));
    for (int b = 1; b <= hPurityDisp->GetNbinsX(); b++) hPurityDisp->SetBinError(b, 0);
    hPurityDisp->SetLineColor(colorPurity);
    hPurityDisp->SetLineWidth(1);
    hPurityDisp->Draw("hist same");
    ls->AddEntry(hPurityDisp, "Purity", "l");

    TH1D * hPurityDispNeg = (TH1D*)hPurityDisp->Clone(Form("hpurity_pt%d_disp_neg%s", ipt, sfx));
    hPurityDispNeg->Scale(-1);
    hPurityDispNeg->Draw("hist same");

    // Up and down are independent; down is negated for display.
    TH1D * hUnfoldingUpDisp = (TH1D*)unfoldingUncUp[ipt]->Clone(Form("hunfolding_up_pt%d_disp%s", ipt, sfx));
    for (int b = 1; b <= hUnfoldingUpDisp->GetNbinsX(); b++) hUnfoldingUpDisp->SetBinError(b, 0);
    hUnfoldingUpDisp->SetLineColor(colorUnfolding);
    hUnfoldingUpDisp->SetLineWidth(1);
    hUnfoldingUpDisp->Draw("hist same");
    ls->AddEntry(hUnfoldingUpDisp, "Unfolding", "l");

    TH1D * hUnfoldingDownDisp = (TH1D*)unfoldingUncDown[ipt]->Clone(Form("hunfolding_down_pt%d_disp%s", ipt, sfx));
    for (int b = 1; b <= hUnfoldingDownDisp->GetNbinsX(); b++) hUnfoldingDownDisp->SetBinError(b, 0);
    hUnfoldingDownDisp->Scale(-1);
    hUnfoldingDownDisp->SetLineColor(colorUnfolding);
    hUnfoldingDownDisp->SetLineWidth(1);
    hUnfoldingDownDisp->Draw("hist same");

    // Up and down are independent; down is negated for display.
    TH1D * hTotalUpDisp = (TH1D*)totalUncUp[ipt]->Clone(Form("hquadsum_up_pt%d_disp%s", ipt, sfx));
    for (int b = 1; b <= hTotalUpDisp->GetNbinsX(); b++) hTotalUpDisp->SetBinError(b, 0);
    hTotalUpDisp->SetLineColor(kBlack);
    hTotalUpDisp->SetLineWidth(2);
    hTotalUpDisp->SetLineStyle(2);

    TH1D * hTotalDownDisp = (TH1D*)totalUncDown[ipt]->Clone(Form("hquadsum_down_pt%d_disp%s", ipt, sfx));
    for (int b = 1; b <= hTotalDownDisp->GetNbinsX(); b++) hTotalDownDisp->SetBinError(b, 0);
    hTotalDownDisp->Scale(-1);
    hTotalDownDisp->SetLineColor(kBlack);
    hTotalDownDisp->SetLineWidth(2);
    hTotalDownDisp->SetLineStyle(2);

    hTotalDownDisp->Draw("hist same");
    hTotalUpDisp->Draw("hist same");
    ls->AddEntry(hTotalUpDisp, "Total", "l");
    if (standalone || idisplay == 1) ls->Draw();

    if (standalone) {
      dLabel.drawAll({ana::dataLabel()},
                     {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
                      Form("Jet R=%.1f",ana::JetRs[ir])}, .19, .33, 18, 600);
      return;
    }
    if (idisplay == 0)
      dLabel.drawAll({ana::dataLabel()}, {Form("Jet R=%.1f",ana::JetRs[ir])}, .18, .86, 16, 300);
    // pT range at the bottom left of every panel
    TLatex * t = new TLatex(.18,gPad->GetBottomMargin()+.06,Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]));
    t->SetNDC();
    t->SetTextFont(43);
    t->SetTextSize(16);
    t->Draw();
  };
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    int idisplay = ipt - ana::firstUsedPtBin;
    pads[idisplay]->cd();
    pads[idisplay]->SetLeftMargin(.15);
    pads[idisplay]->SetBottomMargin(idisplay == nPtBinsUsed-1 ? 0.2 : 0.02);
    pads[idisplay]->SetTopMargin(0.05);
    drawSummaryPanel(ipt, false);
  }
  c->SaveAs(pdfPath);

  // Talk version: the 20-25 GeV panel alone.
  {
    const int iptSingle = ana::findPtBin(22.5);
    TCanvas * cs = new TCanvas("csingle","",800,600);
    cs->SetLeftMargin(.15);
    cs->SetRightMargin(.05);
    cs->SetTopMargin(.05);
    cs->SetBottomMargin(.12);
    drawSummaryPanel(iptSingle, true);
    cs->SaveAs(Form("%s/pdfs/syst_total_%s_pt%.0f_%.0f.pdf", ana::dir(), ana::rnames[ir],
                    ana::ptBins[iptSingle], ana::ptBins[iptSingle+1]));
    c->cd();
  }

  c->SaveAs(Form("%s]", pdfPath));
  cout << "Wrote " << pdfPath << endl;
}
