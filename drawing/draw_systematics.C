#include "../src/ana.h"
#include "../src/drawer.h"
#include "../src/unfold_utility.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Systematic uncertainties: Data's purity-corrected xJ is unfolded through each variation's own
// response and purity and compared with nominal (both shape-normalized, so a change in the
// accepted event count is not read as a shape change). herwig keeps Data at nominal and swaps
// the response generator; niterLow/High change only the iteration count; priorSensitivity is
// read from draw_prior_sensitivity.C's output (run that first, same radius).
//
// Total (per bin):
//   - JER/JES/emscale/EMR high/low: per-source, per-bin sign split (asymmetricSystematics).
//   - Purity: the five ABCD boundary sources in one symmetric quadrature sum (PPG12 Sec. 5.3).
//   - Unfolding: up = niterHigh (+) prior, down = niterLow (+) prior.
//   - threejet, herwig: symmetric.
// Writes the ratio and fractional difference per (source, pT bin) to systematics.root.

// Jet radius index, set by draw_systematics(int); default R=0.4.
int ir = 2;
const int nPtBinsUsed = ana::nPtBinsUsed; // reported bins start at ana::firstUsedPtBin
const int niterate = 2;
// ana::systags minus nominal, plus the sources that are not reprocessings. A new systag needs
// systColors/systSources entries (.at() throws if missing).
const vector<string> systematics = [] {
  vector<string> v;
  for (const string & s : ana::systags) if (s != "nominal") v.push_back(s);
  v.insert(v.end(), {"herwig", "niterLow", "niterHigh", "priorSensitivity"});
  return v;
}();
const map<string,int> systColors = {
  {"JERhigh",      kRed+1},
  {"JERlow",       kOrange+7},
  {"emscale_high", kAzure+1},
  {"emscale_low",  kTeal+2},
  {"EMRhigh",      kGray+2},
  {"EMRlow",       kGray+1},
  {"jes_high",     kPink+6},
  {"jes_low",      kYellow+2},
  {"threejet",     kGreen+2},
  {"narrowBDT",    kBlue+1},
  {"narrowISO",    kMagenta+1},
  {"narrowBDTbkg", kBlue+3},
  {"narrowISObkg", kMagenta+3},
  {"wideISObkg",   kMagenta-4},
  {"herwig",       kSpring+5},
  {"niterLow",     kCyan+2},
  {"niterHigh",    kViolet+1},
  {"priorSensitivity", kOrange+2},
};

// Per source: response generator (sim), Data reprocessing (dataSystag) and niter.
// priorSensitivity has no entry (special-cased below).
struct SystSource { string sim; string dataSystag; int niter; };
const map<string, SystSource> systSources = {
  {"JERhigh",      {"pythia", "JERhigh",      niterate}},
  {"JERlow",       {"pythia", "JERlow",       niterate}},
  {"emscale_high", {"pythia", "emscale_high", niterate}},
  {"emscale_low",  {"pythia", "emscale_low",  niterate}},
  {"EMRhigh",      {"pythia", "EMRhigh",      niterate}},
  {"EMRlow",       {"pythia", "EMRlow",       niterate}},
  {"jes_high",     {"pythia", "jes_high",     niterate}},
  {"jes_low",      {"pythia", "jes_low",      niterate}},
  {"threejet",     {"pythia", "threejet",     niterate}},
  {"narrowBDT",    {"pythia", "narrowBDT",    niterate}},
  {"narrowISO",    {"pythia", "narrowISO",    niterate}},
  {"narrowBDTbkg", {"pythia", "narrowBDTbkg", niterate}},
  {"narrowISObkg", {"pythia", "narrowISObkg", niterate}},
  {"wideISObkg",   {"pythia", "wideISObkg",   niterate}},
  {"herwig",       {"herwig", "nominal",      niterate}},
  {"niterLow",     {"pythia", "nominal",      niterate-1}},
  {"niterHigh",    {"pythia", "nominal",      niterate+1}},
};

// High/low sources split per bin by sign (from ana::asymmetricSystagPairs).
const set<string> asymmetricSystematics = [] {
  set<string> s;
  for (const auto & pr : ana::asymmetricSystagPairs) { s.insert(pr.first); s.insert(pr.second); }
  return s;
}();

// Colorblind-safe 8-hue palette, fixed order.
const int colorBlue      = TColor::GetColor("#2a78d6");
const int colorOrange    = TColor::GetColor("#eb6834");
const int colorAqua      = TColor::GetColor("#1baf7a");
const int colorUnfolding = TColor::GetColor("#eda100");
const int colorMagenta   = TColor::GetColor("#e87ba4");
const int colorPurity    = TColor::GetColor("#008300");
const int colorRed       = TColor::GetColor("#e34948");
const int colorGrey      = TColor::GetColor("#767676");

// Final-page grouping: a high/low pair shares one color and legend entry; a symmetric source is
// drawn with its negation. Purity and Unfolding members are drawn only as the combined curves.
struct DisplayGroup { string label; int color; vector<string> members; bool symmetric; };
const vector<DisplayGroup> displayGroups = {
  {"JER",       colorBlue,    {"JERhigh", "JERlow"},          false},
  {"JES",       colorOrange,  {"jes_high", "jes_low"},        false},
  {"EM scale",  colorAqua,    {"emscale_high", "emscale_low"},false},
  {"EMR",       colorGrey,    {"EMRhigh", "EMRlow"},          false},
  {"threejet",  colorMagenta, {"threejet"},                   true},
  {"herwig",    colorRed,     {"herwig"},                     true},
};
// Combined into one symmetric Purity source.
const set<string> purityMembers = {
  "narrowBDT", "narrowISO", "narrowBDTbkg", "narrowISObkg", "wideISObkg"
};
// Combined into the asymmetric Unfolding source (see unfoldingUncUp/Down).
const set<string> unfoldingMembers = {
  "niterLow", "niterHigh", "priorSensitivity"
};
// R=0.4 keeps the unsuffixed filenames; other radii get _<rname>.
string pdfPathStr, rootPathStr;
const char * pdfPath;
const char * rootPath;

// Data's purity-corrected xJ unfolded through systag's response (sim's generator) with niter
// iterations. Purity always comes from dataSystag.
TH1D * getUnfoldedData(string sim, string systag, int niter, const char * name) {
  drawer d(sim, systag);
  TH1D * respRecoTemplate  = d.get(Form("hrecoxj%i",ir), 1);
  TH1D * respTruthTemplate = d.get(Form("htruthxj%i",ir), 1);
  TH2D * respMatrix2D      = d.get2d(Form("hxjresponse%i",ir), 1);
  TH1D * flatA = d.get(Form("hrecoxj%i_0",ir), 0);
  TH1D * flatC = d.get(Form("hrecoxj%i_2",ir), 0);
  TH1D * flatCorrected = unfold_utility::buildFullyCorrected(flatA, flatC, systag.c_str(), systag, ir);
  return unfold_utility::unfoldOnce(respRecoTemplate, respTruthTemplate, respMatrix2D, flatCorrected, niter, name);
}

void draw_systematics(int jetRadiusIndex = 2) {
  gStyle->SetOptStat(0);

  ir = jetRadiusIndex;
  pdfPathStr  = (ir == 2) ? ana::path("pdfs/draw_systematics.pdf")
                          : Form("%s/pdfs/draw_systematics_%s.pdf", ana::dir(), ana::rnames[ir]);
  rootPathStr = (ir == 2) ? ana::path("hists/systematics.root")
                          : Form("%s/hists/systematics_%s.root", ana::dir(), ana::rnames[ir]);
  pdfPath  = pdfPathStr.c_str();
  rootPath = rootPathStr.c_str();

  // Generic drawer, used only for labels.
  drawer dLabel;

  TFile * fout = TFile::Open(rootPath, "RECREATE");
  TCanvas * c = new TCanvas("c","",700,900);
  c->SaveAs(Form("%s[", pdfPath));

  TH1D * flatNominal = getUnfoldedData("pythia", "nominal", niterate, "hUnfoldedNominal");
  vector<TH1D*> nominalDisp(ana::nPtBins);
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hNom = unfold_utility::unflattenXj(flatNominal, ipt, Form("hNominal_pt%d", ipt));
    nominalDisp[ipt] = unfold_utility::densityForDisplay(hNom, Form("hNominalDisp_pt%d", ipt));
    nominalDisp[ipt]->Scale(1./nominalDisp[ipt]->Integral());
    nominalDisp[ipt]->GetYaxis()->SetTitle("Shape-normalized counts / bin width");
  }

  // fracDiff[systag][ipt]: (variation-nominal)/nominal vs xJ.
  map<string, vector<TH1D*>> fracDiff;

  // Radius-specific (ir==2 keeps the unsuffixed name).
  string priorSensPathStr = (ir == 2) ? ana::path("hists/prior_sensitivity_nominal.root")
                                      : Form("%s/hists/prior_sensitivity_nominal_%s.root", ana::dir(), ana::rnames[ir]);
  TFile * fPriorSens = TFile::Open(priorSensPathStr.c_str());
  if (!fPriorSens || fPriorSens->IsZombie())
    cout << "WARNING: couldn't open " << priorSensPathStr << " - run "
         << "drawing/draw_prior_sensitivity.C(\"nominal\", " << ir << ") first, or the "
         << "priorSensitivity source below will be empty." << endl;

  for (const string & systag : systematics) {
    bool isPriorSens = (systag == "priorSensitivity");
    const SystSource * src = isPriorSens ? nullptr : &systSources.at(systag);
    TH1D * flatVar = isPriorSens ? nullptr : getUnfoldedData(src->sim, src->dataSystag, src->niter, Form("hUnfolded_%s", systag.c_str()));
    fracDiff[systag] = vector<TH1D*>(ana::nPtBins);

    for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
      TH1D * hVar;
      if (isPriorSens) {
        // hw_pt<N>: the raw unfolded "w"-prior result, run through the same pipeline as every other
        // source. Not hPriorSensFracDiff_pt<N>, which is already a ratio-1.
        TH1D * hRaw = fPriorSens ? (TH1D*)fPriorSens->Get(Form("hw_pt%d", ipt)) : nullptr;
        if (!hRaw) {
          cout << "WARNING: hw_pt" << ipt << " missing from prior_sensitivity_nominal.root - skipping pt" << ipt << "." << endl;
          continue;
        }
        hVar = (TH1D*)hRaw->Clone(Form("h%s_pt%d", systag.c_str(), ipt));
      } else {
        hVar = unfold_utility::unflattenXj(flatVar, ipt, Form("h%s_pt%d", systag.c_str(), ipt));
      }
      TH1D * hVarDisp = unfold_utility::densityForDisplay(hVar, Form("h%s_pt%d_disp", systag.c_str(), ipt));
      hVarDisp->Scale(1./hVarDisp->Integral());
      hVarDisp->GetYaxis()->SetTitle("Shape-normalized counts / bin width");

      TH1D * hratio = (TH1D*)hVarDisp->Clone(Form("hratio_%s_pt%d", systag.c_str(), ipt));
      hratio->Divide(nominalDisp[ipt]);
      // Divide leaves 0 where the denominator is 0; NaN masking here blanked the page before.

      TH1D * hfrac = (TH1D*)hratio->Clone(Form("hfracdiff_%s_pt%d", systag.c_str(), ipt));
      for (int b = 1; b <= hfrac->GetNbinsX(); b++) {
        if (hratio->GetBinContent(b) != 0) hfrac->SetBinContent(b, hratio->GetBinContent(b) - 1);
      }
      hfrac->GetYaxis()->SetTitle("(Variation - Nominal) / Nominal");
      fracDiff[systag][ipt] = hfrac;

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
      dLabel.drawAll({"p+p Run24 Data"},{Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
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

      fout->cd();
      hratio->Write();
      hfrac->Write();
    }
  }

  // Purity: symmetric quadrature sum of the five members.
  vector<TH1D*> purityUnc(ana::nPtBins);
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hPurity = (TH1D*)fracDiff[systematics[0]][ipt]->Clone(Form("hpurity_pt%d", ipt));
    hPurity->Reset("ICES");
    for (int b = 1; b <= hPurity->GetNbinsX(); b++) {
      double sumsq = 0;
      for (const string & systag : purityMembers) {
        double v = fracDiff[systag][ipt]->GetBinContent(b);
        sumsq += v*v;
      }
      hPurity->SetBinContent(b, sqrt(sumsq));
    }
    hPurity->GetYaxis()->SetTitle("Purity systematic uncertainty");
    purityUnc[ipt] = hPurity;
    fout->cd();
    hPurity->Write();
  }

  // Unfolding: up = niterHigh (+) prior, down = niterLow (+) prior. niter's direction is coherent
  // across xJ, so no per-bin sign test; the prior has no direction and enters both.
  vector<TH1D*> unfoldingUncUp(ana::nPtBins), unfoldingUncDown(ana::nPtBins);
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    // Null only if draw_prior_sensitivity.C has not been run for this radius.
    if (!fracDiff["priorSensitivity"][ipt]) {
      cout << "ERROR: fracDiff[\"priorSensitivity\"][" << ipt << "] is missing - run "
           << "drawing/draw_prior_sensitivity.C(\"nominal\", " << ir << ") before draw_systematics.C(" << ir << ")." << endl;
      exit(1);
    }
    TH1D * hUp   = (TH1D*)fracDiff[systematics[0]][ipt]->Clone(Form("hunfolding_up_pt%d", ipt));
    TH1D * hDown = (TH1D*)fracDiff[systematics[0]][ipt]->Clone(Form("hunfolding_down_pt%d", ipt));
    hUp->Reset("ICES");
    hDown->Reset("ICES");
    for (int b = 1; b <= hUp->GetNbinsX(); b++) {
      double vHigh  = fracDiff["niterHigh"][ipt]->GetBinContent(b);
      double vLow   = fracDiff["niterLow"][ipt]->GetBinContent(b);
      double vPrior = fracDiff["priorSensitivity"][ipt]->GetBinContent(b);
      hUp->SetBinContent(b, sqrt(vHigh*vHigh + vPrior*vPrior));
      hDown->SetBinContent(b, sqrt(vLow*vLow + vPrior*vPrior));
    }
    hUp->GetYaxis()->SetTitle("Unfolding systematic uncertainty (up)");
    hDown->GetYaxis()->SetTitle("Unfolding systematic uncertainty (down)");
    unfoldingUncUp[ipt] = hUp;
    unfoldingUncDown[ipt] = hDown;
    fout->cd();
    hUp->Write();
    hDown->Write();
  }

  // Total: quadrature sum. Symmetric sources (including combined Purity) feed both totals;
  // asymmetricSystematics sources feed the total matching their sign in that bin.
  vector<TH1D*> totalUncUp(ana::nPtBins), totalUncDown(ana::nPtBins);
  for (int ipt = ana::firstUsedPtBin; ipt < ana::firstUsedPtBin+nPtBinsUsed; ipt++) {
    TH1D * hTotalUp   = (TH1D*)fracDiff[systematics[0]][ipt]->Clone(Form("hquadsum_up_pt%d", ipt));
    TH1D * hTotalDown = (TH1D*)fracDiff[systematics[0]][ipt]->Clone(Form("hquadsum_down_pt%d", ipt));
    hTotalUp->Reset("ICES");
    hTotalDown->Reset("ICES");
    for (int b = 1; b <= hTotalUp->GetNbinsX(); b++) {
      double sumsqUp = 0, sumsqDown = 0;
      for (const string & systag : systematics) {
        if (purityMembers.count(systag)) continue; // folded into purityUnc instead
        if (unfoldingMembers.count(systag)) continue; // folded into unfoldingUncUp/Down instead
        double v = fracDiff[systag][ipt]->GetBinContent(b);
        if (asymmetricSystematics.count(systag)) {
          if (v > 0) sumsqUp   += v*v;
          else       sumsqDown += v*v;
        } else {
          sumsqUp   += v*v;
          sumsqDown += v*v;
        }
      }
      double p = purityUnc[ipt]->GetBinContent(b);
      sumsqUp   += p*p;
      sumsqDown += p*p;
      double uUp   = unfoldingUncUp[ipt]->GetBinContent(b);
      double uDown = unfoldingUncDown[ipt]->GetBinContent(b);
      sumsqUp   += uUp*uUp;
      sumsqDown += uDown*uDown;
      hTotalUp->SetBinContent(b, sqrt(sumsqUp));
      hTotalDown->SetBinContent(b, sqrt(sumsqDown));
    }
    hTotalUp->GetYaxis()->SetTitle("Total systematic uncertainty (up)");
    hTotalDown->GetYaxis()->SetTitle("Total systematic uncertainty (down)");
    totalUncUp[ipt] = hTotalUp;
    totalUncDown[ipt] = hTotalDown;
    fout->cd();
    hTotalUp->Write();
    hTotalDown->Write();
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

    TLegend * ls = standalone ? new TLegend(.17,.72,.62,.93) : new TLegend(.4,.65,.68,.93);
    ls->SetLineWidth(0);
    ls->SetFillStyle(0);
    ls->SetTextSize(standalone ? 0.035 : 0.05);
    if (standalone) ls->SetNColumns(2);
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
    ls->AddEntry(hTotalUpDisp, "Total (asym. quad. sum)", "l");
    ls->Draw();

    if (standalone) {
      dLabel.drawAll({"p+p Run24 Data"},
                     {Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]),
                      Form("Jet R=%.1f",ana::JetRs[ir])}, .19, .33, 18, 600);
      return;
    }
    TLatex * t = new TLatex(.18,.85,Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV",ana::ptBins[ipt],ana::ptBins[ipt+1]));
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
  fout->Close();
  cout << "Done. Wrote " << pdfPath << " and " << rootPath << endl;
}
