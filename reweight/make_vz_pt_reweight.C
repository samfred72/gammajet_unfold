#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/treeuser.h"
// The original gammajet project builds its OWN, differently-laid-out drawer/ana classes
// into /home/samson72/root/lib/libgammajet.so, sitting on the same library search path
// as this project's libgammajet_unfold.so. Without forcing which one loads first, ROOT's
// implicit symbol autoload can bind drawer/ana calls to the wrong (mismatched-layout)
// library and segfault - explicit load removes the ambiguity.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Builds a Data/MC reweighting for v_z and photon-cluster p_T, to eventually be applied
// as an event weight when using pythia MC (e.g. jet_object/pho_object filling, or a new
// weight branch downstream). Both quantities are compared as shapes normalized to unit
// area - the ratio is what an event weight needs, not the raw (very different-statistics)
// counts.
//
// MC reco is the three pythia photon-threshold samples (Photon5/10/20) stitched by
// truth_cluster_pt using the same threshmap/threshmap_high thresholds treeuser::
// check_keep_MC applies to the photon leg (src/treeuser.cc:60-73) - each truth-pT range
// is covered by exactly one sample, so summing all three without this cut would
// double-count the overlapping generation ranges.
void make_vz_pt_reweight() {
  const double vzcut = ana::vzcut; // 60 cm, src/ana.h:56
  const int nVzBins = 60; // 2 cm bins across the standard analysis vz window
  const double ptCut = ana::cluster_pt_cut; // 10 GeV, src/ana.h:65
  // Capped at 40 GeV rather than ana::ptBins' 100 GeV edge - checked cluster_pt>40 is
  // <0.3% of entries in both Data and MC (Data: 43/815018 above 40 GeV), so 1 GeV bins
  // out to 100 leave a long near-empty tail with wild bin-to-bin ratio noise instead of
  // finer resolution where the statistics actually are.
  const double ptMax = 40;
  const int nPtBins = 30; // 1 GeV bins, finer than the analysis's own 5-bin ana::ptBins

  TH1D * hVzData = new TH1D("hVzData", ";v_{z} [cm];Events (norm.)", nVzBins, -vzcut, vzcut);
  TH1D * hVzMC   = new TH1D("hVzMC",   ";v_{z} [cm];Events (norm.)", nVzBins, -vzcut, vzcut);
  TH1D * hPtData = new TH1D("hPtData", ";Cluster p_{T} [GeV];Clusters (norm.)", nPtBins, ptCut, ptMax);
  TH1D * hPtMC   = new TH1D("hPtMC",   ";Cluster p_{T} [GeV];Clusters (norm.)", nPtBins, ptCut, ptMax);
  hVzData->Sumw2();
  hVzMC->Sumw2();
  hPtData->Sumw2();
  hPtMC->Sumw2();

  // ---- Data ----
  {
    treeuser tu("Data");
    Long64_t nentries = tu.t->GetEntriesFast();
    for (Long64_t e = 0; e < nentries; e++) {
      tu.t->GetEntry(e);
      hVzData->Fill(tu.vz);
      if (fabs(tu.vz) < vzcut && tu.cluster_pt > ptCut) hPtData->Fill(tu.cluster_pt);
    }
    cout << "Data: " << nentries << " entries" << endl;
  }

  // ---- MC reco (pythia Photon5+Photon10+Photon20, stitched by truth pT) ----
  // Cross-section weights for combining the Photon5/10/20 MC samples - same numbers as
  // drawer.h's scalemap[isphoton=1][sample] and insitu/grid_insitu.C's photon_scale, for
  // sim="pythia". Each sample's raw entries need this per-sample scale (not just the
  // truth-pT stitching cut below) before summing, exactly as buildMCXjByPtBin() does via
  // Fill(value, s.second) in insitu/grid_insitu.C:283.
  const map<string,double> photonScale = {{"Photon5",146359.3},{"Photon10",6944.675},{"Photon20",130.4461}};
  const vector<string> mcSamples = {"Photon5", "Photon10", "Photon20"};
  for (const string & sample : mcSamples) {
    treeuser tu(sample, "pythia");
    double loThresh = tu.threshmap[-1][sample];
    double hiThresh = tu.threshmap_high[-1][sample];
    double scale = photonScale.at(sample);
    Long64_t nentries = tu.t->GetEntriesFast();
    Long64_t nkept = 0;
    for (Long64_t e = 0; e < nentries; e++) {
      tu.t->GetEntry(e);
      if (tu.truth_cluster_pt <= loThresh || tu.truth_cluster_pt >= hiThresh) continue;
      nkept++;
      hVzMC->Fill(tu.vz, scale);
      if (fabs(tu.vz) < vzcut && tu.cluster_pt > ptCut) hPtMC->Fill(tu.cluster_pt, scale);
    }
    cout << sample << ": " << nkept << "/" << nentries
         << " entries kept (truth_cluster_pt in (" << loThresh << "," << hiThresh << ")), scale=" << scale << endl;
  }

  // Normalize to unit area before dividing - the ratio needs to compare shapes, not the
  // very different raw statistics between Data and stitched MC.
  if (hVzData->Integral() > 0) hVzData->Scale(1.0 / hVzData->Integral());
  if (hVzMC->Integral()   > 0) hVzMC->Scale(1.0 / hVzMC->Integral());
  if (hPtData->Integral() > 0) hPtData->Scale(1.0 / hPtData->Integral());
  if (hPtMC->Integral()   > 0) hPtMC->Scale(1.0 / hPtMC->Integral());

  TH1D * hVzWeight = (TH1D*)hVzData->Clone("hVzWeight");
  hVzWeight->SetTitle(";v_{z} [cm];Data / MC weight");
  hVzWeight->Divide(hVzMC);

  TH1D * hPtWeight = (TH1D*)hPtData->Clone("hPtWeight");
  hPtWeight->SetTitle(";Cluster p_{T} [GeV];Data / MC weight");
  hPtWeight->Divide(hPtMC);

  // Sanitize NaN/Inf (empty-MC-bin divisions) before writing/drawing - NaN poisons ROOT's
  // axis auto-ranging and silently produces a blank canvas.
  for (TH1D * h : {hVzWeight, hPtWeight}) {
    for (int b = 0; b <= h->GetNbinsX() + 1; b++) {
      double content = h->GetBinContent(b);
      if (std::isnan(content) || std::isinf(content)) {
        h->SetBinContent(b, 1); // no data to reweight against - leave MC unweighted there
        h->SetBinError(b, 0);
      }
    }
  }

  // Statistics run out above ~25 GeV (the 25-40 GeV bins of hPtWeight swing wildly and
  // some are consistent with zero MC entries) - fit only the well-populated 11-25 GeV
  // range and use the fit, not the raw bins, above that. "expo" (exp(p0+p1*x)) is a
  // natural shape for the ratio of two falling spectra and - unlike a polynomial - stays
  // positive under extrapolation past 25 GeV, which a per-event weight must (checked:
  // pol1 crosses zero around pT=37 GeV; expo's chi2/ndf over the fit range is also ~3x
  // better, 35/12 vs 102/12).
  const double ptFitLow = 11, ptFitHigh = 25;
  TF1 * fPtWeight = new TF1("fPtWeight", "expo", ptFitLow, ptFitHigh);
  hPtWeight->Fit(fPtWeight, "RQ");

  const char * outPath = "/home/samson72/sphnx/gammajet_unfold/reweight/vz_pt_reweight.root";
  TFile * fout = new TFile(outPath, "recreate");
  hVzData->Write();
  hVzMC->Write();
  hVzWeight->Write();
  hPtData->Write();
  hPtMC->Write();
  hPtWeight->Write();
  fPtWeight->Write("fPtWeight");
  fout->Close();
  cout << "Saved " << outPath << endl;

  // Diagnostic plot: shape comparison for both variables, plus the pT weight/fit.
  TCanvas * c = new TCanvas("c", "", 2100, 700);
  c->Divide(3, 1);

  c->cd(1);
  TPad * p1 = new TPad("p1", "", 0, 0, 1, 1);
  p1->Draw();
  p1->cd();
  p1->SetLeftMargin(.15);
  p1->SetBottomMargin(.13);
  gPad->SetTicks(1, 1);
  hVzData->SetLineColor(kBlack);
  hVzData->SetLineWidth(2);
  hVzMC->SetLineColor(kRed);
  hVzMC->SetLineWidth(2);
  hVzData->GetYaxis()->SetRangeUser(0, std::max(hVzData->GetMaximum(), hVzMC->GetMaximum()) * 1.4);
  hVzData->Draw("hist");
  hVzMC->Draw("hist same");
  drawer d1("pythia", "nominal");
  d1.drawAll({"p+p Run24 Data (black), Pythia8 #gamma+jet MC (red)"}, {"|v_{z}| reweight input"}, .35, .85, 14, gPad->GetWh() * 0.8);

  c->cd(2);
  TPad * p2 = new TPad("p2", "", 0, 0, 1, 1);
  p2->Draw();
  p2->cd();
  p2->SetLeftMargin(.15);
  p2->SetBottomMargin(.13);
  gPad->SetTicks(1, 1);
  hPtData->SetLineColor(kBlack);
  hPtData->SetLineWidth(2);
  hPtMC->SetLineColor(kRed);
  hPtMC->SetLineWidth(2);
  hPtData->GetYaxis()->SetRangeUser(0, std::max(hPtData->GetMaximum(), hPtMC->GetMaximum()) * 1.4);
  hPtData->Draw("hist");
  hPtMC->Draw("hist same");
  drawer d2("pythia", "nominal");
  d2.drawAll({"p+p Run24 Data (black), Pythia8 #gamma+jet MC (red)"}, {"Cluster p_{T} reweight input"}, .35, .85, 14, gPad->GetWh() * 0.8);

  c->cd(3);
  TPad * p3 = new TPad("p3", "", 0, 0, 1, 1);
  p3->Draw();
  p3->cd();
  p3->SetLeftMargin(.15);
  p3->SetBottomMargin(.13);
  gPad->SetTicks(1, 1);
  hPtWeight->SetLineColor(kBlack);
  hPtWeight->SetMarkerColor(kBlack);
  hPtWeight->SetMarkerStyle(20);
  hPtWeight->GetYaxis()->SetRangeUser(0, hPtWeight->GetBinContent(1) * 1.6);
  hPtWeight->Draw("e");
  fPtWeight->SetLineColor(kBlue);
  fPtWeight->SetRange(ptFitLow, ptMax); // draw the extrapolation past the fit range too
  fPtWeight->Draw("same");
  TLine * lThresh = new TLine(11, 0, 11, hPtWeight->GetBinContent(1) * 1.6);
  lThresh->SetLineStyle(2);
  lThresh->Draw();
  drawer d3("pythia", "nominal");
  d3.drawAll({"p+p Run24 Data / Pythia8 #gamma+jet MC"},
      {"Cluster p_{T} weight", "black = first bin below 11 GeV, blue = expo fit"}, .35, .85, 14, gPad->GetWh() * 0.8);

  const char * pdfPath = "/home/samson72/sphnx/gammajet_unfold/pdfs/reweight_vz_pt.pdf";
  c->SaveAs(pdfPath);
  cout << "Saved " << pdfPath << endl;
}
