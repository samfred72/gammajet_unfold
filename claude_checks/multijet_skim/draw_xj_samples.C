// x_j shape of each MC sample per leading-pT bin (unit area), with the total MC and Data,
// from analysis.cc's per-sample trees. One page per radius and JER variation, all in one file.
// Usage: root -b -l -q 'draw_xj_samples.C("multijet_analysis_pythia.root", "out.pdf")'
// leadCut/recoilCut: the analysis cuts on the leading jet and the recoil |pT2+pT3| (20/14, or 25/14 for
// analysis --tight), plus jets 2 and 3 >= 7 GeV; applied to Data, whose tree keeps looser events for the
// in-situ scans.
void draw_xj_samples(const char * infile = "multijet_analysis_pythia.root", const char * outpdf = "xj_samples.pdf",
                     double leadCut = 20, double recoilCut = 14) {
  gStyle->SetOptStat(0);
  TFile * f = TFile::Open(infile);
  if (!f || f->IsZombie()) { cout << "Could not open " << infile << endl; return; }

  const int nPt = 4;
  const double ptBins[nPt+1] = {20, 25, 30, 35, 50};
  const int nS = 5;
  const char * samples[nS] = {"Jet5", "Jet8", "Jet12", "Jet20", "Jet30"};
  const int colors[nS] = {kViolet+1, kAzure+1, kGreen+2, kOrange+7, kRed+1};

  TCanvas * c = new TCanvas("c", "", 1500, 900);
  c->SaveAs(Form("%s[", outpdf));
  for (int R : {2, 3, 4, 5, 6, 7, 8}) for (const char * sys : {"RECO", "HIGH", "LOW"}) {
    c->Clear();
    c->Divide(3, 2);
    for (int k = 0; k < nPt; k++) {
      TVirtualPad * p = c->cd(k+1);
      p->SetLeftMargin(0.14); p->SetRightMargin(0.03); p->SetTopMargin(0.07); p->SetBottomMargin(0.13);
      TString sel = Form("weight*(leadingPT>%g && leadingPT<%g)", ptBins[k], ptBins[k+1]);

      TH1D * h[nS] = {nullptr};
      double w[nS] = {0}, wtot = 0;
      long n[nS] = {0};
      TH1D * hMC = new TH1D(Form("hmc_%d_%s_%d", R, sys, k), ";x_{j};Normalized counts", 45, 0.4, 2.65);
      hMC->Sumw2();
      for (int s = 0; s < nS; s++) {
        TTree * t = (TTree*)f->Get(Form("ttree_%s_r%d_%s", samples[s], R, sys));
        if (!t || t->GetEntries() == 0) continue;
        h[s] = new TH1D(Form("h_%d_%s_%d_%d", R, sys, k, s), ";x_{j};Normalized counts", 45, 0.4, 2.65);
        h[s]->Sumw2();
        n[s] = t->Draw(Form("leadingPT/PT23>>h_%d_%s_%d_%d", R, sys, k, s), sel, "goff");
        w[s] = h[s]->Integral();
        wtot += w[s];
        hMC->Add(h[s]);
      }
      TH1D * hD = new TH1D(Form("hd_%d_%s_%d", R, sys, k), ";x_{j};Normalized counts", 45, 0.4, 2.65);
      hD->Sumw2();
      TTree * td = (TTree*)f->Get(Form("ttree_data_r%d", R));
      TString dataSel = Form("(leadingPT>%g && leadingPT<%g && leadingPT>=%g && SLPT>=7 && SSLPT>=7 && PT23>=%g)",
                             ptBins[k], ptBins[k+1], leadCut, recoilCut);
      long nD = td->Draw(Form("leadingPT/PT23>>hd_%d_%s_%d", R, sys, k), dataSel, "goff");

      double ymax = 0;
      auto norm = [&](TH1D * x) {
        for (int b = 0; b <= x->GetNbinsX()+1; b++)
          if (!std::isfinite(x->GetBinContent(b))) { x->SetBinContent(b, 0); x->SetBinError(b, 0); }
        if (x->Integral() > 0) x->Scale(1.0/x->Integral());
        ymax = std::max(ymax, x->GetMaximum());
      };
      for (int s = 0; s < nS; s++) if (h[s] && w[s] > 0) norm(h[s]);
      norm(hMC); norm(hD);

      hMC->SetLineColor(kBlack); hMC->SetLineStyle(2); hMC->SetLineWidth(2);
      hMC->SetMinimum(0); hMC->SetMaximum(1.5*ymax);
      hMC->GetXaxis()->SetTitleSize(0.05); hMC->GetYaxis()->SetTitleSize(0.05); hMC->GetYaxis()->SetTitleOffset(1.35);
      hMC->Draw("hist");
      TLegend * l = new TLegend(0.47, 0.50, 0.97, 0.86);
      l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.037);
      for (int s = 0; s < nS; s++) {
        if (!h[s] || w[s] <= 0) continue;
        h[s]->SetLineColor(colors[s]); h[s]->SetLineWidth(2);
        h[s]->Draw("hist same");
        l->AddEntry(h[s], Form("%s: %.1f%% of MC, %ld ev", samples[s], 100*w[s]/wtot, n[s]), "l");
      }
      hD->SetMarkerStyle(20); hD->SetMarkerSize(0.6); hD->SetLineColor(kBlack);
      hD->Draw("e1 same");
      hMC->Draw("hist same");
      l->AddEntry(hMC, "total MC", "l");
      l->AddEntry(hD, Form("Data, %ld ev", nD), "pl");
      l->Draw();
      TLatex tx; tx.SetNDC(); tx.SetTextSize(0.045);
      tx.DrawLatex(0.16, 0.88, Form("R = 0.%d, %.0f < p_{T}^{lead} < %.0f GeV, JER %s", R, ptBins[k], ptBins[k+1], sys));
    }
    c->SaveAs(outpdf);
  }
  c->SaveAs(Form("%s]", outpdf));
}
