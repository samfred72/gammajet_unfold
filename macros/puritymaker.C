#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
// The original gammajet project builds its OWN, differently-laid-out drawer/ana classes
// into /home/samson72/root/lib/libgammajet.so, sitting on the same library search path
// as this project's libgammajet_unfold.so. Without forcing which one loads first, ROOT's
// implicit symbol autoload can bind drawer/ana calls to the wrong (mismatched-layout)
// library and segfault - explicit load removes the ambiguity. See unfold.C for the same
// pattern; it never had this problem because it already did this.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Writes every object via bare ->Write() (implicit current TDirectory) - the caller is
// responsible for cd()'ing into the right target (a subdirectory of the shared
// per-systag purity file, one per jet radius - see puritymaker()) before calling this,
// so it no longer opens its own output file itself.
TGraphAsymmErrors * combine_hists(TH1D * h[], TH1D * f[], TGraphAsymmErrors ** graphCOut = nullptr) {
  TRandom3 * rand = new TRandom3();
  TH1D * hA = h[0];
  TH1D * hB = h[1];
  TH1D * hC = h[2];
  TH1D * hD = h[3];
  TH1D * ha = f[0];
  TH1D * hb = f[1];
  TH1D * hc = f[2];
  TH1D * hd = f[3];
  TH1D * H[ana::nPtBins];
  TH1D * HC[ana::nPtBins]; // bootstrap distributions of region-C purity (S_C/C = c*S/C)
  TH2D * H2 = new TH2D("bootstrap2D",";bin number;bootstrapped value",ana::nPtBins,ana::ptBins,100,-0.2,1.5);
  TGraphAsymmErrors * oh = new TGraphAsymmErrors(hA->GetNbinsX());
  oh->SetName("combined");
  TGraphAsymmErrors * ohC = new TGraphAsymmErrors(hA->GetNbinsX());
  ohC->SetName("combined_C");
  TH1D * oH = (TH1D*)hA->Clone("hcombined");
  oH->Reset("ICES");
  TH1D * oH_noleak = (TH1D*)hA->Clone("hcombined_noleak");
  oH_noleak->Reset("ICES");
  TH1D * oHC = (TH1D*)hA->Clone("hcombined_C");
  oHC->Reset("ICES");
  for (int i = 0; i < ana::nPtBins; i++) {
    H[i] = new TH1D(Form("bootstrap%i",i),";bootstrapped value; counts",100,-0.2,1.5);
    HC[i] = new TH1D(Form("bootstrapC%i",i),";bootstrapped value; counts",100,-0.2,1.5);
    for (int j = 0; j < 10000; j++) {
      float A = rand->Gaus(hA->GetBinContent(i+1), hA->GetBinError(i+1));
      float B = rand->Gaus(hB->GetBinContent(i+1), hB->GetBinError(i+1));
      float C = rand->Gaus(hC->GetBinContent(i+1), hC->GetBinError(i+1));
      float D = rand->Gaus(hD->GetBinContent(i+1), hD->GetBinError(i+1));
      float a = rand->Gaus(ha->GetBinContent(i+1), ha->GetBinError(i+1));
      float b = rand->Gaus(hb->GetBinContent(i+1), hb->GetBinError(i+1));
      float c = rand->Gaus(hc->GetBinContent(i+1), hc->GetBinError(i+1));
      float d = rand->Gaus(hd->GetBinContent(i+1), hd->GetBinError(i+1));

      float qa = d-b*c;
      float qb = -(A*d+D)+(B*c+C*b);
      float qc = A*D-B*C;

      float S;

      if (A == 0 || B == 0 || C == 0 || D == 0 || fabs(qa) < 1e-10 || qb*qb - 4*qa*qc < 0) continue;

      float Sp = (-qb + TMath::Sqrt(qb*qb - 4*qa*qc))/2/qa;
      float Sm = (-qb - TMath::Sqrt(qb*qb - 4*qa*qc))/2/qa;
      if (Sp < A && Sp > 0) {
        S = Sp;
      }
      else {
        S = Sm;
      }

      H[i]->Fill(S/A);
      H2->Fill(i,S/A);
      // Signal content of region C falls straight out of this same leakage-corrected
      // solve (n_s^C = c*S, by definition of c as the MC leakage fraction of C relative
      // to A) - no independent quadratic/MC template needed for region C's purity.
      if (C != 0) HC[i]->Fill(c*S/C);
    }
    // Non-bootstrap version
    float A = hA->GetBinContent(i+1);
    float B = hB->GetBinContent(i+1);
    float C = hC->GetBinContent(i+1);
    float D = hD->GetBinContent(i+1);
    float a = ha->GetBinContent(i+1);
    float b = hb->GetBinContent(i+1);
    float c = hc->GetBinContent(i+1);
    float d = hd->GetBinContent(i+1);

    float qa = d-b*c;
    float qb = -(A*d+D)+(B*c+C*b);
    float qc = A*D-B*C;

    float S;

    if (A == 0 || B == 0 || C == 0 || D == 0 || fabs(qa) < 1e-10 || qb*qb - 4*qa*qc < 0) continue;

    float Sp = (-qb + TMath::Sqrt(qb*qb - 4*qa*qc))/2/qa;
    float Sm = (-qb - TMath::Sqrt(qb*qb - 4*qa*qc))/2/qa;
    if (Sp < A && Sp > 0) {
      S = Sp;
    }
    else {
      S = Sm;
    }
    oH->SetBinContent(i+1,S/A);
    oH_noleak->SetBinContent(i+1, 1-B*C/A/D);
    oHC->SetBinContent(i+1, c*S/C);

    if (H[i]->GetEntries() > 0) {
      double probs[3] = {0.16, 0.50, 0.84};
      double q[3];
      H[i]->GetQuantiles(3, q, probs);
      cout << q[0] << " " << q[1] << " " << q[2] << endl;

      oh->SetPoint(i, hA->GetBinCenter(i+1),q[1]);
      oh->SetPointError(i, hA->GetBinWidth(i+1)/2.0,hA->GetBinWidth(i+1)/2.0,q[1] - q[0],q[2] - q[1]);
    }
    if (HC[i]->GetEntries() > 0) {
      double probs[3] = {0.16, 0.50, 0.84};
      double qC[3];
      HC[i]->GetQuantiles(3, qC, probs);

      ohC->SetPoint(i, hC->GetBinCenter(i+1),qC[1]);
      ohC->SetPointError(i, hC->GetBinWidth(i+1)/2.0,hC->GetBinWidth(i+1)/2.0,qC[1] - qC[0],qC[2] - qC[1]);
    }
  }

  TF1 * func = new TF1("func","TMath::Erf((x - [1])/[2])",8,100);
  func->SetParameter(0,1);
  func->SetParameter(1,13);
  func->SetParameter(2,5);
  // func itself (fitted here) is what ana::getPurity(val,...) evaluates downstream - no
  // caller reads the fit's own covariance, so the fit isn't asked to return one ("S"
  // dropped from the option string below; used to be captured as a TFitResultPtr and
  // written out as "purityFitResult", but nothing ever read that object back).
  oh->Fit(func,"RIMQ0");

  for (int i = 0; i < ana::nPtBins; i++) {
    H[i]->Write();
    HC[i]->Write();
  }
  H2->Write();
  oh->Write();
  ohC->Write();
  oH->Write();
  oH_noleak->Write();
  oHC->Write();
  func->Write();

  if (graphCOut) *graphCOut = ohC;
  return oh;
}

    
// systag: nominal (default), JERhigh, JERlow, emscale_high, emscale_low, EMRhigh,
// EMRlow, jes_high, jes_low, threejet, narrowBDT, narrowISO, narrowBDTbkg,
// narrowISObkg, wideISObkg - selects which reprocessing of
// hclusterpt_abcd (both Data and the Photon MC leakage fractions)
// this purity curve is derived from. See the unfolder constructor comment in
// src/unfolder.h for what each one means.
//
// Loops every jet radius internally (hclusterpt_abcd%i_%i is already filled per radius,
// gated on ispaired[ir] - see unfolder.cc) and writes all seven into ONE
// ana::purityFilename(systag) file, one ana::rnames[ir] subdirectory per radius - purity
// is "of paired photons", and pairing genuinely differs by jet radius, so it needs its
// own value per radius, not one number reused everywhere (see src/ana.h's getPurity ir
// parameter). drawer/canvases are constructed once and reused/Clear()'d each radius
// rather than rebuilt, since drawer's own file opens and TCanvas's fixed names would
// otherwise be repeated 7x pointlessly (drawer) or warn on collision (TCanvas) within
// one process.
void puritymaker(string systag = "nominal") {
  const char * histname = "hclusterpt_abcd";
  // MC leakage-fraction templates (fp[i] below) are read from the truth-matched (photon
  // deltaR<0.1) subset instead of histname - the method's own definition is
  // f^X=N_sig^X/N_sig^A, a true-signal ratio, not an all-reconstructed-cluster ratio.
  // Data's own ABCD counts (h[] below) still read histname unchanged: Data has no truth
  // info, and its raw ABCD counts are genuinely what the data-driven method has to work
  // with regardless.
  const char * histname_truthmatched = "hclusterpt_abcd_truthmatched";
  drawer d("pythia", systag);
  gStyle->SetOptStat(0);

  TCanvas * cf = new TCanvas("cf","",700,700);
  TCanvas * co = new TCanvas("co","",700,700);
  TCanvas * cu = new TCanvas("cu","",700,700);
  TCanvas * c  = new TCanvas("c","",700,700);

  string purityOutfile = ana::purityFilename(systag);
  TFile * fout = TFile::Open(purityOutfile.c_str(), "RECREATE");

  string purityPdfPath = Form("/home/samson72/sphnx/gammajet_unfold/pdfs/purity_%s.pdf", systag.c_str());
  cu->SaveAs(Form("%s[", purityPdfPath.c_str()));

  // MC leakage-fraction (f^X) plot, one page per radius - previously drawn to canvas
  // `cf` but never saved to disk.
  string leakagePdfPath = Form("/home/samson72/sphnx/gammajet_unfold/pdfs/purity_leakage_%s.pdf", systag.c_str());
  cf->SaveAs(Form("%s[", leakagePdfPath.c_str()));

  for (int ir = 0; ir < ana::nJetR; ir++) {
  TH1D * h[4]; // for ABCD
  TH1D * hp[4];
  TH1D * fp[4];
  for (int i = 0; i < 4; i++) {
    h[i] = d.get(Form("%s%i_%i",histname,ir,i),0);
    hp[i] = d.get(Form("%s%i_%i",histname_truthmatched,ir,i),1);
    fp[i] = (TH1D*)hp[i]->Clone(Form("fp%i",i));
    fp[i]->Divide(hp[i],hp[0]);
  }
  cf->cd();
  cf->Clear();
  gPad->SetTicks();
  gPad->SetLeftMargin(.15);
  // DrawFrame (rather than relying on the first "same"-drawn histogram to create an
  // axis frame) so the axes reliably render on every page of the multi-page save below -
  // same pattern as cu's frameu further down.
  TH1F * framef = cf->DrawFrame(ana::ptBinsUsed[0], 0, ana::ptBinsUsed[ana::nPtBinsUsed], 1.2);
  framef->GetXaxis()->SetTitle("Leading cluster p_{T} [GeV]");
  framef->GetYaxis()->SetTitle("f^{X} = N^{X}_{sig}/N^{A}_{sig}");
  int colors[4] = {kBlack,kRed, kBlue, kOrange};
  const char * letters[4] = {"A","B","C","D"};
  TLegend * lf = new TLegend(0.5, 0.65, 0.8, 0.8);
  for (int i = 1; i < 4; i++) {
    fp[i]->SetLineColor(colors[i]);
    fp[i]->SetLineWidth(2);
    fp[i]->Draw("hist e same");
    lf->AddEntry(fp[i],Form("N_{%s}/N_{A}",letters[i]));
  }
  lf->SetLineWidth(0);
  lf->Draw();
  d.drawAll({"Pythia8 #gamma+jet MC"},{Form("systag: %s",systag.c_str()),Form("Jet R=%.1f",ana::JetRs[ir]),"truth-matched leakage fractions"},.5,.55,16,700);
  cf->SaveAs(leakagePdfPath.c_str());
  //return;


  TGraphAsymmErrors * odC = nullptr;
  fout->cd();
  fout->mkdir(ana::rnames[ir])->cd();
  TGraphAsymmErrors * od = combine_hists(h,fp,&odC);
  fout->cd();

  TF1 * func = new TF1(Form("func_%s",ana::rnames[ir]),"TMath::Erf((x - [1])/[2])",8,30);
  func->SetParameter(0,1);
  func->SetParameter(1,13);
  func->SetParameter(2,5);
  od->Fit(func,"RIMQ0");
  //od->Draw();
  co->cd();
  co->Clear();
  od->SetLineColor(kBlack);
  od->SetMarkerColor(kBlack);
  od->SetMarkerSize(1);
  od->SetMarkerStyle(20);
  od->GetYaxis()->SetRangeUser(0,1.2);
  od->GetXaxis()->SetRangeUser(8,40);
  od->Draw("ap");
  od->GetXaxis()->SetTitle("leading cluster p_{T}");
  func->Draw("same");
  TLegend * lo = new TLegend(.5,.4,.8,.6);
  lo->AddEntry(od,"purity");
  lo->AddEntry(func,"Error function fit");
  lo->SetLineWidth(0);
  lo->Draw();
  d.drawAll({},{Form("systag: %s",systag.c_str()),Form("Jet R=%.1f",ana::JetRs[ir]),"paired clusters","leakage correction applied"},0.15,0.8,20,700);
  d.drawText(Form("P(p_T) = erf((x - %.2f)/%.2f)",func->GetParameter(1), func->GetParameter(2)), .5, .8,1);

  // Purity vs. photon pT, restricted to the ana::nPtBinsUsed reported bins (15-20,
  // 20-25, 25-35 GeV) - `od` above also carries the low-pT migration-only buffer bin
  // (13-15 GeV) and the high-pT overflow bin (35-100 GeV), which are diagnostic only,
  // not reported physics bins (see ana.h's ptBins/ptBinsUsed/firstUsedPtBin comment), so
  // they're excluded here. Same bootstrap points/asymmetric errors as od, just a subset.
  TGraphAsymmErrors * odUsed = new TGraphAsymmErrors(ana::nPtBinsUsed);
  odUsed->SetName("combined_used");
  // Same subset, region C - plotted alongside odUsed below purely as a sanity check
  // that P_C comes out sensible (e.g. much lower than P_A, since C is the background-
  // enriched sideband) before it's used anywhere downstream.
  TGraphAsymmErrors * odCUsed = new TGraphAsymmErrors(ana::nPtBinsUsed);
  odCUsed->SetName("combined_C_used");
  for (int k = 0; k < ana::nPtBinsUsed; k++) {
    int ipt = ana::firstUsedPtBin + k;
    double x, y;
    od->GetPoint(ipt, x, y);
    odUsed->SetPoint(k, x, y);
    odUsed->SetPointError(k, od->GetErrorXlow(ipt), od->GetErrorXhigh(ipt), od->GetErrorYlow(ipt), od->GetErrorYhigh(ipt));
    if (odC) {
      double xc, yc;
      odC->GetPoint(ipt, xc, yc);
      odCUsed->SetPoint(k, xc, yc);
      odCUsed->SetPointError(k, odC->GetErrorXlow(ipt), odC->GetErrorXhigh(ipt), odC->GetErrorYlow(ipt), odC->GetErrorYhigh(ipt));
    }
  }
  cu->cd();
  cu->Clear();
  gPad->SetTicks();
  gPad->SetLeftMargin(.15);
  TH1F * frameu = cu->DrawFrame(ana::ptBinsUsed[0], 0, ana::ptBinsUsed[ana::nPtBinsUsed], 1.1);
  frameu->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  frameu->GetYaxis()->SetTitle("Purity");
  odUsed->SetLineColor(kBlack);
  odUsed->SetMarkerColor(kBlack);
  odUsed->SetMarkerSize(1);
  odUsed->SetMarkerStyle(20);
  odUsed->SetLineWidth(2);
  odUsed->Draw("p same");
  odCUsed->SetLineColor(kAzure+2);
  odCUsed->SetMarkerColor(kAzure+2);
  odCUsed->SetMarkerSize(1);
  odCUsed->SetMarkerStyle(21);
  odCUsed->SetLineWidth(2);
  odCUsed->Draw("p same");
  TLegend * lu = new TLegend(.2,.7,.5,.85);
  lu->SetLineWidth(0);
  lu->AddEntry(odUsed,  "P_{A} (region A)");
  lu->AddEntry(odCUsed, "P_{C} (region C)");
  lu->Draw();
  d.drawAll({"p+p Run24 Data"},{Form("systag: %s",systag.c_str()),Form("Jet R=%.1f",ana::JetRs[ir]),"paired clusters","leakage correction applied"},.18,.3,16,700);
  cu->SaveAs(purityPdfPath.c_str());

  c->cd();
  c->Clear();
  gPad->SetTicks();
  gPad->SetLogy();
  h[0]->SetLineColor(kBlack);
  h[1]->SetLineColor(kBlue);
  h[2]->SetLineColor(kOrange);
  h[3]->SetLineColor(kRed);
  TLegend * l = new TLegend(.5,.4,.8,.6);
  l->SetLineWidth(0);
  string text[4] = {"Region A","Region B","Region C","Region D"};
  for (int i = 0; i < 4; i++) {
    h[i]->SetLineWidth(2);
    h[i]->GetYaxis()->SetRangeUser(0.5,1e5);
    h[i]->GetXaxis()->SetRangeUser(8,40);
    h[i]->Draw("hist same");
    l->AddEntry(h[i],text[i].c_str());
  }
  l->Draw();
  d.drawAll({"p+p Run24"},{Form("systag: %s",systag.c_str()),"Paired clusters"},.5,.8,20,700);
  } // end of ir loop

  cu->SaveAs(Form("%s]", purityPdfPath.c_str()));
  cout << "Wrote " << purityPdfPath << endl;
  cf->SaveAs(Form("%s]", leakagePdfPath.c_str()));
  cout << "Wrote " << leakagePdfPath << endl;
  fout->Close();
  cout << "Wrote " << purityOutfile << endl;
}
