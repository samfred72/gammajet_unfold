#include "../src/drawer.h"
#include "../src/ana.h"
#include "../src/purity_utility.h"
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);


    
// Builds the purity curves for one systag, every jet radius (pairing differs per radius), into
// ana::purityFilename(systag), one ana::rnames[ir] subdirectory each.
void puritymaker(string systag = "nominal") {
  const char * histname = "hclusterpt_abcd";
  // MC leakage fractions f^X = N_sig^X/N_sig^A use the truth-matched subset; Data uses all clusters.
  const char * histname_truthmatched = "hclusterpt_abcd_truthmatched";
  drawer d("pythia", systag);
  gStyle->SetOptStat(0);

  TCanvas * cf = new TCanvas("cf","",700,700);
  TCanvas * co = new TCanvas("co","",700,700);
  TCanvas * cu = new TCanvas("cu","",700,700);
  TCanvas * c  = new TCanvas("c","",700,700);

  string purityOutfile = ana::purityFilename(systag);
  TFile * fout = TFile::Open(purityOutfile.c_str(), "RECREATE");

  string purityPdfPath = Form("%s/pdfs/purity_%s.pdf", ana::dir(), systag.c_str());
  cu->SaveAs(Form("%s[", purityPdfPath.c_str()));

  // MC leakage-fraction plot, one page per radius.
  string leakagePdfPath = Form("%s/pdfs/purity_leakage_%s.pdf", ana::dir(), systag.c_str());
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

  TGraphAsymmErrors * odC = nullptr;
  fout->cd();
  fout->mkdir(ana::rnames[ir])->cd();
  TGraphAsymmErrors * od = purity_utility::combine(h,fp,&odC); // writes into the radius subdirectory
  fout->cd();

  TF1 * func = new TF1(Form("func_%s",ana::rnames[ir]),"TMath::Erf((x - [1])/[2])",8,30);
  func->SetParameter(0,1);
  func->SetParameter(1,13);
  func->SetParameter(2,5);
  od->Fit(func,"RIMQ0");
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

  // Purity vs. photon pT, reported bins only (ana::nPtBinsUsed).
  TGraphAsymmErrors * odUsed = new TGraphAsymmErrors(ana::nPtBinsUsed);
  odUsed->SetName("combined_used");
  // Region C alongside, as a sanity check (P_C well below P_A).
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
  // Label in the empty band between P_C and P_A.
  d.drawAll({"p+p Run24 Data"},{Form("systag: %s",systag.c_str()),Form("Jet R=%.1f",ana::JetRs[ir]),"paired clusters","leakage correction applied"},.55,.6,16,700);
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
