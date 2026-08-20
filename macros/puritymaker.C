#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "TFitResult.h"
// The original gammajet project builds its OWN, differently-laid-out drawer/ana classes
// into /home/samson72/root/lib/libgammajet.so, sitting on the same library search path
// as this project's libgammajet_unfold.so. Without forcing which one loads first, ROOT's
// implicit symbol autoload can bind drawer/ana calls to the wrong (mismatched-layout)
// library and segfault - explicit load removes the ambiguity. See unfold.C for the same
// pattern; it never had this problem because it already did this.
R__LOAD_LIBRARY(libgammajet_unfold.so);

TGraphAsymmErrors * combine_hists(TH1D * h[], TH1D * f[], const char * outFile) {
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
  TH2D * H2 = new TH2D("bootstrap2D",";bin number;bootstrapped value",ana::nPtBins,ana::ptBins,100,-0.2,1.5);
  TGraphAsymmErrors * oh = new TGraphAsymmErrors(hA->GetNbinsX());
  oh->SetName("combined");
  TH1D * oH = (TH1D*)hA->Clone("hcombined");
  oH->Reset("ICES");
  TH1D * oH_noleak = (TH1D*)hA->Clone("hcombined_noleak");
  oH_noleak->Reset("ICES");
  for (int i = 0; i < ana::nPtBins; i++) {
    H[i] = new TH1D(Form("bootstrap%i",i),";bootstrapped value; counts",100,-0.2,1.5);
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

    if (H[i]->GetEntries() > 0) {
      double probs[3] = {0.16, 0.50, 0.84};
      double q[3];
      H[i]->GetQuantiles(3, q, probs);
      cout << q[0] << " " << q[1] << " " << q[2] << endl;

      oh->SetPoint(i, hA->GetBinCenter(i+1),q[1]);
      oh->SetPointError(i, hA->GetBinWidth(i+1)/2.0,hA->GetBinWidth(i+1)/2.0,q[1] - q[0],q[2] - q[1]);
    }
  }
  
  TF1 * func = new TF1("func","TMath::Erf((x - [1])/[2])",8,100);
  func->SetParameter(0,1);
  func->SetParameter(1,13);
  func->SetParameter(2,5);
  // "S" stores the full fit result (incl. covariance matrix) so ana::getPurityError
  // can propagate the fit uncertainty into the purity-corrected plots downstream.
  TFitResultPtr fitResult = oh->Fit(func,"RIMQS0");

  TFile * of = TFile::Open(outFile,"RECREATE");
  for (int i = 0; i < ana::nPtBins; i++) {
    H[i]->Write();
  }
  H2->Write();
  oh->Write();
  oH->Write();
  oH_noleak->Write();
  func->Write();
  fitResult->Write("purityFitResult");

  return oh;
}

    
// systag: nominal (default), JERhigh, JERlow, emscale_high, emscale_low, jes_high,
// jes_low, threejet, narrowBDT, narrowISO - selects which reprocessing of
// hclusterpt_abcd (both Data and the Photon MC leakage fractions)
// this purity curve is derived from. See the unfolder constructor comment in
// src/unfolder.h for what each one means. Output: hists/purity_<systag>.root (matching
// the <Trigger>_..._<systag>_unfolding.root convention the unfolding stage uses).
void puritymaker(string systag = "nominal") {
  int ir = 2;
  const char * histname = "hclusterpt_abcd";
  drawer d("pythia", systag);
  TH1D * h[4]; // for ABCD
  TH1D * hp[4];
  TH1D * fp[4];
  for (int i = 0; i < 4; i++) {
    h[i] = d.get(Form("%s%i_%i",histname,ir,i),0);
    hp[i] = d.get(Form("%s%i_%i",histname,ir,i),1);
    fp[i] = (TH1D*)hp[i]->Clone(Form("fp%i",i));
    fp[i]->Divide(hp[i],hp[0]);
  }
  TCanvas * cf = new TCanvas("cf","",700,700);
  int colors[4] = {kBlack,kRed, kBlue, kOrange};
  const char * letters[4] = {"A","B","C","D"};
  TLegend * lf = new TLegend(0.5, 0.65, 0.8, 0.8);
  for (int i = 1; i < 4; i++) {
    fp[i]->GetXaxis()->SetRangeUser(10,35);
    fp[i]->GetYaxis()->SetRangeUser(0,1.2);
    fp[i]->SetLineColor(colors[i]);
    fp[i]->SetLineWidth(2);
    fp[i]->Draw("hist e same");
    lf->AddEntry(fp[i],Form("N_{%s}/N_{A}",letters[i]));
  }
  lf->SetLineWidth(0);
  lf->Draw();
  //return;


  TGraphAsymmErrors * od = combine_hists(h,fp,Form("/home/samson72/sphnx/gammajet_unfold/hists/purity_%s.root", systag.c_str()));
  
  TF1 * func = new TF1("func","TMath::Erf((x - [1])/[2])",8,30);
  func->SetParameter(0,1);
  func->SetParameter(1,13);
  func->SetParameter(2,5);
  od->Fit(func,"RIMQ0");
  //od->Draw();
  TCanvas * co = new TCanvas("co","",700,700);
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
  d.drawAll({},{Form("systag: %s",systag.c_str()),"R=0.4","paired clusters","leakage correction applied"},0.15,0.8,20,700);
  d.drawText(Form("P(p_T) = erf((x - %.2f)/%.2f)",func->GetParameter(1), func->GetParameter(2)), .5, .8,1);

  // Purity vs. photon pT, restricted to the ana::nPtBinsUsed reported bins (15-20,
  // 20-25, 25-35 GeV) - `od` above also carries the low-pT migration-only buffer bin
  // (13-15 GeV) and the high-pT overflow bin (35-100 GeV), which are diagnostic only,
  // not reported physics bins (see ana.h's ptBins/ptBinsUsed/firstUsedPtBin comment), so
  // they're excluded here. Same bootstrap points/asymmetric errors as od, just a subset.
  TGraphAsymmErrors * odUsed = new TGraphAsymmErrors(ana::nPtBinsUsed);
  odUsed->SetName("combined_used");
  for (int k = 0; k < ana::nPtBinsUsed; k++) {
    int ipt = ana::firstUsedPtBin + k;
    double x, y;
    od->GetPoint(ipt, x, y);
    odUsed->SetPoint(k, x, y);
    odUsed->SetPointError(k, od->GetErrorXlow(ipt), od->GetErrorXhigh(ipt), od->GetErrorYlow(ipt), od->GetErrorYhigh(ipt));
  }
  TCanvas * cu = new TCanvas("cu","",700,700);
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
  d.drawAll({"p+p Run24 Data"},{Form("systag: %s",systag.c_str()),"R=0.4","paired clusters","leakage correction applied"},.18,.3,16,700);
  const char * purityPdfPath = Form("/home/samson72/sphnx/gammajet_unfold/pdfs/purity_%s.pdf", systag.c_str());
  cu->SaveAs(purityPdfPath);
  cout << "Wrote " << purityPdfPath << endl;

  //TFile * fout = TFile::Open("hists/purity.root","RECREATE");
  //od->Write();
  //func->Write();
  TCanvas * c = new TCanvas("c","",700,700);
  gPad->SetTicks();
  gPad->SetLogy();
  gStyle->SetOptStat(0);
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
}
