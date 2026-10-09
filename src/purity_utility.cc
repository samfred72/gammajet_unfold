#include "purity_utility.h"
#include "insitu_utility.h"
#include "drawer.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TMath.h"
#include <iostream>

TGraphAsymmErrors * purity_utility::combine(TH1D * h[], TH1D * f[], TGraphAsymmErrors ** graphCOut, bool write) {
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
      // Region-C signal: n_s^C = c*S.
      if (C != 0) HC[i]->Fill(c*S/C);
    }
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
      if (write) cout << q[0] << " " << q[1] << " " << q[2] << endl;

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

  if (write) {
    TF1 * func = new TF1("func","TMath::Erf((x - [1])/[2])",8,100);
    func->SetParameter(0,1);
    func->SetParameter(1,13);
    func->SetParameter(2,5);
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
  } else {
    for (int i = 0; i < ana::nPtBins; i++) { delete H[i]; delete HC[i]; }
    delete H2; delete oH; delete oH_noleak; delete oHC;
  }
  delete rand;

  if (graphCOut) *graphCOut = ohC;
  return oh;
}

void purity_utility::leakageFractions(const string & systag, int ir, TH1D * f[]) {
  drawer d("pythia", systag);
  TH1D * hp[4];
  for (int i = 0; i < 4; i++) hp[i] = d.get(Form("hclusterpt_abcd_truthmatched%i_%i",ir,i),1);
  for (int i = 0; i < 4; i++) {
    f[i] = (TH1D*)hp[i]->Clone(Form("fp_%s_%d_%i", systag.c_str(), ir, i));
    f[i]->SetDirectory(nullptr);
    f[i]->Divide(hp[i],hp[0]);
  }
}

void purity_utility::dataCountsFromTree(const string & filename, int ir, float pa, TH1D * h[]) {
  for (int j = 0; j < 4; j++) h[j]->Reset("ICES");
  TFile * f = TFile::Open(filename.c_str(), "READ");
  if (!f || f->IsZombie()) { cout << "WARNING: purity_utility - could not open " << filename << endl; return; }
  TTree * t = (TTree*)f->Get("insitutree");
  Float_t pho_pt, jet_pt, third_pt = -1;
  Int_t abcd, evIr;
  t->SetBranchAddress("pho_pt", &pho_pt);
  t->SetBranchAddress("jet_pt", &jet_pt);
  t->SetBranchAddress("abcd", &abcd);
  t->SetBranchAddress("ir", &evIr);
  if (t->GetBranch("thirdjet_pt")) t->SetBranchAddress("thirdjet_pt", &third_pt);
  Long64_t nentries = t->GetEntries();
  for (Long64_t e = 0; e < nentries; e++) {
    t->GetEntry(e);
    if (evIr != ir || abcd < 0 || abcd > 3) continue;
    float recoJetPt = jet_pt / pa;                              // unfolder.cc: jet_pt_calib / p_a
    if (recoJetPt <= ana::jet_calib_pt_cut[ir]) continue;       // ispaired jet-pT pre-filter
    int ptbin = ana::findPtBin(pho_pt);
    if (ptbin == -1) continue;
    float lowbin = insitu_utility::lowXjFloor(ir, ana::ptBins[ptbin]);  // float, as check_pair compares
    if (recoJetPt/pho_pt < lowbin) continue;
    if (third_pt >= 0 && third_pt/pa > ana::thirdJetPtCut) continue;                       // deferred threejet veto
    h[abcd]->Fill(pho_pt);
  }
  f->Close();
}
