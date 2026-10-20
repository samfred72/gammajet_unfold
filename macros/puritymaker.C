#include "../src/drawer.h"
#include "../src/ana.h"
#include "../src/purity_utility.h"
// Load explicitly: the sibling gammajet project's libgammajet.so has same-named classes.
R__LOAD_LIBRARY(libgammajet_unfold.so);


    
// Builds the purity curves for one systag, every jet radius (pairing differs per radius), into
// ana::purityFilename(systag), one ana::rnames[ir] subdirectory each. Calculation only:
// drawing/draw_purity.C draws the purity and leakage-fraction pages from this file.
void puritymaker(string systag = "nominal") {
  const char * histname = "hclusterpt_abcd";
  // MC leakage fractions f^X = N_sig^X/N_sig^A use the truth-matched subset; Data uses all clusters.
  const char * histname_truthmatched = "hclusterpt_abcd_truthmatched";
  drawer d("pythia", systag);

  string purityOutfile = ana::purityFilename(systag);
  TFile * fout = TFile::Open(purityOutfile.c_str(), "RECREATE");

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
  TGraphAsymmErrors * odC = nullptr;
  fout->cd();
  TDirectory * rdir = fout->mkdir(ana::rnames[ir]);
  rdir->cd();
  TGraphAsymmErrors * od = purity_utility::combine(h,fp,&odC); // writes into the radius subdirectory
  // MC leakage fractions f^B, f^C, f^D, for drawing/draw_purity.C.
  rdir->cd();
  const char * letters[4] = {"A","B","C","D"};
  for (int i = 1; i < 4; i++) fp[i]->Write(Form("leakage_%s", letters[i]));

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
  rdir->cd();
  odUsed->Write();
  odCUsed->Write();
  fout->cd();
  } // end of ir loop

  fout->Close();
  cout << "Wrote " << purityOutfile << endl;
}
