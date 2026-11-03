// Writes the in-situ JES for the JER iteration hand-off: the nominal p_a per jet radius from
// insitu/output/grid_insitu_nominal.root, with its statistical errors.
// Convention: Data jet pT is corrected as pT_corr = jet_pt_calib / p_a (unfolder.cc).
// Objects: h_pa (one bin per radius, bin label "R=0.x"; bin error = mean of the -/+ stat errors),
// g_pa (p_a vs R, asymmetric stat errors), TNamed "provenance" (MC JER used, git commit, date).
#include "../src/ana.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);

void write_jes_file(string outname, string jerDescription) {
  TFile fin(ana::path("insitu/output/grid_insitu_nominal.root"));
  TFile fout(outname.c_str(), "RECREATE");
  TH1D * h = new TH1D("h_pa", ";Jet R;p_{a} (pT_{corr} = pT_{calib}/p_{a})", ana::nJetR, 0, ana::nJetR);
  TGraphAsymmErrors * g = new TGraphAsymmErrors(ana::nJetR);
  g->SetName("g_pa");
  g->SetTitle(";Jet R;p_{a}");
  for (int ir = 0; ir < ana::nJetR; ir++) {
    TTree * t = (TTree*)fin.Get(Form("%s/results", ana::rnames[ir]));
    float pa, lo, hi;
    t->SetBranchAddress("pa_puritycorrected", &pa);
    t->SetBranchAddress("errLow_puritycorrected", &lo);
    t->SetBranchAddress("errHigh_puritycorrected", &hi);
    t->GetEntry(0);
    h->GetXaxis()->SetBinLabel(ir+1, Form("R=%.1f", ana::JetRs[ir]));
    h->SetBinContent(ir+1, pa);
    h->SetBinError(ir+1, 0.5*(lo+hi));
    g->SetPoint(ir, ana::JetRs[ir], pa);
    g->SetPointError(ir, 0, 0, lo, hi);
    printf("R=%.1f  p_a = %.4f  -%.4f/+%.4f (stat)\n", ana::JetRs[ir], pa, lo, hi);
  }
  // gSystem->pwd() is JER_iter/ (run_iter.sh runs this macro from there): main's checkout
  string git = gSystem->GetFromPipe("git rev-parse --short HEAD; git diff --quiet HEAD -- ../src ../insitu ../macros ../reweight || echo '(dirty)'").Data();
  string date = gSystem->GetFromPipe("date -Iseconds").Data();
  const char * table = gSystem->Getenv("GAMMAJET_JES_TABLE");
  TNamed prov("provenance", Form("in-situ gamma+jet JES, nominal only\nMC JER: %s\ngit: %s\nlast JES table in: %s\ndate: %s\n"
                                 "pT_corr = jet_pt_calib / p_a for Data jets",
                                 jerDescription.c_str(), git.c_str(), table ? table : "ana.h", date.c_str()));
  fout.cd();
  h->Write();
  g->Write();
  prov.Write();
  fout.Close();
  cout << "wrote " << outname << endl;
}
