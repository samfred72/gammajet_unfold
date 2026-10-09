#include "../../src/ana.h"
#include "../../src/insitu_utility.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Why does the Jet12_long region-A in-situ cross-check (grid_insitu_jet12.C) give p_a ~2% below the
// region-A-alone fit against the Photon reference? Both fits use raw Data region A; they differ only in
// the MC reference: Photon5/10/20 region A (nearly all signal) vs Jet12_long region A (signal + fake
// photons). Per radius and reported pT bin: Jet12 region-A truth-matched fraction vs the Data purity,
// and reco <x_J> (in-situ trees, low-x_J floor) of each sample/region, Data divided by nominal p_a.
const char * in = ana::path("insitu/inputs");
const char * hists = ana::path("hists");
const char * out = ana::path("insitu/output");
map<int,double> phoXsec = {{5,146359.3},{10,6944.675},{20,130.4461}};

// mean x_J per reported pT bin, region abcd, from files (weight * scale), Data jets / pa
void means(vector<pair<string,double>> files, int abcd, int ir, double pa, double m[], double e[], double n[]) {
  int np = ana::nPtBinsUsed; vector<double> sw(np,0), swx(np,0), swx2(np,0), sw2(np,0);
  for (auto & fs : files) {
    TFile * f = TFile::Open(fs.first.c_str()); TTree * t = (TTree*)f->Get("insitutree");
    Float_t pho, jet, w; Int_t a, r;
    t->SetBranchAddress("pho_pt",&pho); t->SetBranchAddress("jet_pt",&jet); t->SetBranchAddress("abcd",&a); t->SetBranchAddress("weight",&w); t->SetBranchAddress("ir",&r);
    for (Long64_t i = 0; i < t->GetEntries(); i++) {
      t->GetEntry(i); if (a != abcd || r != ir) continue;
      int ip = ana::findPtBin(pho) - ana::firstUsedPtBin; if (ip < 0 || ip >= np) continue;
      double x = jet/pa/pho; if (x < insitu_utility::lowXjFloor(ir, ana::ptBinsUsed[ip])) continue;
      double ww = w*fs.second; sw[ip] += ww; sw2[ip] += ww*ww; swx[ip] += ww*x; swx2[ip] += ww*x*x;
    }
    f->Close();
  }
  for (int ip = 0; ip < np; ip++) { m[ip] = swx[ip]/sw[ip]; double var = swx2[ip]/sw[ip] - m[ip]*m[ip]; double neff = sw[ip]*sw[ip]/sw2[ip]; e[ip] = sqrt(var/neff); n[ip] = neff; }
}
double readPa(const char * file, int ir, const char * br) {
  TFile * f = TFile::Open(Form("%s/%s", out, file)); TTree * t = (TTree*)f->Get(Form("%s/results", ana::rnames[ir]));
  float v; t->SetBranchAddress(br, &v); t->GetEntry(0); f->Close(); return v;
}
void jet12_composition() {
  const int np = ana::nPtBinsUsed;
  vector<pair<string,double>> pho; for (int s : {5,10,20}) pho.push_back({insitu_utility::insituFilename(in, Form("Photon%d",s), "pythia", "nominal"), phoXsec[s]});
  vector<pair<string,double>> j12 = {{insitu_utility::insituFilename(in, "Jet12_long", "pythia", "nominal"), 1.0}};
  vector<pair<string,double>> dat = {{insitu_utility::insituFilename(in, "Data", "", "nominal"), 1.0}};
  printf("%-4s %-6s | %6s %6s | %7s %7s %7s %7s %7s | %7s %7s\n", "R", "pT", "P_data", "f_J12", "sigPho", "J12_A", "J12_C", "DataA", "DataC", "pred", "J12/sig");
  for (int ir = 0; ir < ana::nJetR; ir++) {
    double paNom = readPa("grid_insitu_nominal.root", ir, "pa_puritycorrected");
    double mS[3], eS[3], nS[3], mJA[3], eJA[3], nJA[3], mJC[3], eJC[3], nJC[3], mDA[3], eDA[3], nDA[3], mDC[3], eDC[3], nDC[3];
    means(pho, 0, ir, 1.0, mS, eS, nS); means(j12, 0, ir, 1.0, mJA, eJA, nJA); means(j12, 2, ir, 1.0, mJC, eJC, nJC);
    means(dat, 0, ir, paNom, mDA, eDA, nDA); means(dat, 2, ir, paNom, mDC, eDC, nDC);
    TFile * f = TFile::Open(Form("%s/Jet12_long_pythia_nominal_unfolding.root", hists));
    TH1D * all = (TH1D*)f->Get(Form("hclusterpt_abcd%d_0", ir)), * mat = (TH1D*)f->Get(Form("hclusterpt_abcd_truthmatched%d_0", ir));
    double sumR = 0, sumW = 0;
    for (int ip = 0; ip < np; ip++) {
      double lo = ana::ptBinsUsed[ip], hi = ana::ptBinsUsed[ip+1];
      int b1 = all->FindBin(lo + 1e-3), b2 = all->FindBin(hi - 1e-3);
      double fJ = mat->Integral(b1, b2)/all->Integral(b1, b2);
      double Pd = ana::getPurity(lo, hi, "nominal", ir);
      // Jet12 region A predicted from its matched fraction, signal = Photon ref, background = Jet12 region C
      double pred = fJ*mS[ip] + (1-fJ)*mJC[ip];
      printf("%-4.1f %2.0f-%-3.0f | %6.3f %6.3f | %7.4f %7.4f %7.4f %7.4f %7.4f | %7.4f %7.4f+-%.4f\n", ana::JetRs[ir], lo, hi, Pd, fJ,
          mS[ip], mJA[ip], mJC[ip], mDA[ip], mDC[ip], pred, mJA[ip]/mS[ip], (mJA[ip]/mS[ip])*hypot(eJA[ip]/mJA[ip], eS[ip]/mS[ip]));
      double w = 1/pow(eJA[ip]/mJA[ip], 2); sumR += w*mJA[ip]/mS[ip]; sumW += w;
    }
    f->Close();
    double paJ = readPa("grid_insitu_jet12_nominal.root", ir, "pa_regionA_jet12ref"), paA = readPa("grid_insitu_nominal.root", ir, "pa_regionA");
    printf("     -> <x_J> J12/Photon (weighted over pT) %.4f ; p_a Photon-ref regionA %.4f, Jet12-ref %.4f, ratio Jet12/Photon %.4f ; purity-corrected nominal %.4f\n",
        sumR/sumW, paA, paJ, paJ/paA, paNom);
  }
}
