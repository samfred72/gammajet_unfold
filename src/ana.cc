#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
using namespace std;
ana::ana() {
}

Bool_t ana::PassEtaCut(float eta, float vz = 0)
{
  //float loweta = GetShiftedEta(vz,etamin);
  //float higheta = GetShiftedEta(vz,etamax);
  if (eta < etamin || eta > etamax) return false;
  else return true;
}

Double_t ana::GetShiftedEta(float _vz, float _eta)
{
  double theta = 2*atan(exp(-_eta));
  double z = radius / tan(theta);
  double zshifted = z - _vz;
  double thetashifted = atan2(radius,zshifted);
  double etashifted = -log(tan(thetashifted/2.0));
  return etashifted;
}

Int_t ana::findPtBin(double value)
{
  for (int i = 0; i < nPtBins; ++i) {
    if (value >= ptBins[i] && value < ptBins[i + 1]) {
      return  i;
    }
  }
  return -1;
}
Int_t ana::findUnfoldXjBin(double value)
{
  for (int i = 0; i < nUnfoldXjBins; ++i) {
    if (value >= unfoldXjBins[i] && value < unfoldXjBins[i + 1]) {
      return  i;
    }
  }
  return -1;
}
Int_t ana::findxjBin(double value)
{
  for (int i = 0; i < nxjBins; ++i) {
    if (value >= xjBins[i] && value < xjBins[i + 1]) {
      return  i;
    }
  }
  return -1;
}
Int_t ana::findBdtBin(double value)
{
  for (int i = 0; i < nBdtBins; ++i) {
    if (value >= bdtBins[i] && value < bdtBins[i + 1]) {
      return  i;
    }
  }
  return -1;
}

Int_t ana::findUnfoldBin(double xj, double pt)
{
  int ipt = findPtBin(pt);
  int ixj = findUnfoldXjBin(xj);
   
  if (xj < 0) return -1;
  if (ipt < 0) return -1;
  
  if (xj >= 2.0) ixj = nUnfoldXjBins;
 
  return ipt*(nUnfoldXjBins+2) + ixj + 1; // +2 for underflow and overflow bins, +1 for the undeflow bin
}
Int_t ana::findabcdBin(double iso, double bdt, float pt, int bin)
{
  int isiso;
  int isbdt;
  if (iso <= -999) { // PhotonClusterBuilder's iso_topo_valid==0 sentinel: topocluster iso not computed
    isiso = -1;
  }
  else if (iso < isoBins[bin]) {
    isiso = 1;
  }
  else if (iso > isoBinsHigh[bin]) {
    isiso = 0;
  }
  else {
    isiso = -1;
  }
  if (bdt > bdtGoodLow[bin] - 0.00156 * pt && bdt < bdtGoodHigh[bin]) {
    isbdt = 1;
  }
  else if (bdt > bdtBadLow[bin] && bdt < bdtBadHigh[bin]) {
    isbdt = 0;
  }
  else {
    isbdt = -1;
  }

  if (isiso == -1 || isbdt == -1) {
    return -1;
  }
  else {
    bool b_isiso = isiso;
    bool b_isbdt = isbdt;
    int iabcd = (((b_isbdt << 0b1) | b_isiso) ^ 0b11); // silly bitwise operations to map isiso+isbdt->A,B,C,D (index 0,1,2,3)
    return iabcd;
  }
}

Int_t ana::findabcdBin(double iso, double bdt, int bin)
{
  int isiso;
  int isbdt;
  if (iso <= -999) { // PhotonClusterBuilder's iso_topo_valid==0 sentinel: topocluster iso not computed
    isiso = -1;
  }
  else if (iso < isoBins[bin]) {
    isiso = 1;
  }
  else if (iso > isoBinsHigh[bin]) {
    isiso = 0;
  }
  else {
    isiso = -1;
  }
  if (bdt > bdtGoodLow[bin] && bdt < bdtGoodHigh[bin]) {
    isbdt = 1;
  }
  else if (bdt > bdtBadLow[bin] && bdt < bdtBadHigh[bin]) {
    isbdt = 0;
  }
  else {
    isbdt = -1;
  }

  if (isiso == -1 || isbdt == -1) {
    return -1;
  }
  else {
    bool b_isiso = isiso;
    bool b_isbdt = isbdt;
    int iabcd = (((b_isbdt << 0b1) | b_isiso) ^ 0b11); // silly bitwise operations to map isiso+isbdt->A,B,C,D (index 0,1,2,3)
    return iabcd;
  }
}

Int_t ana::findabcdBin(double iso, int showershape, int bin) {
  int isiso;
  int istight;
  if (iso <= -999) { // PhotonClusterBuilder's iso_topo_valid==0 sentinel: topocluster iso not computed
    isiso = -1;
  }
  else if (iso < isoBins[bin]) {
    isiso = 1;
  }
  else if (iso > isoBinsHigh[bin]) {
    isiso = 0;
  }
  else {
    isiso = -1;
  }
  if (showershape == 2) {
    istight = 1;
  }
  else if (showershape == 1) {
    istight = 0;
  }
  else {
    istight = -1;
  }
  if (isiso == -1 || istight == -1) {
    return -1;
  }
  else {
    bool b_isiso = isiso;
    bool b_istight = istight;
    int iabcd = (((b_istight << 0b1) | b_isiso) ^ 0b11); // silly bitwise operations to map isiso+isbdt->A,B,C,D (index 0,1,2,3)
    return iabcd;
  }
}

Int_t ana::findHadronBin(double value) {
  for (int i = 0; i < nHadronBins; i++) {
    float binlow = hadronBins[i][0];
    float binhigh = hadronBins[i][1];
    if (value > binlow && value < binhigh) {
      return i;
    }
  }
  return -1;
}
Int_t ana::findEmfracBin(double value) {
  for (int i = 0; i < nEmfracBins; i++) {
    float binlow = emfracBins[i];
    float binhigh = emfracBins[i+1];
    if (value > binlow && value < binhigh) {
      return i;
    }
  }
  return -1;
}

// low/high are expected to be the edges of one ana::ptBins bin (that's how every caller
// invokes this), so (low+high)/2 lands on the bin center and findPtBin recovers the bin
// index directly - this reads the actual puritymaker.C point/error for that bin rather
// than a smooth fit evaluated/integrated over the range.
float ana::getPurity(float low, float high, string systag) {
  TFile * f = TFile::Open(Form("/home/samson72/sphnx/gammajet_unfold/hists/purity_%s.root", systag.c_str()));
  TGraphAsymmErrors * oh = (TGraphAsymmErrors*)f->Get("combined");
  int ipt = findPtBin((low+high)/2.0);
  if (ipt < 0) {
    f->Close();
    cout << "WARNING: ana::getPurity - (low+high)/2 = " << (low+high)/2.0
         << " doesn't fall within any ana::ptBins bin. Returning 0." << endl;
    return 0;
  }
  double x, y;
  oh->GetPoint(ipt, x, y);
  f->Close();
  return y;
}
float ana::getPurity(float val, string systag) {
  TFile * f = TFile::Open(Form("/home/samson72/sphnx/gammajet_unfold/hists/purity_%s.root", systag.c_str()));
  TF1 * func = (TF1*)f->Get("func");
  float ret = func->Eval(val);
  f->Close();
  return ret;
}
// puritymaker.C's bootstrap errors are asymmetric (16th/84th percentile around the
// median) - keep them that way rather than collapsing to one symmetric number.
float ana::getPurityErrorLow(float low, float high, string systag) {
  TFile * f = TFile::Open(Form("/home/samson72/sphnx/gammajet_unfold/hists/purity_%s.root", systag.c_str()));
  TGraphAsymmErrors * oh = (TGraphAsymmErrors*)f->Get("combined");
  int ipt = findPtBin((low+high)/2.0);
  if (ipt < 0) {
    f->Close();
    cout << "WARNING: ana::getPurityErrorLow - (low+high)/2 = " << (low+high)/2.0
         << " doesn't fall within any ana::ptBins bin. Returning 0." << endl;
    return 0;
  }
  float err = oh->GetErrorYlow(ipt);
  f->Close();
  return err;
}
float ana::getPurityErrorHigh(float low, float high, string systag) {
  TFile * f = TFile::Open(Form("/home/samson72/sphnx/gammajet_unfold/hists/purity_%s.root", systag.c_str()));
  TGraphAsymmErrors * oh = (TGraphAsymmErrors*)f->Get("combined");
  int ipt = findPtBin((low+high)/2.0);
  if (ipt < 0) {
    f->Close();
    cout << "WARNING: ana::getPurityErrorHigh - (low+high)/2 = " << (low+high)/2.0
         << " doesn't fall within any ana::ptBins bin. Returning 0." << endl;
    return 0;
  }
  float err = oh->GetErrorYhigh(ipt);
  f->Close();
  return err;
}
