#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
using namespace std;
ana::ana() {
}

// See ana.h's comment above emResolutionSigma for the prescription/provenance.
static float emSigmaOverE(float E, float p0, float p1, float p2) {
  return sqrt(p0*p0/E + p1*p1/(E*E) + p2*p2);
}
float ana::emResolutionSigma(float truthPt, int emrVariant) {
  if (emrVariant == emrLow || truthPt <= 0) return 0.0;
  const float mc = emSigmaOverE(truthPt, 0.185, 0.0, 0.040);
  const float data = (emrVariant == emrHigh) ? emSigmaOverE(truthPt, 0.13, 0.08, 0.08)
                                             : emSigmaOverE(truthPt, 0.15, 0.05, 0.05);
  return sqrt(max(0.0f, data*data - mc*mc));
}

// See ana.h's comment above the declaration - the single place to add/remove a systag.
const vector<string> ana::systags = {
  "nominal", "JERhigh", "JERlow", "emscale_high", "emscale_low",
  "jes_high", "jes_low", "threejet", "narrowBDT", "narrowISO",
  "EMRhigh", "EMRlow", "narrowBDTbkg", "narrowISObkg", "wideISObkg"
};
const vector<pair<string,string>> ana::asymmetricSystagPairs = {
  {"JERhigh", "JERlow"}, {"emscale_high", "emscale_low"}, {"jes_high", "jes_low"},
  {"EMRhigh", "EMRlow"}
};
// narrowISObkg/wideISObkg are a genuine two-sided variation of the same boundary (see
// isoBinsHigh below) but are deliberately NOT listed above - unlike JER/emscale/jes/EMR,
// they're each treated as their own independent symmetrized source (full magnitude to
// both up and down), not sign-split against each other. See drawing/draw_systematics.C's
// purityMembers comment for why.

double ana::jesForSystag(const string & systag, int ir) {
  for (int i = 0; i < nJesSystags; i++) {
    if (systag == jesSystagNames[i]) return jesBySystag[i][ir];
  }
  cout << "WARNING: ana::jesForSystag - no in-situ p_a for systag \"" << systag
       << "\" in ana.h's jesBySystag (rerun insitu/draw_jes_summary.C); using jesNominal." << endl;
  return jesNominal[ir];
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

string ana::purityFilename(const string & systag) {
  return string(Form("/home/samson72/sphnx/gammajet_unfold/hists/purity_%s.root", systag.c_str()));
}

// getPurity/getPurityC and their ErrorLow/ErrorHigh siblings below used to
// TFile::Open() the same purity_<systag>.root fresh on every single call (and never
// `delete f` after Close() - Close() alone doesn't free the TFile object, only the OS
// file handle, so every call also leaked one small TFile object). Calling these 6
// functions 9 times each (once per pT bin) is negligible at the couple-of-calls-per-
// macro-run rate they were designed for, but a toy bootstrap loop that re-derives the
// full purity-corrected spectrum on every toy (draw_toy_vs_analytic.C's data-side toys,
// via unfold_utility::buildFullyCorrected) calls this 9*6=54 times PER TOY - 540,000
// file-opens (plus 540,000 leaked TFile objects) over a 10,000-toy run, which is real
// I/O and allocation overhead dominating the macro's runtime. The purity value for a
// given (systag, ir) never changes between toys (it only depends on the fixed purity
// curve, not the toyed A/C counts), so cache the open file per systag instead of
// reopening it - correctness is unaffected, every caller just gets the same file back.
static TFile * cachedPurityFile(const string & systag) {
  static map<string, TFile*> cache;
  string fname = ana::purityFilename(systag);
  auto it = cache.find(fname);
  if (it != cache.end() && it->second && !it->second->IsZombie()) return it->second;
  TFile * f = TFile::Open(fname.c_str());
  cache[fname] = f;
  return f;
}

// low/high are expected to be the edges of one ana::ptBins bin (that's how every caller
// invokes this), so (low+high)/2 lands on the bin center and findPtBin recovers the bin
// index directly - this reads the actual puritymaker.C point/error for that bin rather
// than a smooth fit evaluated/integrated over the range.
//
// purityFilename(systag) now holds every jet radius's purity curve in its own
// ana::rnames[ir] subdirectory (see puritymaker.C) - Get() reaches into that
// subdirectory via a "<rname>/objname" path instead of opening a radius-suffixed file.
float ana::getPurity(float low, float high, string systag, int ir) {
  TFile * f = cachedPurityFile(systag);
  TGraphAsymmErrors * oh = (TGraphAsymmErrors*)f->Get(Form("%s/combined", rnames[ir]));
  int ipt = findPtBin((low+high)/2.0);
  if (ipt < 0) {
    cout << "WARNING: ana::getPurity - (low+high)/2 = " << (low+high)/2.0
         << " doesn't fall within any ana::ptBins bin. Returning 0." << endl;
    return 0;
  }
  double x, y;
  oh->GetPoint(ipt, x, y);
  return y;
}
float ana::getPurity(float val, string systag, int ir) {
  TFile * f = cachedPurityFile(systag);
  TF1 * func = (TF1*)f->Get(Form("%s/func", rnames[ir]));
  float ret = func->Eval(val);
  return ret;
}
// puritymaker.C's bootstrap errors are asymmetric (16th/84th percentile around the
// median) - keep them that way rather than collapsing to one symmetric number.
float ana::getPurityErrorLow(float low, float high, string systag, int ir) {
  TFile * f = cachedPurityFile(systag);
  TGraphAsymmErrors * oh = (TGraphAsymmErrors*)f->Get(Form("%s/combined", rnames[ir]));
  int ipt = findPtBin((low+high)/2.0);
  if (ipt < 0) {
    cout << "WARNING: ana::getPurityErrorLow - (low+high)/2 = " << (low+high)/2.0
         << " doesn't fall within any ana::ptBins bin. Returning 0." << endl;
    return 0;
  }
  float err = oh->GetErrorYlow(ipt);
  return err;
}
float ana::getPurityErrorHigh(float low, float high, string systag, int ir) {
  TFile * f = cachedPurityFile(systag);
  TGraphAsymmErrors * oh = (TGraphAsymmErrors*)f->Get(Form("%s/combined", rnames[ir]));
  int ipt = findPtBin((low+high)/2.0);
  if (ipt < 0) {
    cout << "WARNING: ana::getPurityErrorHigh - (low+high)/2 = " << (low+high)/2.0
         << " doesn't fall within any ana::ptBins bin. Returning 0." << endl;
    return 0;
  }
  float err = oh->GetErrorYhigh(ipt);
  return err;
}
// Region-C analogues of getPurity/getPurityErrorLow/getPurityErrorHigh above - same
// puritymaker.C bootstrap, read from the "combined_C" graph instead of "combined".
float ana::getPurityC(float low, float high, string systag, int ir) {
  TFile * f = cachedPurityFile(systag);
  TGraphAsymmErrors * oh = (TGraphAsymmErrors*)f->Get(Form("%s/combined_C", rnames[ir]));
  int ipt = findPtBin((low+high)/2.0);
  if (ipt < 0) {
    cout << "WARNING: ana::getPurityC - (low+high)/2 = " << (low+high)/2.0
         << " doesn't fall within any ana::ptBins bin. Returning 0." << endl;
    return 0;
  }
  double x, y;
  oh->GetPoint(ipt, x, y);
  return y;
}
float ana::getPurityCErrorLow(float low, float high, string systag, int ir) {
  TFile * f = cachedPurityFile(systag);
  TGraphAsymmErrors * oh = (TGraphAsymmErrors*)f->Get(Form("%s/combined_C", rnames[ir]));
  int ipt = findPtBin((low+high)/2.0);
  if (ipt < 0) {
    cout << "WARNING: ana::getPurityCErrorLow - (low+high)/2 = " << (low+high)/2.0
         << " doesn't fall within any ana::ptBins bin. Returning 0." << endl;
    return 0;
  }
  float err = oh->GetErrorYlow(ipt);
  return err;
}
float ana::getPurityCErrorHigh(float low, float high, string systag, int ir) {
  TFile * f = cachedPurityFile(systag);
  TGraphAsymmErrors * oh = (TGraphAsymmErrors*)f->Get(Form("%s/combined_C", rnames[ir]));
  int ipt = findPtBin((low+high)/2.0);
  if (ipt < 0) {
    cout << "WARNING: ana::getPurityCErrorHigh - (low+high)/2 = " << (low+high)/2.0
         << " doesn't fall within any ana::ptBins bin. Returning 0." << endl;
    return 0;
  }
  float err = oh->GetErrorYhigh(ipt);
  return err;
}
