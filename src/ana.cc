#include "ana.h"
#include <cstdlib>
#include <deque>
using namespace std;

const char * ana::dir() {
  static const string d = [] {
    const char * e = getenv("GAMMAJET_UNFOLD");
    if (!e || !*e) { cerr << "GAMMAJET_UNFOLD is not set (see README.md)" << endl; exit(1); }
    return string(e);
  }();
  return d.c_str();
}

const char * ana::path(const string & rel) {
  static deque<string> store; // push_back never moves existing elements, so returned pointers stay valid
  store.push_back(string(dir()) + "/" + rel);
  return store.back().c_str();
}
ana::ana() {
}

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

const vector<string> ana::systags = {
  "nominal", "JERhigh", "JERlow", "emscale_high", "emscale_low",
  "jes_high", "jes_low", "threejet", "narrowBDT", "narrowISO",
  "EMRhigh", "EMRlow", "narrowBDTbkg", "narrowISObkg", "wideISObkg"
};
const vector<pair<string,string>> ana::asymmetricSystagPairs = {
  {"JERhigh", "JERlow"}, {"emscale_high", "emscale_low"}, {"jes_high", "jes_low"},
  {"EMRhigh", "EMRlow"}
};
// narrowISObkg/wideISObkg are symmetrized independently, not sign-split (see draw_systematics.C).

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
 
  return ipt*(nUnfoldXjBins+2) + ixj + 1; // +2 for underflow and overflow
}
Int_t ana::findabcdBin(double iso, double bdt, int bin)
{
  int isiso;
  int isbdt;
  if (iso <= -999) { // iso_topo_valid==0: topo iso not computed
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
    int iabcd = (((b_isbdt << 0b1) | b_isiso) ^ 0b11); // isiso+isbdt -> A,B,C,D (0..3)
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
  return string(Form("%s/hists/purity_%s.root", ana::dir(), systag.c_str()));
}

// Cache the purity file per systag: toy loops call these thousands of times.
// threejet at R = 0.2 uses the nominal purity: there the veto removes the non-isolated regions B/D about
// twice as often as A/C (claude_checks/insitu_variations/threejet_dr.C, threejet_runaway.C), which leaves
// the 25-35 GeV purity at 0.66 against 0.76-0.89 for every other variation.
static TFile * cachedPurityFile(const string & systag, int ir) {
  static map<string, TFile*> cache;
  string fname = ana::purityFilename((systag == "threejet" && ir == 0) ? "nominal" : systag);
  auto it = cache.find(fname);
  if (it != cache.end() && it->second && !it->second->IsZombie()) return it->second;
  TFile * f = TFile::Open(fname.c_str());
  cache[fname] = f;
  return f;
}

// low/high are one ana::ptBins bin's edges, so the center gives the bin index. Each radius's curves
// live in the ana::rnames[ir] subdirectory.
float ana::getPurity(float low, float high, string systag, int ir) {
  TFile * f = cachedPurityFile(systag, ir);
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
  TFile * f = cachedPurityFile(systag, ir);
  TF1 * func = (TF1*)f->Get(Form("%s/func", rnames[ir]));
  float ret = func->Eval(val);
  return ret;
}
// Bootstrap errors are asymmetric (16th/84th percentiles); keep them so.
float ana::getPurityErrorLow(float low, float high, string systag, int ir) {
  TFile * f = cachedPurityFile(systag, ir);
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
  TFile * f = cachedPurityFile(systag, ir);
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
// Region-C analogues of the above ("combined_C" graph).
float ana::getPurityC(float low, float high, string systag, int ir) {
  TFile * f = cachedPurityFile(systag, ir);
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
  TFile * f = cachedPurityFile(systag, ir);
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
  TFile * f = cachedPurityFile(systag, ir);
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
