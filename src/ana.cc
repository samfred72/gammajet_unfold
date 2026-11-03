#include "ana.h"
#include <cstdlib>
#include <TFile.h>
#include <TH1.h>
#include <deque>
#include <fstream>
#include <map>
#include <sstream>
using namespace std;

const char * ana::dir() {
  static const string d = [] {
    const char * e = getenv("GAMMAJET_UNFOLD");
    if (!e || !*e) { cerr << "GAMMAJET_UNFOLD is not set (see README.md)" << endl; exit(1); }
    return string(e);
  }();
  return d.c_str();
}

bool ana::paperLabels() {
  const char * e = getenv("GJ_PAPER_LABELS");
  return e && string(e) == "1";
}

const char * ana::sphenixLabel() {
  static string label;
  const char * e = getenv("GJ_LABEL_STATUS");
  string status = e ? string(e) : (paperLabels() ? "" : "Internal");
  label = "#bf{#it{sPHENIX}}";
  if (!status.empty()) label += " #kern[0.5]{" + status + "}";
  return label.c_str();
}

ana::JerMode ana::mcJerMode() {
  static const JerMode mode = [] {
    const char * e = getenv("GAMMAJET_JER");
    string m = (e && *e) ? string(e) : "tree";
    return m == "tree" ? jerTree : m == "none" ? jerNone : jerFile;
  }();
  return mode;
}

float ana::mcJetPt(int ir, float calib, float treeSmeared, float z, float truthmatch) {
  if (mcJerMode() == jerTree) return treeSmeared;
  if (mcJerMode() == jerNone) return calib;
  static const vector<TH1*> sigma = [] {
    const char * e = getenv("GAMMAJET_JER_HIST");
    string pattern = (e && *e) ? string(e) : "h_jer_smear_r04_pileup_EMfracJES_nominal";
    TFile f(getenv("GAMMAJET_JER"));
    if (f.IsZombie()) { cerr << "GAMMAJET_JER=" << getenv("GAMMAJET_JER") << ": cannot open" << endl; exit(1); }
    vector<TH1*> h(nJetR, nullptr);
    vector<string> names(nJetR);
    for (int jr = 0; jr < nJetR; jr++) {
      names[jr] = pattern;
      size_t pos = names[jr].find("{R}");
      if (pos != string::npos) names[jr].replace(pos, 3, Form("%02d", (int)lround(10*JetRs[jr])));
      h[jr] = (TH1*)f.Get(names[jr].c_str());
      if (h[jr]) { h[jr] = (TH1*)h[jr]->Clone(); h[jr]->SetDirectory(nullptr); }
    }
    vector<TH1*> out(nJetR, nullptr);
    for (int jr = 0; jr < nJetR; jr++) {
      int use = -1;
      for (int kr = 0; kr < nJetR; kr++)
        if (h[kr] && (use < 0 || fabs(JetRs[kr] - JetRs[jr]) < fabs(JetRs[use] - JetRs[jr]))) use = kr;
      if (use < 0) { cerr << "GAMMAJET_JER=" << getenv("GAMMAJET_JER") << ": no TH1 " << pattern << endl; exit(1); }
      out[jr] = h[use];
      cout << "MC JER templates R=" << JetRs[jr] << ": " << getenv("GAMMAJET_JER") << " : " << names[use]
           << (use != jr ? "  (WARNING: no " + names[jr] + ", using the nearest radius)" : "") << endl;
    }
    return out;
  }();
  const float ref = truthmatch > 0 ? truthmatch : calib;
  return calib + z * ref * sigma[ir]->Interpolate(ref);
}

const char * ana::dataLabel() { return paperLabels() ? "p+p #sqrt{s}=200 GeV" : "p+p Run24 Data"; }

const char * ana::mcLabel() { return paperLabels() ? "Pythia8 #gamma+jet" : "Pythia8 #gamma+jet MC"; }

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
  "EMRhigh", "EMRlow", "narrowBDTbkg", "narrowISObkg", "wideISObkg",
  "timingwide"
};
const vector<pair<string,string>> ana::asymmetricSystagPairs = {
  {"JERhigh", "JERlow"}, {"emscale_high", "emscale_low"}, {"jes_high", "jes_low"},
  {"EMRhigh", "EMRlow"}
};
// narrowISObkg/wideISObkg are symmetrized independently, not sign-split (see draw_systematics.C).

// GAMMAJET_JES_TABLE: optional text file replacing the compiled p_a table, one line per systag,
// "<systag> <p_a R=0.2> ... <p_a R=0.8>" ('#' starts a comment). Read once.
static const map<string, vector<double>> & jesTableOverride() {
  static const map<string, vector<double>> table = [] {
    map<string, vector<double>> t;
    const char * e = getenv("GAMMAJET_JES_TABLE");
    if (!e || !*e) return t;
    ifstream f(e);
    if (!f) { cerr << "GAMMAJET_JES_TABLE=" << e << ": cannot open" << endl; exit(1); }
    string line;
    while (getline(f, line)) {
      line = line.substr(0, line.find('#'));
      istringstream iss(line);
      string tag;
      if (!(iss >> tag)) continue;
      vector<double> v(ana::nJetR);
      for (auto & x : v) if (!(iss >> x)) { cerr << "GAMMAJET_JES_TABLE=" << e << ": short row " << tag << endl; exit(1); }
      t[tag] = v;
    }
    if (!t.count("nominal")) { cerr << "GAMMAJET_JES_TABLE=" << e << ": no nominal row" << endl; exit(1); }
    cout << "ana: in-situ p_a from GAMMAJET_JES_TABLE=" << e << " (" << t.size() << " systags), not ana.h" << endl;
    return t;
  }();
  return table;
}

double ana::jesForSystag(const string & systag, int ir) {
  const auto & ov = jesTableOverride();
  if (!ov.empty()) {
    auto it = ov.find(systag);
    if (it != ov.end()) return it->second[ir];
    static bool warned = false;
    if (!warned) cout << "WARNING: ana::jesForSystag - systag \"" << systag << "\" not in GAMMAJET_JES_TABLE; using its nominal row." << endl;
    warned = true;
    return ov.at("nominal")[ir];
  }
  for (int i = 0; i < nJesSystags; i++) {
    if (systag == jesSystagNames[i]) return jesBySystag[i][ir];
  }
  cout << "WARNING: ana::jesForSystag - no in-situ p_a for systag \"" << systag
       << "\" in ana.h's jesBySystag (rerun insitu/compute_jes_summary.C); using jesNominal." << endl;
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
