#include "insitu_utility.h"
#include "unfold_utility.h"
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TLatex.h"
using namespace std;

vector<DataEvent> insitu_utility::cacheDataEvents(const char * filename, int abcdSelect, int ir, bool restrictToUsed) {
  vector<DataEvent> events;
  TFile * f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie()) {
    cout << "WARNING: could not open " << filename << endl;
    return events;
  }
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
    if (abcd != abcdSelect) continue;
    if (evIr != ir) continue;
    int ipt = ana::findPtBin(pho_pt);
    if (restrictToUsed) {
      if (ipt < ana::firstUsedPtBin || ipt >= ana::firstUsedPtBin + ana::nPtBinsUsed) continue;
      ipt -= ana::firstUsedPtBin;
    }
    else {
      if (ipt < 0) continue;
    }
    events.push_back({pho_pt, jet_pt, ipt, third_pt});
  }
  f->Close();
  return events;
}

double insitu_utility::lowXjFloor(int ir, double ptLow) {
  double lowval = ana::jet_calib_pt_cut[ir]/ptLow;
  return ana::unfoldXjBins[ana::findUnfoldXjBin(lowval)+1];
}

string insitu_utility::insituFilename(const char * insitu_dir, const char * trigger,
    const char * sim, const string & systag) {
  if (sim && sim[0] != '\0') {
    return string(Form("%s/%s_%s_%s_insitu.root", insitu_dir, trigger, sim, systag.c_str()));
  }
  return string(Form("%s/%s_%s_insitu.root", insitu_dir, trigger, systag.c_str()));
}

void insitu_utility::findError(TGraph * g, int ibest, float minchisq, float & errLow, float & errHigh) {
  double xbest, ytmp;
  g->GetPoint(ibest, xbest, ytmp);
  errLow = xbest - g->GetX()[0];
  errHigh = g->GetX()[g->GetN()-1] - xbest;
  double x, y;
  for (int i = ibest; i >= 0; i--) {
    g->GetPoint(i, x, y);
    if (y - minchisq > 1.0) { errLow = xbest - x; break; }
  }
  for (int i = ibest; i < g->GetN(); i++) {
    g->GetPoint(i, x, y);
    if (y - minchisq > 1.0) { errHigh = x - xbest; break; }
  }
}

TGraphErrors * insitu_utility::meanGraph(const float mean[], const float err[], const char * name) {
  TGraphErrors * g = new TGraphErrors(ana::nPtBinsUsed);
  g->SetName(name);
  for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) {
    float lo = ana::ptBinsUsed[ipt], hi = ana::ptBinsUsed[ipt+1];
    g->SetPoint(ipt, (lo+hi)/2.0, mean[ipt]);
    g->SetPointError(ipt, (hi-lo)/2.0, err[ipt]);
  }
  return g;
}

TGraphErrors * insitu_utility::ratioGraph(const float meanNum[], const float errNum[],
    const float meanDen[], const float errDen[], const char * name) {
  TGraphErrors * g = new TGraphErrors(ana::nPtBinsUsed);
  g->SetName(name);
  for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) {
    float lo = ana::ptBinsUsed[ipt], hi = ana::ptBinsUsed[ipt+1];
    float num = meanNum[ipt], den = meanDen[ipt];
    if (num <= 0 || den <= 0) { g->SetPoint(ipt, (lo+hi)/2.0, 0); g->SetPointError(ipt, (hi-lo)/2.0, 0); continue; }
    float ratio = num/den;
    float err = ratio*sqrt(pow(errNum[ipt]/num,2) + pow(errDen[ipt]/den,2));
    g->SetPoint(ipt, (lo+hi)/2.0, ratio);
    g->SetPointError(ipt, (hi-lo)/2.0, err);
  }
  return g;
}

void insitu_utility::drawSPhenixLabel(vector<string> samples, vector<string> features, float drawx, float drawy, int fontsize, float csize) {
  float titlescale = 1.25;
  float subtitlescale = 1.25;
  float ydiff = fontsize * 0.0017 * 700.0/csize;
  auto drawOne = [&](const char * text, float xp, float yp, int size) {
    TLatex * tex = new TLatex(xp, yp, text);
    tex->SetTextFont(43);
    tex->SetTextSize(size);
    tex->SetTextColor(kBlack);
    tex->SetLineWidth(1);
    tex->SetNDC();
    tex->Draw();
  };
  drawOne("#bf{#it{sPHENIX}} #kern[0.5]{Internal}", drawx, drawy, (int)(fontsize*titlescale));
  for (unsigned i = 0; i < samples.size(); i++) {
    drawOne(samples[i].c_str(), drawx, drawy-ydiff*subtitlescale*(i+1), (int)(fontsize*subtitlescale));
  }
  for (unsigned i = 0; i < features.size(); i++) {
    drawOne(features[i].c_str(), drawx, drawy-ydiff*subtitlescale*samples.size()-ydiff*(i+1)*subtitlescale, fontsize);
  }
}

// ----- Shared helpers for the grid_insitu*.C JES scans -----

vector<double> insitu_utility::coarsenSum(const vector<double> & fine, int startBin, int nBinsForFit, int groupSize) {
  int nGroups = (nBinsForFit - startBin + groupSize - 1) / groupSize;
  vector<double> coarse(nGroups, 0.);
  for (int ixj = startBin; ixj < nBinsForFit; ixj++) coarse[(ixj - startBin)/groupSize] += fine[ixj];
  return coarse;
}

vector<double> insitu_utility::coarsenQuadrature(const vector<double> & fineErr, int startBin, int nBinsForFit, int groupSize) {
  int nGroups = (nBinsForFit - startBin + groupSize - 1) / groupSize;
  vector<double> coarse(nGroups, 0.);
  for (int ixj = startBin; ixj < nBinsForFit; ixj++) coarse[(ixj - startBin)/groupSize] += fineErr[ixj]*fineErr[ixj];
  for (auto & v : coarse) v = sqrt(v);
  return coarse;
}

void insitu_utility::referenceMeans(const vector<pair<string,double>> & samples, int abcdSelect, int ir,
    float refMean[], float refMeanErr[], const float lowXj[]) {
  vector<double> sumw(ana::nPtBinsUsed,0), sumw2(ana::nPtBinsUsed,0), sumwx(ana::nPtBinsUsed,0), sumwx2(ana::nPtBinsUsed,0);
  for (auto & s : samples) {
    TFile * f = TFile::Open(s.first.c_str(), "READ");
    if (!f || f->IsZombie()) {
      cout << "WARNING: could not open " << s.first << endl;
      continue;
    }
    TTree * t = (TTree*)f->Get("insitutree");
    Float_t pho_pt, jet_pt, mcWeight;
    Int_t abcd, evIr;
    t->SetBranchAddress("pho_pt", &pho_pt);
    t->SetBranchAddress("jet_pt", &jet_pt);
    t->SetBranchAddress("abcd", &abcd);
    t->SetBranchAddress("weight", &mcWeight);
    t->SetBranchAddress("ir", &evIr);
    Long64_t nentries = t->GetEntries();
    for (Long64_t e = 0; e < nentries; e++) {
      t->GetEntry(e);
      if (abcd != abcdSelect) continue;
      if (evIr != ir) continue;
      int ipt = ana::findPtBin(pho_pt);
      if (ipt < ana::firstUsedPtBin || ipt >= ana::firstUsedPtBin + ana::nPtBinsUsed) continue;
      ipt -= ana::firstUsedPtBin;
      double x = jet_pt/pho_pt;
      if (x < lowXj[ipt]) continue;
      double w = s.second*mcWeight;
      sumw[ipt]   += w;
      sumw2[ipt]  += w*w;
      sumwx[ipt]  += w*x;
      sumwx2[ipt] += w*x*x;
    }
    f->Close();
  }
  for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) {
    if (sumw[ipt] <= 0) { refMean[ipt] = 0; refMeanErr[ipt] = 0; continue; }
    double mean = sumwx[ipt]/sumw[ipt];
    double var  = sumwx2[ipt]/sumw[ipt] - mean*mean;
    double neff = sumw[ipt]*sumw[ipt]/sumw2[ipt]; // Kish effective sample size
    refMean[ipt] = mean;
    refMeanErr[ipt] = sqrt(std::max(var,0.)/neff);
  }
}

void insitu_utility::referenceShape(const vector<pair<string,double>> & samples, int abcdSelect, int ir,
    vector<vector<double>> & refFrac, vector<vector<double>> & refFracErr, const float lowXj[]) {
  vector<vector<double>> sumw(ana::nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
  vector<vector<double>> sumw2(ana::nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
  for (auto & s : samples) {
    TFile * f = TFile::Open(s.first.c_str(), "READ");
    if (!f || f->IsZombie()) {
      cout << "WARNING: could not open " << s.first << endl;
      continue;
    }
    TTree * t = (TTree*)f->Get("insitutree");
    Float_t pho_pt, jet_pt, mcWeight;
    Int_t abcd, evIr;
    t->SetBranchAddress("pho_pt", &pho_pt);
    t->SetBranchAddress("jet_pt", &jet_pt);
    t->SetBranchAddress("abcd", &abcd);
    t->SetBranchAddress("weight", &mcWeight);
    t->SetBranchAddress("ir", &evIr);
    Long64_t nentries = t->GetEntries();
    for (Long64_t e = 0; e < nentries; e++) {
      t->GetEntry(e);
      if (abcd != abcdSelect) continue;
      if (evIr != ir) continue;
      int ipt = ana::findPtBin(pho_pt);
      if (ipt < ana::firstUsedPtBin || ipt >= ana::firstUsedPtBin + ana::nPtBinsUsed) continue;
      ipt -= ana::firstUsedPtBin;
      double x = jet_pt/pho_pt;
      if (x < lowXj[ipt]) continue;
      int ixj = ana::findUnfoldXjBin(x);
      if (ixj < 0 || ixj >= ana::nUnfoldXjBins) continue;
      double w = s.second*mcWeight;
      sumw[ipt][ixj]  += w;
      sumw2[ipt][ixj] += w*w;
    }
    f->Close();
  }
  refFrac.assign(ana::nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
  refFracErr.assign(ana::nPtBinsUsed, vector<double>(ana::nUnfoldXjBins, 0.));
  for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) {
    double N = 0;
    for (double w : sumw[ipt]) N += w;
    if (N <= 0) continue;
    for (int ixj = 0; ixj < ana::nUnfoldXjBins; ixj++) {
      refFrac[ipt][ixj]    = sumw[ipt][ixj]/N;
      refFracErr[ipt][ixj] = sqrt(sumw2[ipt][ixj])/N;
    }
  }
}

void insitu_utility::computeRegionAMeans(const vector<DataEvent> & dataA, float pa, float mean[], float err[], const float lowXj[]) {
  vector<double> sum(ana::nPtBinsUsed,0), sum2(ana::nPtBinsUsed,0);
  vector<int> count(ana::nPtBinsUsed,0);
  for (auto & ev : dataA) {
    if (vetoed(ev, pa)) continue;
    float x = (ev.jet_pt/pa)/ev.pho_pt;
    if (x < lowXj[ev.ptbin]) continue;
    sum[ev.ptbin]  += x;
    sum2[ev.ptbin] += x*x;
    count[ev.ptbin]++;
  }
  for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) {
    if (count[ipt] == 0) { mean[ipt] = 0; err[ipt] = 0; continue; }
    double m   = sum[ipt]/count[ipt];
    double var = sum2[ipt]/count[ipt] - m*m;
    mean[ipt] = m;
    err[ipt]  = sqrt(std::max(var,0.)/count[ipt]);
  }
}

void insitu_utility::computeCorrectedMeans(const vector<DataEvent> & dataA, const vector<DataEvent> & dataC,
    float pa, const float purity[], const float purityC[], float mean[], float err[], const float lowXj[]) {
  vector<double> sumA(ana::nPtBinsUsed,0), sumA2(ana::nPtBinsUsed,0);
  vector<int> countA(ana::nPtBinsUsed,0);
  vector<double> sumC(ana::nPtBinsUsed,0), sumC2(ana::nPtBinsUsed,0);
  vector<int> countC(ana::nPtBinsUsed,0);
  for (auto & ev : dataA) {
    if (vetoed(ev, pa)) continue;
    float x = (ev.jet_pt/pa)/ev.pho_pt;
    if (x < lowXj[ev.ptbin]) continue;
    sumA[ev.ptbin] += x; sumA2[ev.ptbin] += x*x; countA[ev.ptbin]++;
  }
  for (auto & ev : dataC) {
    if (vetoed(ev, pa)) continue;
    float x = (ev.jet_pt/pa)/ev.pho_pt;
    if (x < lowXj[ev.ptbin]) continue;
    sumC[ev.ptbin] += x; sumC2[ev.ptbin] += x*x; countC[ev.ptbin]++;
  }
  for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) {
    if (countA[ipt] == 0 || countC[ipt] == 0) { mean[ipt] = 0; err[ipt] = 0; continue; }
    double NA = countA[ipt], NC = countC[ipt];
    float coeffA, coeffC;
    unfold_utility::purityCorrectCoeffs(purity[ipt], purityC[ipt], NA, NC, coeffA, coeffC);
    double sumXcorr  = coeffA*sumA[ipt]  - coeffC*sumC[ipt];
    double sumX2corr = coeffA*sumA2[ipt] - coeffC*sumC2[ipt];
    double Ncorr = coeffA*NA - coeffC*NC;
    if (Ncorr <= 0) { mean[ipt] = 0; err[ipt] = 0; continue; }
    double m   = sumXcorr/Ncorr;
    double var = sumX2corr/Ncorr - m*m;
    mean[ipt] = m;
    err[ipt]  = sqrt(std::max(var,0.)/Ncorr);
  }
}

vector<TH1D*> insitu_utility::buildXjByPtBin(const vector<DataEvent> & data, float pa, int nBins,
    const char * prefix, const float lowXj[]) {
  vector<TH1D*> h(nBins);
  for (int ipt = 0; ipt < nBins; ipt++) {
    h[ipt] = new TH1D(Form("%s_pt%d", prefix, ipt), ";x_{J#gamma};Counts", ana::nUnfoldXjBins, ana::unfoldXjBins);
  }
  for (auto & ev : data) {
    if (vetoed(ev, pa)) continue;
    float x = (ev.jet_pt/pa)/ev.pho_pt;
    if (x < lowXj[ev.ptbin]) continue;
    h[ev.ptbin]->Fill(x);
  }
  return h;
}

vector<TH1D*> insitu_utility::buildMCXjByPtBin(const vector<pair<string,double>> & samples, int abcdSelect, int ir,
    const char * prefix, const float lowXj[]) {
  vector<TH1D*> h(ana::nPtBinsUsed);
  for (int ipt = 0; ipt < ana::nPtBinsUsed; ipt++) {
    h[ipt] = new TH1D(Form("%s_pt%d", prefix, ipt), ";x_{J#gamma};Counts", ana::nUnfoldXjBins, ana::unfoldXjBins);
  }
  for (auto & s : samples) {
    TFile * f = TFile::Open(s.first.c_str(), "READ");
    if (!f || f->IsZombie()) continue;
    TTree * t = (TTree*)f->Get("insitutree");
    Float_t pho_pt, jet_pt, mcWeight; Int_t abcd, evIr;
    t->SetBranchAddress("pho_pt", &pho_pt);
    t->SetBranchAddress("jet_pt", &jet_pt);
    t->SetBranchAddress("abcd", &abcd);
    t->SetBranchAddress("weight", &mcWeight);
    t->SetBranchAddress("ir", &evIr);
    Long64_t nentries = t->GetEntries();
    for (Long64_t e = 0; e < nentries; e++) {
      t->GetEntry(e);
      if (abcd != abcdSelect) continue;
      if (evIr != ir) continue;
      int ipt = ana::findPtBin(pho_pt);
      if (ipt < ana::firstUsedPtBin || ipt >= ana::firstUsedPtBin + ana::nPtBinsUsed) continue;
      ipt -= ana::firstUsedPtBin;
      float x = jet_pt/pho_pt;
      if (x < lowXj[ipt]) continue;
      h[ipt]->Fill(x, s.second*mcWeight);
    }
    f->Close();
  }
  return h;
}

vector<TH1D*> insitu_utility::purityCorrectByPtBin(const vector<TH1D*> & hA, const vector<TH1D*> & hC, int nBins,
    const float purity[], const float purityErrLow[], const float purityErrHigh[],
    const float purityC[], const float purityCErrLow[], const float purityCErrHigh[],
    const char * prefix) {
  vector<TH1D*> h(nBins);
  for (int ipt = 0; ipt < nBins; ipt++) {
    // Last pT bin (35-100 GeV migration buffer) has an empty region B.
    bool quiet = (ipt == ana::nPtBins - 1);
    TH1D * hcorr = unfold_utility::purityCorrect(hA[ipt], hC[ipt],
        purity[ipt], purityErrLow[ipt], purityErrHigh[ipt],
        purityC[ipt], purityCErrLow[ipt], purityCErrHigh[ipt],
        Form("%s_pt%d", prefix, ipt), nullptr, quiet);
    h[ipt] = hcorr ? hcorr : (TH1D*)hA[ipt]->Clone(Form("%s_pt%d", prefix, ipt));
  }
  return h;
}

// ----- Multijet balance -----

int insitu_utility::findMultijetPtBin(double leadPt) {
  for (int i = 0; i < nMultijetPtBins; i++)
    if (leadPt >= multijetPtBins[i] && leadPt < multijetPtBins[i+1]) return i;
  return -1;
}

string insitu_utility::multijetAnalysisFilename(const char * multijet_dir, const char * sim) {
  return Form("%s/multijet_analysis_%s.root", multijet_dir, sim);
}

string insitu_utility::multijetSysName(const string & systag) {
  if (systag == "JERhigh") return "HIGH";
  if (systag == "JERlow")  return "LOW";
  return "RECO";
}

vector<MultijetEvent> insitu_utility::cacheMultijetEvents(const char * filename, const string & treename) {
  vector<MultijetEvent> events;
  TFile * f = TFile::Open(filename, "READ");
  if (!f || f->IsZombie()) {
    cout << "WARNING: could not open multijet analysis file " << filename << endl;
    return events;
  }
  TTree * t = (TTree*)f->Get(treename.c_str());
  if (!t) {
    cout << "WARNING: no tree " << treename << " in " << filename << endl;
    f->Close();
    return events;
  }
  float lead, sl, slphi, ssl, sslphi, w;
  t->SetBranchAddress("leadingPT", &lead);
  t->SetBranchAddress("SLPT", &sl);
  t->SetBranchAddress("SLphi", &slphi);
  t->SetBranchAddress("SSLPT", &ssl);
  t->SetBranchAddress("SSLphi", &sslphi);
  t->SetBranchAddress("weight", &w);
  Long64_t n = t->GetEntries();
  for (Long64_t e = 0; e < n; e++) {
    t->GetEntry(e);
    events.push_back({lead, sl, slphi, ssl, sslphi, w});
  }
  f->Close();
  return events;
}

double insitu_utility::multijetBalance(const MultijetEvent & ev, double pa, double pb, double & leadCorr, bool applyCuts) {
  double lead = ev.lead / (pa + pb*ev.lead);
  double s1 = ev.sl / (pa + pb*ev.sl), s2 = ev.ssl / (pa + pb*ev.ssl);
  leadCorr = lead;
  double px = s1*std::cos(ev.slphi) + s2*std::cos(ev.sslphi);
  double py = s1*std::sin(ev.slphi) + s2*std::sin(ev.sslphi);
  double recoil = std::sqrt(px*px + py*py);
  if (applyCuts && (lead < multijetLeadCut || s1 < multijetRecoilJetCut || s2 < multijetRecoilJetCut ||
                    recoil < multijetRecoilCut)) return -1;
  return recoil > 0 ? lead/recoil : -1;
}

void insitu_utility::multijetMeans(const vector<MultijetEvent> & events, double pa, double pb,
    float mean[], float err[], bool applyCuts) {
  vector<double> sw(nMultijetPtBins,0), sw2(nMultijetPtBins,0), swx(nMultijetPtBins,0), swx2(nMultijetPtBins,0);
  for (auto & ev : events) {
    double lead;
    double b = multijetBalance(ev, pa, pb, lead, applyCuts);
    int bin = findMultijetPtBin(lead);
    if (bin < 0 || b < multijetBalanceLow || b >= multijetBalanceHigh) continue;
    sw[bin] += ev.w; sw2[bin] += ev.w*ev.w; swx[bin] += ev.w*b; swx2[bin] += ev.w*b*b;
  }
  for (int i = 0; i < nMultijetPtBins; i++) {
    if (sw[i] <= 0) { mean[i] = 0; err[i] = 0; continue; }
    double m = swx[i]/sw[i];
    double var = swx2[i]/sw[i] - m*m;
    double neff = sw[i]*sw[i]/sw2[i];
    mean[i] = m;
    err[i] = std::sqrt(std::max(var, 0.)/neff);
  }
}
