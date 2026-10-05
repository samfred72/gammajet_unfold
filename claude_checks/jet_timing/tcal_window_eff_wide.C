#include "../../src/ana.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);
// Run interpreted: root -l -b -q tcal_window_eff_wide.C - same as tcal_window_eff.C with +-5/7/10 ns windows
// Window efficiency on t_MBD - t_jet (fixed half-width around each version's overall median),
// vs jet pT and EM fraction: standard, corrected, corrected with |t_tower| > 9 ns towers dropped.
void tcal_eff_wide() {
  TChain T("T"); T.Add(ana::path("claude_checks/jet_timing/timingana/tcal/*.root"));
  float vz, mbd, pt, eta, emfrac;
  std::vector<int> *tw_calo=0,*em_status=0; std::vector<float> *tw_e=0,*tw_t=0,*tw_ts=0,*em_e=0,*em_t=0,*em_ts=0; std::vector<bool> *em_good=0;
  T.SetBranchAddress("vz",&vz); T.SetBranchAddress("mbd_time",&mbd); T.SetBranchAddress("jet_pt_calib",&pt); T.SetBranchAddress("jet_eta",&eta); T.SetBranchAddress("jet_emfrac",&emfrac);
  T.SetBranchAddress("tw_calo",&tw_calo); T.SetBranchAddress("tw_e",&tw_e); T.SetBranchAddress("tw_time",&tw_t); T.SetBranchAddress("tw_time_std",&tw_ts);
  T.SetBranchAddress("em_e",&em_e); T.SetBranchAddress("em_time",&em_t); T.SetBranchAddress("em_time_std",&em_ts); T.SetBranchAddress("em_isgood",&em_good); T.SetBranchAddress("em_status",&em_status);
  const int NV=3; const char* vn[NV]={"std","corr","corr+rej"};
  struct J { float pt, ef, d[NV]; }; std::vector<J> jets;
  long nrej=0, ntw=0;
  for (Long64_t i=0;i<T.GetEntries();i++){ T.GetEntry(i);
    if(!std::isfinite(mbd)||!std::isfinite(vz)||fabs(vz)>60||fabs(eta)>0.7) continue;
    float jpt=pt/0.9265; if(jpt<5) continue;
    double tw[NV]={0,0,0}, te[NV]={0,0,0};
    auto add=[&](float e,float tc,float ts){ if(e<=0.5||!std::isfinite(tc)||!std::isfinite(ts)) return; ntw++;
      tw[0]+=e*ts; te[0]+=e; tw[1]+=e*tc; te[1]+=e; if(fabs(tc)<=9){tw[2]+=e*tc; te[2]+=e;} else nrej++; };
    for(size_t j=0;j<em_e->size();j++) if(em_good->at(j)&&!((em_status->at(j)>>5)&1)) add(em_e->at(j),em_t->at(j),em_ts->at(j));
    for(size_t j=0;j<tw_e->size();j++) if(tw_calo->at(j)>0) add(tw_e->at(j),tw_t->at(j),tw_ts->at(j));
    if(te[0]<=0||te[2]<=0) continue;
    J j; j.pt=jpt; j.ef=emfrac; for(int v=0;v<NV;v++) j.d[v]=mbd-tw[v]/te[v]; jets.push_back(j);
  }
  printf("jets %zu, towers rejected by |t|>9 ns: %.2f%%\n", jets.size(), 100.*nrej/ntw);
  double med[NV];
  for(int v=0;v<NV;v++){ std::vector<double> x; for(auto&j:jets) x.push_back(j.d[v]); std::sort(x.begin(),x.end()); med[v]=x[x.size()/2];
    double lo=x[size_t(0.16*x.size())], hi=x[size_t(0.84*x.size())]; printf("%-9s median %6.2f  half 16-84 %.2f ns\n",vn[v],med[v],0.5*(hi-lo)); }
  double ptb[]={5,7,9,11,13,15,18,22,28,40}; double efb[]={0,0.2,0.4,0.6,0.8,1.0001};
  for(double hw: {5.0,7.0,10.0}){
    printf("\nwindow |t_MBD - t_jet - median| < %.0f ns: efficiency\n   pT:     ",hw); for(int b=0;b<9;b++) printf(" %5.0f",ptb[b]); printf("   | EMfrac: "); for(int b=0;b<5;b++) printf(" %4.1f",efb[b]); printf("\n");
    for(int v=0;v<NV;v++){ printf("  %-9s",vn[v]);
      for(int b=0;b<9;b++){ int n=0,k=0; for(auto&j:jets) if(j.pt>=ptb[b]&&j.pt<ptb[b+1]){n++; if(fabs(j.d[v]-med[v])<hw)k++;} printf(" %5.3f",n?double(k)/n:-1); }
      printf("   |        ");
      for(int b=0;b<5;b++){ int n=0,k=0; for(auto&j:jets) if(j.ef>=efb[b]&&j.ef<efb[b+1]){n++; if(fabs(j.d[v]-med[v])<hw)k++;} printf(" %4.2f",n?double(k)/n:-1); }
      printf("\n"); }
  }
}
