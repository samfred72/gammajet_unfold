#ifndef DRAWER_H
#define DRAWER_H

#include "ana.h"
#include <string>
#include <vector>
#include "TLatex.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLine.h"
using namespace std;

class drawer {
  public :
    // Reads unfolder.cc's output (<Trigger>_<sim>_<systag>_unfolding.root for MC,
    // <Trigger>_<systag>_unfolding.root for Data).
    drawer(string sim = "pythia", string systag = "nominal") : systag(systag) {
      TFile * f = TFile::Open(Form("%s/hists/Data_%s_unfolding.root", ana::dir(), systag.c_str()));
      TFile * f05_p = TFile::Open(Form("%s/hists/%s_%s_%s_unfolding.root", ana::dir(),"Photon5" , sim.c_str(), systag.c_str()));
      TFile * f10_p = TFile::Open(Form("%s/hists/%s_%s_%s_unfolding.root", ana::dir(),"Photon10", sim.c_str(), systag.c_str()));
      TFile * f20_p = TFile::Open(Form("%s/hists/%s_%s_%s_unfolding.root", ana::dir(),"Photon20", sim.c_str(), systag.c_str()));
      TFile * f08_j = TFile::Open(Form("%s/hists/%s_%s_%s_unfolding.root", ana::dir(),"Jet8"    , sim.c_str(), systag.c_str()));
      TFile * f12_j = TFile::Open(Form("%s/hists/%s_%s_%s_unfolding.root", ana::dir(),"Jet12"   , sim.c_str(), systag.c_str()));
      TFile * f20_j = TFile::Open(Form("%s/hists/%s_%s_%s_unfolding.root", ana::dir(),"Jet20"   , sim.c_str(), systag.c_str()));
      TFile * f30_j = TFile::Open(Form("%s/hists/%s_%s_%s_unfolding.root", ana::dir(),"Jet30"   , sim.c_str(), systag.c_str()));
      TFile * f50_j = TFile::Open(Form("%s/hists/%s_%s_%s_unfolding.root", ana::dir(),"Jet50"   , sim.c_str(), systag.c_str()));
      TFile * f60_j = TFile::Open(Form("%s/hists/%s_%s_%s_unfolding.root", ana::dir(),"Jet60"   , sim.c_str(), systag.c_str()));
      TFile * f80_j = TFile::Open(Form("%s/hists/%s_%s_%s_unfolding.root", ana::dir(),"Jet80"   , sim.c_str(), systag.c_str()));

      dfiles[0] = f;

      pfiles[0] = f05_p;
      pfiles[1] = f10_p;
      pfiles[2] = f20_p;
      
      jfiles[0] = f08_j;
      jfiles[1] = f12_j;
      jfiles[2] = f20_j;
      jfiles[3] = f30_j;
      jfiles[4] = f50_j;
      jfiles[5] = f60_j;
      jfiles[6] = f80_j;
      
      scalemap = (sim == "pythia" ?
        map<bool,map<int,double>> {
          {0,{{5,1.369e+08},{8,1.3013e+07},{12,3.997e+06},{20,6.218e+04},{30,2.502e+03},{50,7.2695},{60,0},{80,0}}},
          {1,{{5,146359.3},{10,6944.675},{20,130.4461}}},
        } :
           map<bool,map<int,double>> {
          {0,{{5,1.369e+08},{8,1.3013e+07},{12,3.997e+06},{20,6.218e+04},{30,2.502e+03},{50,7.2695},{60,0},{80,0}}},
          {1,{{5,6.48487e+05},{10,3.62808e+02},{20,5.34010e+01}}}
        }
      );

    }
    ~drawer(); 

    void drawLine(float x1, float y1, float x2, float y2);
    void drawText(const char *text, float xp, float yp, int textColor=kBlack, int textSize=18);
    void drawAll(vector<string> samples, vector<string> features, float drawx, float drawy, int fontsize, float csize);
    void drawMany(vector<string> features, float drawx, float drawy, int fontsize, float csize);
    void scale(TH1D * h, float low = 0, float high = 2);
    void format(TH1D * h, int type);
    void format(TF1 * h, int type);
    TF1 * fit(TH1D * h, float low, float high, const char * options);
    vector<vector<vector<vector<vector<vector<TH1D*>>>>>> collect_hists(const char * histname, int type);
    vector<vector<vector<vector<vector<vector<TH1D*>>>>>> get_empty_TH1D();
    vector<vector<vector<vector<vector<vector<TF1*>>>>>> get_empty_TF1();
    vector<vector<vector<vector<vector<vector<float>>>>>> get_empty_float();
    TH1D * combine_hists(TH1D * A, TH1D * B, TH1D * C, TH1D * D, int ipt, string name); 
    TH1D * combineMC(const char * histname, bool isphoton); 
    TH2D * combineMC2d(const char * histname, bool isphoton); 
    TH1D * get(const char * histname, int type, int ihist = -1);
    TH2D * get2d(const char * histname, int type, int ihist = -1);
    double getScale(bool isphoton, int sample) { return scalemap[isphoton][sample]; }

  private :

    string systag;
    const static int ndfiles = 1;
    const static int npfiles = 3;
    const static int njfiles = 7;
    vector<TFile*> dfiles{ndfiles};
    vector<TFile*> pfiles{npfiles};
    vector<TFile*> jfiles{njfiles};

    vector<int> psamples = {5,10,20};
    vector<int> jsamples = {8,12,20,30,50,60,80};

    string simulation = "";
    // 0: pythia jets, 1: pythia photon, 2: herwig photon
    map<bool,map<int,double>> scalemap;

    TH1D * empty_hist = new TH1D("empty_hist","",25,0,2);
};

#endif // DRAWER_H
