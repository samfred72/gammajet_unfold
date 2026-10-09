#ifndef UNFOLDER_H
#define UNFOLDER_H

#include "ana.h"
#include "object.h"
#include "pho_object.h"
#include "jet_object.h"
#include "treeuser.h"
#include "reweight_utility.h"
#include <string>
#include <vector>
#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TMarker.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "RooUnfoldTUnfold.h"
using namespace std;

class unfolder : public treeuser {
  public:
    // systags: variations of the same input tree, all filled in one pass (fill_matrix reads each entry
    // once and fills every systag's histogram set). MC only: JERhigh/low, emscale_high/low, EMRhigh/low.
    // Data only: jes_high/low, timingwide. Both: threejet and the five ABCD-boundary variations. sim selects a
    // different input tree (pythia/herwig).
    // All systags share histogram names, so TH1::AddDirectory(kFALSE) is required.
    unfolder(string trigger, string sim, vector<string> systags) : treeuser(trigger, sim), systags(systags) {
      disableBranchesUnusedByUnfolder();
      gErrorIgnoreLevel = kWarning;
      TH1::AddDirectory(kFALSE);
      // MC histograms are filled with the vz/cluster-pT weight, so Sumw2 everywhere.
      TH1::SetDefaultSumw2();
      int nsys = (int)systags.size();

      pho_response.resize(nsys, vector<RooUnfoldResponse*>(ana::nJetR));
      jet_response.resize(nsys, vector<RooUnfoldResponse*>(ana::nJetR));
      pho_response_half.resize(nsys, vector<RooUnfoldResponse*>(ana::nJetR));
      jet_response_half.resize(nsys, vector<RooUnfoldResponse*>(ana::nJetR));

      hphodr.resize(nsys, vector<TH1D*>(ana::nJetR));
      hjetdr.resize(nsys, vector<TH1D*>(ana::nJetR));

      hphoeffnum.resize(nsys, vector<TH1D*>(ana::nJetR));
      hphoeffden.resize(nsys, vector<TH1D*>(ana::nJetR));
      hphoeff.resize(nsys, vector<TEfficiency*>(ana::nJetR));
      hjeteffnum.resize(nsys, vector<TH1D*>(ana::nJetR));
      hjeteffden.resize(nsys, vector<TH1D*>(ana::nJetR));
      hjeteff.resize(nsys, vector<TEfficiency*>(ana::nJetR));

      hphopurnum.resize(nsys, vector<TH1D*>(ana::nJetR));
      hphopurden.resize(nsys, vector<TH1D*>(ana::nJetR));
      hphopur.resize(nsys, vector<TEfficiency*>(ana::nJetR));
      hjetpurnum.resize(nsys, vector<TH1D*>(ana::nJetR));
      hjetpurden.resize(nsys, vector<TH1D*>(ana::nJetR));
      hjetpur.resize(nsys, vector<TEfficiency*>(ana::nJetR));

      hpaireffnum.resize(nsys, vector<TH1D*>(ana::nJetR));
      hpaireffden.resize(nsys, vector<TH1D*>(ana::nJetR));
      hpaireff.resize(nsys, vector<TEfficiency*>(ana::nJetR));
      hpairpurnum.resize(nsys, vector<TH1D*>(ana::nJetR));
      hpairpurden.resize(nsys, vector<TH1D*>(ana::nJetR));
      hpairpur.resize(nsys, vector<TEfficiency*>(ana::nJetR));

      hphomissfake.resize(nsys, vector<TH2D*>(ana::nJetR));
      hjetmissfake.resize(nsys, vector<TH2D*>(ana::nJetR));
      hpairmissfake.resize(nsys, vector<TH2D*>(ana::nJetR));

      hrecojetpt.resize(nsys, vector<TH1D*>(ana::nJetR));
      htruthjetpt.resize(nsys, vector<TH1D*>(ana::nJetR));
      hjetresponse.resize(nsys, vector<TH2D*>(ana::nJetR));
      hrecophopt.resize(nsys, vector<TH1D*>(ana::nJetR));
      htruthphopt.resize(nsys, vector<TH1D*>(ana::nJetR));
      hphoresponse.resize(nsys, vector<TH2D*>(ana::nJetR));

      hrecojetpt_half.resize(nsys, vector<TH1D*>(ana::nJetR));
      htruthjetpt_half.resize(nsys, vector<TH1D*>(ana::nJetR));
      hjetresponse_half.resize(nsys, vector<TH2D*>(ana::nJetR));
      hrecophopt_half.resize(nsys, vector<TH1D*>(ana::nJetR));
      htruthphopt_half.resize(nsys, vector<TH1D*>(ana::nJetR));
      hphoresponse_half.resize(nsys, vector<TH2D*>(ana::nJetR));

      jet_response2D.resize(nsys, vector<RooUnfoldResponse*>(ana::nJetR));
      jet_response_half2D.resize(nsys, vector<RooUnfoldResponse*>(ana::nJetR));

      hrecoxj_abcd.resize(nsys, vector<vector<TH1D*>>(ana::nJetR, vector<TH1D*>(4)));
      htruthxj_abcd.resize(nsys, vector<vector<TH1D*>>(ana::nJetR, vector<TH1D*>(4)));
      hclusterpt_abcd.resize(nsys, vector<vector<TH1D*>>(ana::nJetR, vector<TH1D*>(4)));
      hclusterpt_abcd_truthmatched.resize(nsys, vector<vector<TH1D*>>(ana::nJetR, vector<TH1D*>(4)));
      hrecoxj.resize(nsys, vector<TH1D*>(ana::nJetR));
      htruthxj.resize(nsys, vector<TH1D*>(ana::nJetR));
      hxjresponse.resize(nsys, vector<TH2D*>(ana::nJetR));

      hrecoxj_half.resize(nsys, vector<TH1D*>(ana::nJetR));
      htruthxj_half.resize(nsys, vector<TH1D*>(ana::nJetR));
      hxjresponse_half.resize(nsys, vector<TH2D*>(ana::nJetR));

      hpurity_num.resize(nsys);
      hpurity_den.resize(nsys);
      hpurity_num_1D.resize(nsys);
      hpurity_den_1D.resize(nsys);
      hphoIDeff_bdt.resize(nsys);
      hphoIDeff_iso.resize(nsys);

      insitu_file.resize(nsys);
      insitu_tree.resize(nsys);
      insitu_pho_pt.resize(nsys);
      insitu_jet_pt.resize(nsys);
      insitu_abcd.resize(nsys);
      insitu_weight.resize(nsys);
      insitu_ir.resize(nsys);
      insitu_third_pt.resize(nsys);

      for (int isys = 0; isys < nsys; isys++) {
        for (int i = 0; i < ana::nJetR; i++) {
          hphodr[isys][i] = new TH1D(Form("hphodr%i",i),";#DeltaR_{truth,reco};counts", 100,0,0.4);
          hjetdr[isys][i] = new TH1D(Form("hjetdr%i",i),";#DeltaR_{truth,reco};counts",100,0,0.4);

          hphoeffnum[isys][i] = new TH1D(Form("hphoeffnum%i",i),";jet p_{T};counts",100,0,100);
          hphoeffden[isys][i] = new TH1D(Form("hphoeffden%i",i),";jet p_{T};counts",100,0,100);
          hjeteffnum[isys][i] = new TH1D(Form("hjeteffnum%i",i),";jet p_{T};counts",100,0,100);
          hjeteffden[isys][i] = new TH1D(Form("hjeteffden%i",i),";jet p_{T};counts",100,0,100);

          hphopurnum[isys][i] = new TH1D(Form("hphopurnum%i",i),";jet p_{T};counts",100,0,100);
          hphopurden[isys][i] = new TH1D(Form("hphopurden%i",i),";jet p_{T};counts",100,0,100);
          hjetpurnum[isys][i] = new TH1D(Form("hjetpurnum%i",i),";jet p_{T};counts",100,0,100);
          hjetpurden[isys][i] = new TH1D(Form("hjetpurden%i",i),";jet p_{T};counts",100,0,100);

          hpaireffnum[isys][i] = new TH1D(Form("hpaireffnum%i",i),";jet p_{T};counts",100,0,100);
          hpaireffden[isys][i] = new TH1D(Form("hpaireffden%i",i),";jet p_{T};counts",100,0,100);
          hpairpurnum[isys][i] = new TH1D(Form("hpairpurnum%i",i),";jet p_{T};counts",100,0,100);
          hpairpurden[isys][i] = new TH1D(Form("hpairpurden%i",i),";jet p_{T};counts",100,0,100);

          hphomissfake[isys][i] = new TH2D(Form("hphomissfake%i",i),";reco p_{T};truth p_{T}",101,0,101,101,0,101);
          hjetmissfake[isys][i] = new TH2D(Form("hjetmissfake%i",i),";reco p_{T};truth p_{T}",101,0,101,101,0,101);
          hpairmissfake[isys][i] = new TH2D(Form("hpairmissfake%i",i),";reco p_{T};truth p_{T}",nbins+1,0,nbins+1,nbins+1,0,nbins+1);

          hrecojetpt[isys][i] = new TH1D(Form("hrecojetpt%i",i),";jet p_{T,max};counts",100,0,100);
          htruthjetpt[isys][i] = new TH1D(Form("htruthjetpt%i",i),";jet p_{T,max};counts",100,0,100);
          jet_response[isys][i] = new RooUnfoldResponse(hrecojetpt[isys][i], htruthjetpt[isys][i]);

          hrecophopt[isys][i] = new TH1D(Form("hrecophopt%i",i),";pho p_{T,max};counts",100,0,100);
          htruthphopt[isys][i] = new TH1D(Form("htruthphopt%i",i),";pho p_{T,max};counts",100,0,100);
          pho_response[isys][i] = new RooUnfoldResponse(hrecophopt[isys][i], htruthphopt[isys][i]);

          hrecojetpt_half[isys][i] = new TH1D(Form("hrecojetpt_half%i",i),";jet p_{T,max};counts",100,0,100);
          htruthjetpt_half[isys][i] = new TH1D(Form("htruthjetpt_half%i",i),";jet p_{T,max};counts",100,0,100);
          jet_response_half[isys][i] = new RooUnfoldResponse(hrecojetpt_half[isys][i], htruthjetpt_half[isys][i]);

          hrecophopt_half[isys][i] = new TH1D(Form("hrecophopt_half%i",i),";pho p_{T,max};counts",100,0,100);
          htruthphopt_half[isys][i] = new TH1D(Form("htruthphopt_half%i",i),";pho p_{T,max};counts",100,0,100);
          pho_response_half[isys][i] = new RooUnfoldResponse(hrecophopt_half[isys][i], htruthphopt_half[isys][i]);

          for (int j = 0; j < 4; j++) { // 4  for ABCD
            hrecoxj_abcd[isys][i][j] = new TH1D(Form("hrecoxj%i_%i",i,j),";reco cluster p_{T}; x_{J#gamma}"      ,nbins,0,nbins);
            htruthxj_abcd[isys][i][j] = new TH1D(Form("htruthxj%i_%i",i,j),";truth cluster p_{T}; x_{J#gamma}"      ,nbins,0,nbins);
            hclusterpt_abcd[isys][i][j] = new TH1D(Form("hclusterpt_abcd%i_%i",i,j),";p_{T}^{lead cluster};Counts",ana::nPtBins,ana::ptBins);
            // MC only: truth-matched (dR < 0.1) subset of hclusterpt_abcd, for puritymaker.C's leakage
            // fractions, which are defined for true signal.
            hclusterpt_abcd_truthmatched[isys][i][j] = new TH1D(Form("hclusterpt_abcd_truthmatched%i_%i",i,j),";p_{T}^{lead cluster};Counts",ana::nPtBins,ana::ptBins);
          }
          hrecoxj[isys][i] = new TH1D(Form("hrecoxj%i",i),";reco cluster p_{T}; x_{J#gamma}"      ,nbins,0,nbins);
          htruthxj[isys][i] = new TH1D(Form("htruthxj%i",i),";truth cluster p_{T}; x_{J#gamma}"   ,nbins,0,nbins);
          jet_response2D[isys][i] = new RooUnfoldResponse(hrecoxj[isys][i], htruthxj[isys][i],Form("response_full_jetR%d",i),Form("response_%d",i));

          hrecoxj_half[isys][i] = new TH1D(Form("hrecoxj_half%i",i),";reco cluster p_{T}; x_{J#gamma}"      ,nbins,0,nbins);
          htruthxj_half[isys][i] = new TH1D(Form("htruthxj_half%i",i),";truth cluster p_{T}; x_{J#gamma}"   ,nbins,0,nbins);
          jet_response_half2D[isys][i] = new RooUnfoldResponse(hrecoxj_half[isys][i], htruthxj_half[isys][i],Form("response_half_jetR%d",i),Form("response_half_%d",i));
        }

        // In-situ tree: one per systag; the ir branch marks the radius (an event fills once per radius it
        // pairs at).
        const char * insitu_filename = isMC ?
            Form("%s/insitu/inputs/%s_%s_%s_insitu.root", ana::dir(), trigger.c_str(), sim.c_str(), systags[isys].c_str()) :
            Form("%s/insitu/inputs/%s_%s_insitu.root", ana::dir(), trigger.c_str(), systags[isys].c_str());
        insitu_file[isys] = TFile::Open(insitu_filename, "RECREATE");
        insitu_tree[isys] = new TTree("insitutree", "photon-jet pairs (all jet radii) for in-situ test");
        insitu_tree[isys]->Branch("pho_pt", &insitu_pho_pt[isys]);
        insitu_tree[isys]->Branch("jet_pt", &insitu_jet_pt[isys]);
        insitu_tree[isys]->Branch("abcd", &insitu_abcd[isys]);
        insitu_tree[isys]->Branch("weight", &insitu_weight[isys]);
        insitu_tree[isys]->Branch("ir", &insitu_ir[isys]);
        // thirdjet_pt: uncorrected third-jet pT where the threejet veto is left to the reader (Data,
        // threejet systag: the veto depends on p_a); -1 elsewhere.
        insitu_tree[isys]->Branch("thirdjet_pt", &insitu_third_pt[isys]);

        hpurity_num[isys]    = new TH2D("hpurity_num", ";p_{T};x_{j}", ana::nPtBinsUsed, ana::ptBinsUsed, ana::nUnfoldXjBins, ana::unfoldXjBins);
        hpurity_den[isys]    = new TH2D("hpurity_den", ";p_{T};x_{j}", ana::nPtBinsUsed, ana::ptBinsUsed, ana::nUnfoldXjBins, ana::unfoldXjBins);
        hpurity_num_1D[isys] = new TH1D("hpurity_num_1D", ";p_{T};counts", ana::nPtBinsUsed, ana::ptBinsUsed);
        hpurity_den_1D[isys] = new TH1D("hpurity_den_1D", ";p_{T};counts", ana::nPtBinsUsed, ana::ptBinsUsed);
        hpurity_num[isys]->Sumw2();
        hpurity_den[isys]->Sumw2();
        hpurity_num_1D[isys]->Sumw2();
        hpurity_den_1D[isys]->Sumw2();

        hphoIDeff_bdt[isys] = new TH2D("hphoIDeff_bdt", ";truth p_{T}^{#gamma};BDT score", ana::nPtBins, ana::ptBins, 100, 0.0, 1.0);
        hphoIDeff_iso[isys] = new TH2D("hphoIDeff_iso", ";truth p_{T}^{#gamma};isolation E_{T} [GeV]", ana::nPtBins, ana::ptBins, 120, -2.0, 10.0);
      }
    }
    // Single-systag constructor (macros/unfold.C).
    unfolder(string trigger, string sim, string systag) : unfolder(trigger, sim, vector<string>{systag}) {}

    ~unfolder();
    void fill_matrix();
    void unfold();
    void end();
    void savehists(TH1D * h[], int n);
    void savehists(TH2D * h[], int n);
    void savehists(RooUnfoldResponse * h[], int n);
    void savehists(TEfficiency * h[], int n);
    bool check_pair(jet_object jet, int ir, pho_object pho, bool isreco, float testPt = -1, float floorScale = 1.0);
    bool check_match(pho_object p1, pho_object p2, jet_object j1, jet_object j2);
    bool check_match(pho_object p1, pho_object p2);
    bool check_match(jet_object j1, jet_object j2);
    void set_dodraw(bool draw) {dodraw = draw; }
    // Progress line every 1000 entries; off for batch runs.
    void set_progress(bool show) {showProgress = show; }
    bool showProgress = true;

    template <typename T>
      float findmaxpt(const vector<T>& objs)
      {
        float maxpt = -1;
        for (const auto& obj : objs) {
          if (obj.pt > maxpt)
            maxpt = obj.pt;
        }
        return maxpt;
      }
  private:
    vector<string> systags;
    float dodraw = false;
    int count_isc = 0;
    int count_isj[ana::nJetR] = { 0 };
    int nentries = 0;
    TRandom rand;
    // Data/MC vz and cluster-pT reweighting (reweight/make_vz_pt_reweight.C), applied to all MC.
    Reweighter rw;

    // Members below are indexed [isys][ir] (some [isys][ir][abcd] or [isys]).
    vector<vector<RooUnfoldResponse*>> pho_response;
    vector<vector<RooUnfoldResponse*>> jet_response;
    vector<vector<RooUnfoldResponse*>> pho_response_half;
    vector<vector<RooUnfoldResponse*>> jet_response_half;

    vector<vector<TH1D*>> hphodr;
    vector<vector<TH1D*>> hjetdr;

    vector<vector<TH1D*>> hphoeffnum;
    vector<vector<TH1D*>> hphoeffden;
    vector<vector<TEfficiency*>> hphoeff;
    vector<vector<TH1D*>> hjeteffnum;
    vector<vector<TH1D*>> hjeteffden;
    vector<vector<TEfficiency*>> hjeteff;

    vector<vector<TH1D*>> hphopurnum;
    vector<vector<TH1D*>> hphopurden;
    vector<vector<TEfficiency*>> hphopur;
    vector<vector<TH1D*>> hjetpurnum;
    vector<vector<TH1D*>> hjetpurden;
    vector<vector<TEfficiency*>> hjetpur;

    vector<vector<TH1D*>> hpaireffnum;
    vector<vector<TH1D*>> hpaireffden;
    vector<vector<TEfficiency*>> hpaireff;
    vector<vector<TH1D*>> hpairpurnum;
    vector<vector<TH1D*>> hpairpurden;
    vector<vector<TEfficiency*>> hpairpur;

    vector<vector<TH2D*>> hphomissfake;
    vector<vector<TH2D*>> hjetmissfake;
    vector<vector<TH2D*>> hpairmissfake;

    vector<vector<TH1D*>> hrecojetpt;
    vector<vector<TH1D*>> htruthjetpt;
    vector<vector<TH2D*>> hjetresponse;
    vector<vector<TH1D*>> hrecophopt;
    vector<vector<TH1D*>> htruthphopt;
    vector<vector<TH2D*>> hphoresponse;

    vector<vector<TH1D*>> hrecojetpt_half;
    vector<vector<TH1D*>> htruthjetpt_half;
    vector<vector<TH2D*>> hjetresponse_half;
    vector<vector<TH1D*>> hrecophopt_half;
    vector<vector<TH1D*>> htruthphopt_half;
    vector<vector<TH2D*>> hphoresponse_half;

    int nbins = (ana::nPtBins) * (ana::nUnfoldXjBins+2); // +2 bins per pT bin for overflow and underflow
    vector<vector<RooUnfoldResponse*>> jet_response2D;
    vector<vector<RooUnfoldResponse*>> jet_response_half2D;

    vector<vector<vector<TH1D*>>> hrecoxj_abcd;    // [isys][ir][4] for ABCD
    vector<vector<vector<TH1D*>>> htruthxj_abcd;   // [isys][ir][4]
    vector<vector<vector<TH1D*>>> hclusterpt_abcd; // [isys][ir][abcd] reco cluster pT, for puritymaker.C
    vector<vector<vector<TH1D*>>> hclusterpt_abcd_truthmatched; // [isys][ir][abcd] truth-matched subset, for the leakage fractions
    vector<vector<TH1D*>> hrecoxj;
    vector<vector<TH1D*>> htruthxj;
    vector<vector<TH2D*>> hxjresponse;

    vector<vector<TH1D*>> hrecoxj_half;
    vector<vector<TH1D*>> htruthxj_half;
    vector<vector<TH2D*>> hxjresponse_half;

    vector<TH2D*> hpurity_num;
    vector<TH2D*> hpurity_den;
    vector<TH1D*> hpurity_num_1D;
    vector<TH1D*> hpurity_den_1D;

    // Truth-matched photons' raw BDT/isolation score vs truth pT, no ID cut: gives the efficiency of
    // any threshold vs pT without rerunning.
    vector<TH2D*> hphoIDeff_bdt;
    vector<TH2D*> hphoIDeff_iso;

    vector<TFile*> insitu_file;
    vector<TTree*> insitu_tree;
    vector<Float_t> insitu_pho_pt;
    vector<Float_t> insitu_jet_pt;
    vector<Int_t>   insitu_abcd;
    vector<Float_t> insitu_weight;
    vector<Int_t>   insitu_ir;
    vector<Float_t> insitu_third_pt;

};

#endif // UNFOLDER_H
