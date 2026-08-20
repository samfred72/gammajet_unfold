#ifndef UNFOLDER_H
#define UNFOLDER_H

#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/object.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/pho_object.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/jet_object.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/treeuser.h"
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
    // systag selects a systematic-variation reprocessing of the SAME input tree (as
    // opposed to `sim`, which selects a genuinely different input tree/sample - pythia
    // vs herwig). Recognized values: "nominal" (default), "JERhigh", "JERlow",
    // "emscale_high", "emscale_low" (all four MC-only - Data has no JER-smearing or
    // EM-scale-shifted variant of itself, since smearing is applied to MC to match Data's
    // resolution, not something Data itself has; silently falls back to nominal
    // JES-corrected reco jet pT/cluster pT for Data regardless of systag), "jes_high",
    // "jes_low" (the reverse - Data-only, varying the +-0.03 uncertainty on Data's 0.977
    // in-situ JES correction; MC is already on-scale and falls back to nominal regardless
    // of systag), "threejet", "narrowBDT", "narrowISO". See fill_matrix() for exactly what
    // each one changes.
    unfolder(string trigger, string sim = "pythia", string systag = "nominal") : treeuser(trigger, sim), systag(systag) {
      treesetup();
      for (int i = 0; i < ana::nJetR; i++) {
        hphodr[i] = new TH1D(Form("hphodr%i",i),";#DeltaR_{truth,reco};counts", 100,0,0.4);
        hjetdr[i] = new TH1D(Form("hjetdr%i",i),";#DeltaR_{truth,reco};counts",100,0,0.4);
    
        hphoeffnum[i] = new TH1D(Form("hphoeffnum%i",i),";jet p_{T};counts",100,0,100);
        hphoeffden[i] = new TH1D(Form("hphoeffden%i",i),";jet p_{T};counts",100,0,100);
        hjeteffnum[i] = new TH1D(Form("hjeteffnum%i",i),";jet p_{T};counts",100,0,100);
        hjeteffden[i] = new TH1D(Form("hjeteffden%i",i),";jet p_{T};counts",100,0,100);
        
        hphopurnum[i] = new TH1D(Form("hphopurnum%i",i),";jet p_{T};counts",100,0,100);
        hphopurden[i] = new TH1D(Form("hphopurden%i",i),";jet p_{T};counts",100,0,100);
        hjetpurnum[i] = new TH1D(Form("hjetpurnum%i",i),";jet p_{T};counts",100,0,100);
        hjetpurden[i] = new TH1D(Form("hjetpurden%i",i),";jet p_{T};counts",100,0,100);
        
        hpaireffnum[i] = new TH1D(Form("hpaireffnum%i",i),";jet p_{T};counts",100,0,100);
        hpaireffden[i] = new TH1D(Form("hpaireffden%i",i),";jet p_{T};counts",100,0,100);
        hpairpurnum[i] = new TH1D(Form("hpairpurnum%i",i),";jet p_{T};counts",100,0,100);
        hpairpurden[i] = new TH1D(Form("hpairpurden%i",i),";jet p_{T};counts",100,0,100);

        hphomissfake[i] = new TH2D(Form("hphomissfake%i",i),";reco p_{T};truth p_{T}",101,0,101,101,0,101);
        hjetmissfake[i] = new TH2D(Form("hjetmissfake%i",i),";reco p_{T};truth p_{T}",101,0,101,101,0,101);
        hpairmissfake[i] = new TH2D(Form("hpairmissfake%i",i),";reco p_{T};truth p_{T}",nbins+1,0,nbins+1,nbins+1,0,nbins+1);
        
        hrecojetpt[i] = new TH1D(Form("hrecojetpt%i",i),";jet p_{T,max};counts",100,0,100);
        htruthjetpt[i] = new TH1D(Form("htruthjetpt%i",i),";jet p_{T,max};counts",100,0,100);
        jet_response[i] = new RooUnfoldResponse(hrecojetpt[i], htruthjetpt[i]);
        
        hrecophopt[i] = new TH1D(Form("hrecophopt%i",i),";pho p_{T,max};counts",100,0,100);
        htruthphopt[i] = new TH1D(Form("htruthphopt%i",i),";pho p_{T,max};counts",100,0,100);
        pho_response[i] = new RooUnfoldResponse(hrecophopt[i], htruthphopt[i]);
        
        hrecojetpt_half[i] = new TH1D(Form("hrecojetpt_half%i",i),";jet p_{T,max};counts",100,0,100);
        htruthjetpt_half[i] = new TH1D(Form("htruthjetpt_half%i",i),";jet p_{T,max};counts",100,0,100);
        jet_response_half[i] = new RooUnfoldResponse(hrecojetpt_half[i], htruthjetpt_half[i]);
        
        hrecophopt_half[i] = new TH1D(Form("hrecophopt_half%i",i),";pho p_{T,max};counts",100,0,100);
        htruthphopt_half[i] = new TH1D(Form("htruthphopt_half%i",i),";pho p_{T,max};counts",100,0,100);
        pho_response_half[i] = new RooUnfoldResponse(hrecophopt_half[i], htruthphopt_half[i]);
       
        for (int j = 0; j < 4; j++) { // 4  for ABCD
          hrecoxj_abcd[i][j] = new TH1D(Form("hrecoxj%i_%i",i,j),";reco cluster p_{T}; x_{J#gamma}"      ,nbins,0,nbins);
          htruthxj_abcd[i][j] = new TH1D(Form("htruthxj%i_%i",i,j),";truth cluster p_{T}; x_{J#gamma}"      ,nbins,0,nbins);
          hclusterpt_abcd[i][j] = new TH1D(Form("hclusterpt_abcd%i_%i",i,j),";p_{T}^{lead cluster};Counts",ana::nPtBins,ana::ptBins);
        }
        hrecoxj[i] = new TH1D(Form("hrecoxj%i",i),";reco cluster p_{T}; x_{J#gamma}"      ,nbins,0,nbins);
        htruthxj[i] = new TH1D(Form("htruthxj%i",i),";truth cluster p_{T}; x_{J#gamma}"   ,nbins,0,nbins);
        hunfoldxj[i] = new TH1D(Form("hunfoldxj%i",i),";unfold cluster p_{T}; x_{J#gamma}",nbins,0,nbins);
        jet_response2D[i] = new RooUnfoldResponse(hrecoxj[i], htruthxj[i],Form("response_full_jetR%d",i),Form("response_%d",i));

        hrecoxj_half[i] = new TH1D(Form("hrecoxj_half%i",i),";reco cluster p_{T}; x_{J#gamma}"      ,nbins,0,nbins);
        htruthxj_half[i] = new TH1D(Form("htruthxj_half%i",i),";truth cluster p_{T}; x_{J#gamma}"   ,nbins,0,nbins);
        hunfoldxj_half[i] = new TH1D(Form("hunfoldxj_half%i",i),";unfold cluster p_{T}; x_{J#gamma}",nbins,0,nbins);
        jet_response_half2D[i] = new RooUnfoldResponse(hrecoxj_half[i], htruthxj_half[i],Form("response_half_jetR%d",i),Form("response_half_%d",i));

      }
      hpurity_num->   Sumw2();
      hpurity_den->   Sumw2();
      hpurity_num_1D->Sumw2();
      hpurity_den_1D->Sumw2();

      // In-situ test tree: photon pt, jet pt, and ABCD region for R=0.4 events with a
      // valid jet pair. Opened last so it doesn't inherit ownership of the histograms
      // created above (ROOT associates new TH1/TH2 objects with whatever file is
      // currently open, i.e. the raw input file opened by treeuser).
      const char * insitu_filename = isMC ?
          Form("/home/samson72/sphnx/gammajet_unfold/insitu/%s_%s_%s_insitu.root", trigger.c_str(), sim.c_str(), systag.c_str()) :
          Form("/home/samson72/sphnx/gammajet_unfold/insitu/%s_%s_insitu.root", trigger.c_str(), systag.c_str());
      insitu_file = TFile::Open(insitu_filename, "RECREATE");
      insitu_tree = new TTree("insitutree", "photon-jet pairs (R=0.4) for in-situ test");
      insitu_tree->Branch("pho_pt", &insitu_pho_pt);
      insitu_tree->Branch("jet_pt", &insitu_jet_pt);
      insitu_tree->Branch("abcd", &insitu_abcd);
    }
    ~unfolder(); 
    void fill_matrix();
    void unfold();
    void end();
    void savehists(TH1D * h[], int n);
    void savehists(TH2D * h[], int n);
    void savehists(RooUnfoldResponse * h[], int n);
    void savehists(TEfficiency * h[], int n);
    bool check_pair(jet_object jet, int ir, pho_object pho, bool isreco);
    bool check_match(pho_object p1, pho_object p2, jet_object j1, jet_object j2);
    bool check_match(pho_object p1, pho_object p2);
    bool check_match(jet_object j1, jet_object j2);
    void set_dodraw(bool draw) {dodraw = draw; }
    
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
    string systag;
    float dodraw = false;
    int count_isc = 0;
    int count_isj[ana::nJetR] = { 0 };
    int nentries = 0;
    int niterate = 1;
    TRandom rand;
    
    RooUnfoldResponse * pho_response[ana::nJetR];
    RooUnfoldResponse * jet_response[ana::nJetR];
    RooUnfoldResponse * pho_response_half[ana::nJetR];
    RooUnfoldResponse * jet_response_half[ana::nJetR];

    TH1D * hphodr[ana::nJetR];
    TH1D * hjetdr[ana::nJetR];

    TH1D * hphoeffnum[ana::nJetR];
    TH1D * hphoeffden[ana::nJetR];
    TEfficiency * hphoeff[ana::nJetR];
    TH1D * hjeteffnum[ana::nJetR];
    TH1D * hjeteffden[ana::nJetR];
    TEfficiency * hjeteff[ana::nJetR];
    
    TH1D * hphopurnum[ana::nJetR];
    TH1D * hphopurden[ana::nJetR];
    TEfficiency * hphopur[ana::nJetR];
    TH1D * hjetpurnum[ana::nJetR];
    TH1D * hjetpurden[ana::nJetR];
    TEfficiency * hjetpur[ana::nJetR];
  
    TH1D * hpaireffnum[ana::nJetR];
    TH1D * hpaireffden[ana::nJetR];
    TEfficiency * hpaireff[ana::nJetR];
    TH1D * hpairpurnum[ana::nJetR];
    TH1D * hpairpurden[ana::nJetR];
    TEfficiency * hpairpur[ana::nJetR];

    TH2D * hphomissfake[ana::nJetR];
    TH2D * hjetmissfake[ana::nJetR];
    TH2D * hpairmissfake[ana::nJetR];

    TH1D * hrecojetpt[ana::nJetR];
    TH1D * htruthjetpt[ana::nJetR];
    TH1D * hunfoldjetpt[ana::nJetR];
    TH2D * hjetresponse[ana::nJetR];
    TH1D * hrecophopt[ana::nJetR];
    TH1D * htruthphopt[ana::nJetR];
    TH1D * hunfoldphopt[ana::nJetR];
    TH2D * hphoresponse[ana::nJetR];
    
    TH1D * hrecojetpt_half[ana::nJetR];
    TH1D * htruthjetpt_half[ana::nJetR];
    TH1D * hunfoldjetpt_half[ana::nJetR];
    TH2D * hjetresponse_half[ana::nJetR];
    TH1D * hrecophopt_half[ana::nJetR];
    TH1D * htruthphopt_half[ana::nJetR];
    TH1D * hunfoldphopt_half[ana::nJetR];
    TH2D * hphoresponse_half[ana::nJetR];
    
   
    int nbins = (ana::nPtBins) * (ana::nUnfoldXjBins+2); // +2 bins per pT bin for overflow and underflow
    RooUnfoldResponse * jet_response2D[ana::nJetR];
    RooUnfoldResponse * jet_response_half2D[ana::nJetR];
    
    TH1D * hrecoxj_abcd    [ana::nJetR][4];
    TH1D * htruthxj_abcd   [ana::nJetR][4];
    TH1D * hclusterpt_abcd [ana::nJetR][4]; // reco cluster pT per ABCD region, binned in ana::ptBins - what puritymaker.C needs
    TH1D * hrecoxj         [ana::nJetR];
    TH1D * htruthxj        [ana::nJetR];
    TH1D * hunfoldxj       [ana::nJetR];
    TH2D * hxjresponse     [ana::nJetR];
    
    TH1D * hrecoxj_half    [ana::nJetR];
    TH1D * htruthxj_half   [ana::nJetR];
    TH1D * hunfoldxj_half  [ana::nJetR];
    TH2D * hxjresponse_half[ana::nJetR];

    // ana::ptBinsUsed (not ana::ptBins from index 0): the reported pT bins no longer start
    // at index 0 now that ana::ptBins has a low-pT migration-only buffer bin at index 0 -
    // see ana.h.
    TH2D * hpurity_num = new TH2D("hpurity_num", ";p_{T};x_{j}", ana::nPtBinsUsed, ana::ptBinsUsed, ana::nUnfoldXjBins, ana::unfoldXjBins);
    TH2D * hpurity_den = new TH2D("hpurity_den", ";p_{T};x_{j}", ana::nPtBinsUsed, ana::ptBinsUsed, ana::nUnfoldXjBins, ana::unfoldXjBins);
    TH1D * hpurity_num_1D = new TH1D("hpurity_num_1D", ";p_{T};counts", ana::nPtBinsUsed, ana::ptBinsUsed);
    TH1D * hpurity_den_1D = new TH1D("hpurity_den_1D", ";p_{T};counts", ana::nPtBinsUsed, ana::ptBinsUsed);

    // Truth-matched photons' raw BDT/isolation score vs. truth pT, with no ID cut applied
    // (unlike hphoeffnum/den, which only track reco+match efficiency, also with no ID cut).
    // Slicing this at any score threshold per pT-bin gives the true-photon efficiency of
    // that threshold as a function of pT, without rerunning the pipeline per threshold -
    // used to check whether the fixed nominal cut (ana::bdtGoodLow[0]/ana::isoBins[0])
    // maps to a pT-dependent efficiency.
    TH2D * hphoIDeff_bdt = new TH2D("hphoIDeff_bdt", ";truth p_{T}^{#gamma};BDT score", ana::nPtBins, ana::ptBins, 100, 0.0, 1.0);
    TH2D * hphoIDeff_iso = new TH2D("hphoIDeff_iso", ";truth p_{T}^{#gamma};isolation E_{T} [GeV]", ana::nPtBins, ana::ptBins, 120, -2.0, 10.0);

    TFile * insitu_file;
    TTree * insitu_tree;
    Float_t insitu_pho_pt;
    Float_t insitu_jet_pt;
    Int_t   insitu_abcd;

};

#endif // UNFOLDER_H
