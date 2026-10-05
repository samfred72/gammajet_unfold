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
    // systags: one or more systematic-variation reprocessings of the SAME input tree
    // (as opposed to `sim`, which selects a genuinely different input tree/sample -
    // pythia vs herwig), processed together in a single pass over the tree - fill_matrix()
    // reads each entry once and, per systag, recomputes that systag's reco-level
    // quantities and fills that systag's own histogram set (see fill_matrix()'s inner
    // systag loop). This replaces re-running the whole analysis once per systag (each
    // re-reading/re-deserializing the same tree from disk) with one read plus cheap
    // per-systag arithmetic - every systag switch below is a closed-form shift or a
    // selection among branches already present in a single tree entry (e.g. the JER
    // high/low variants are pre-computed branches, not fresh random draws), so nothing
    // requires a second pass. Recognized values: "nominal" (default), "JERhigh", "JERlow",
    // "emscale_high", "emscale_low", "EMRhigh", "EMRlow" (all six MC-only - Data has no
    // JER-smearing, EM-scale-shifted, or EM-resolution-smeared variant of itself, since
    // smearing is applied to MC to match Data's resolution, not something Data itself has;
    // silently falls back to nominal JES-corrected reco jet pT/cluster pT for Data
    // regardless of systag), "jes_high", "jes_low" (the reverse - Data-only, varying the
    // per-radius correction (ana::jesNominal[ir]) by the nominal in-situ fit's statistical
    // uncertainty (ana::jesStatErrLow/High[ir]) - every other systag instead corrects Data
    // with its own in-situ p_a (ana::jesForSystag), see src/ana.h and unfolder.cc; MC is
    // already on-scale and falls back to nominal regardless of systag),
    // "threejet", "narrowBDT",
    // "narrowISO", "narrowBDTbkg", "narrowISObkg", "wideISObkg" (the last five each shift
    // one ABCD sideband boundary - narrowBDT/narrowISO on the signal-side cut, the other
    // three on the background-side cut - applied identically to Data and MC, unlike the
    // MC-only/Data-only pairs above; see ana.h's isoBins/isoBinsHigh/bdtGoodLow/bdtBadLow
    // comment). See fill_matrix() for exactly what each one changes.
    //
    // TH1::AddDirectory(kFALSE) below is essential: every systag's histogram set uses the
    // SAME names (e.g. "hphodr0") as every other systag's, by design - so that end()'s
    // per-systag Write() calls reproduce exactly today's per-systag output file format
    // with no renaming. Multiple same-named TH1/TH2 objects living in memory at once is
    // only safe because AddDirectory(kFALSE) stops them from being auto-registered into
    // (and colliding within) whatever TFile/TDirectory happens to be "current" at
    // construction time - each is a free-standing object we explicitly Write() ourselves.
    unfolder(string trigger, string sim, vector<string> systags) : treeuser(trigger, sim), systags(systags) {
      // See treeuser::disableBranchesUnusedByUnfolder()'s comment - skips I/O/
      // decompression for branches this pipeline never reads.
      disableBranchesUnusedByUnfolder();
      gErrorIgnoreLevel = kWarning;
      TH1::AddDirectory(kFALSE);
      // Every histogram below gets filled with the per-event vz/cluster-pT MC weight
      // (fill_matrix()'s mcWeight, isMC only - Data always uses weight 1) - Sumw2 makes
      // GetBinError() reflect that (sqrt(sum w_i^2) instead of sqrt(N)) instead of
      // silently understating errors on reweighted MC. Applies to every TH1D/TH2D made in
      // this constructor, so no need to call Sumw2() on each one individually.
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
            // MC-only, truth-matched (photon deltaR<0.1) subset of hclusterpt_abcd above -
            // see unfolder.cc's fill for why: macros/puritymaker.C's leakage-fraction
            // templates (fp[i]=hp[i]/hp[0]) need the true-signal-only ABCD ratio the method
            // is actually defined with, not an all-reconstructed-cluster ratio. Data's own
            // hclusterpt_abcd (unmatched, all clusters) is untouched and still correct as
            // the genuine ABCD counts a data-driven method has to work with.
            hclusterpt_abcd_truthmatched[isys][i][j] = new TH1D(Form("hclusterpt_abcd_truthmatched%i_%i",i,j),";p_{T}^{lead cluster};Counts",ana::nPtBins,ana::ptBins);
          }
          hrecoxj[isys][i] = new TH1D(Form("hrecoxj%i",i),";reco cluster p_{T}; x_{J#gamma}"      ,nbins,0,nbins);
          htruthxj[isys][i] = new TH1D(Form("htruthxj%i",i),";truth cluster p_{T}; x_{J#gamma}"   ,nbins,0,nbins);
          jet_response2D[isys][i] = new RooUnfoldResponse(hrecoxj[isys][i], htruthxj[isys][i],Form("response_full_jetR%d",i),Form("response_%d",i));

          hrecoxj_half[isys][i] = new TH1D(Form("hrecoxj_half%i",i),";reco cluster p_{T}; x_{J#gamma}"      ,nbins,0,nbins);
          htruthxj_half[isys][i] = new TH1D(Form("htruthxj_half%i",i),";truth cluster p_{T}; x_{J#gamma}"   ,nbins,0,nbins);
          jet_response_half2D[isys][i] = new RooUnfoldResponse(hrecoxj_half[isys][i], htruthxj_half[isys][i],Form("response_half_jetR%d",i),Form("response_half_%d",i));
        }

        // In-situ test tree: photon pt, jet pt, ABCD region, and jet radius for events
        // with a valid jet pair at ANY ana::nJetR radius - one tree per systag (not per
        // radius; an "ir" branch distinguishes rows instead), same naming convention as
        // the original single-systag file (the systag component of the filename already
        // makes it unique on disk) - keeps the insitu/ file count at nsys instead of
        // nsys*ana::nJetR. A single event can Fill() up to ana::nJetR times (once per
        // radius it pairs at), so this tree is up to ana::nJetR times longer than the
        // old R=0.4-only tree, not wider.
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
    // Convenience overload matching the original single-systag call signature (still
    // used by macros/unfold.C) - delegates to the vector<string> constructor above with
    // a single-element list, so a one-systag run goes through exactly the same code path
    // as a multi-systag one (nsys=1 is the validation baseline for the refactor).
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
    // Per-event progress line (every 1000 entries). Off for batch/pipeline runs, where it
    // only bloats the logs (and has previously been mistaken for a slowdown).
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
    // Data/MC vz and cluster-pT reweighting (reweight/make_vz_pt_reweight.C) - applied to
    // every isMC event in fill_matrix(), regardless of trigger or sim ("pythia"/"herwig"):
    // it corrects a Data-vs-detector-simulation vertex/threshold mismatch, not something
    // specific to one generator or trigger sample.
    Reweighter rw;

    // Every member below is indexed [isys][ir] (or [isys][ir][abcdRegion], or just
    // [isys] for the handful that aren't per-jet-radius) - one full copy per systag in
    // `systags`, since fill_matrix() computes each systag's reco-level quantities
    // independently per event (see that function's inner systag loop). This is the
    // direct multi-systag generalization of what used to be single, unindexed members
    // when systag was fixed at construction.
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
    vector<vector<vector<TH1D*>>> hclusterpt_abcd; // [isys][ir][4] - reco cluster pT per ABCD region, what puritymaker.C needs
    vector<vector<vector<TH1D*>>> hclusterpt_abcd_truthmatched; // [isys][ir][4] - MC-only, truth-matched subset of the above; what puritymaker.C's leakage fractions actually need
    vector<vector<TH1D*>> hrecoxj;
    vector<vector<TH1D*>> htruthxj;
    vector<vector<TH2D*>> hxjresponse;

    vector<vector<TH1D*>> hrecoxj_half;
    vector<vector<TH1D*>> htruthxj_half;
    vector<vector<TH2D*>> hxjresponse_half;

    // ana::ptBinsUsed (not ana::ptBins from index 0): the reported pT bins no longer start
    // at index 0 now that ana::ptBins has a low-pT migration-only buffer bin at index 0 -
    // see ana.h.
    vector<TH2D*> hpurity_num;
    vector<TH2D*> hpurity_den;
    vector<TH1D*> hpurity_num_1D;
    vector<TH1D*> hpurity_den_1D;

    // Truth-matched photons' raw BDT/isolation score vs. truth pT, with no ID cut applied
    // (unlike hphoeffnum/den, which only track reco+match efficiency, also with no ID cut).
    // Slicing this at any score threshold per pT-bin gives the true-photon efficiency of
    // that threshold as a function of pT, without rerunning the pipeline per threshold -
    // used to check whether the fixed nominal cut (ana::bdtGoodLow[0]/ana::isoBins[0])
    // maps to a pT-dependent efficiency.
    vector<TH2D*> hphoIDeff_bdt;
    vector<TH2D*> hphoIDeff_iso;

    vector<TFile*> insitu_file;
    vector<TTree*> insitu_tree;
    vector<Float_t> insitu_pho_pt;
    vector<Float_t> insitu_jet_pt;
    vector<Int_t>   insitu_abcd;
    vector<Float_t> insitu_weight;
    vector<Int_t>   insitu_ir;

};

#endif // UNFOLDER_H
