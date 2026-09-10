#ifndef ana_h
#define ana_h

#include <TROOT.h>
#include "TClonesArray.h"
#include <TLorentzVector.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <map>
#include <TEfficiency.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TFitResult.h>
#include <TGraphAsymmErrors.h>
#include <string>
using namespace std;

class ana {
  public :
    ana();
    ~ana() = default;

    static Bool_t   PassEtaCut(float eta, float vz);
    static Double_t GetShiftedEta(float _vz, float _eta);
    // ir (default 2 = R=0.4) selects which jet radius's purity to read - purity is
    // "of paired photons" (see puritymaker.C's hclusterpt_abcd%i_%i input, gated on
    // ispaired[ir]/pairing status that genuinely differs by jet radius), not a pure
    // photon-ID quantity independent of the jet, so it needs its own value per radius,
    // not one number reused everywhere. Defaulting to ir=2 keeps every caller that
    // doesn't pass ir reading exactly the pre-existing hists/purity_<systag>.root (see
    // purityFilename below) with zero behavior change.
    static float    getPurity(float low, float high, string systag = "nominal", int ir = 2);
    static float    getPurity(float val, string systag = "nominal", int ir = 2);
    static float    getPurityErrorLow(float low, float high, string systag = "nominal", int ir = 2);
    static float    getPurityErrorHigh(float low, float high, string systag = "nominal", int ir = 2);
    // Purity of the ABCD sideband region C (not region A) - same puritymaker.C bootstrap
    // (S/A quadratic solve), read back from the "combined_C" graph. n_s^C = c*S falls out
    // of the same leakage-corrected fit used for region A's purity, so this needs no
    // separate MC template or independent solve - see puritymaker.C's combine_hists().
    static float    getPurityC(float low, float high, string systag = "nominal", int ir = 2);
    static float    getPurityCErrorLow(float low, float high, string systag = "nominal", int ir = 2);
    static float    getPurityCErrorHigh(float low, float high, string systag = "nominal", int ir = 2);
    // hists/purity_<systag>.root - one file per systag holding every jet radius's
    // purity curve in its own ana::rnames[ir] subdirectory (see puritymaker.C, the only
    // writer of this file). getPurity/getPurityC etc. below select the radius back out
    // via that subdirectory.
    static string   purityFilename(const string & systag);
    static Int_t findPtBin(double value);
    static Int_t findUnfoldXjBin(double value);
    // Single overload only - a second (double,double,float,int) and third
    // (double,int,int) overload used to exist here but were unused dead code, and their
    // presence made ROOT/Cling's overload resolution for this call non-deterministic
    // across process runs (the interpreter would occasionally bind the float bdt argument
    // to the int showershape parameter instead of the intended double bdt parameter,
    // truncating it to 0 and silently zeroing the ABCD signal-region count for the whole
    // run). Do not re-add another findabcdBin overload without renaming it.
    static Int_t findabcdBin(double iso, double bdt, int bin = 0);
    static Int_t findxjBin(double value);
    static Int_t findBdtBin(double value);
    static Int_t findUnfoldBin(double xj, double pt);
    static Int_t findHadronBin(double value);
    static Int_t findEmfracBin(double value);

    static constexpr float sPHENIX_posx = 0.6;
    static constexpr float sPHENIX_posy = 0.85;
    static constexpr float posy_diff = 0.05;
    static constexpr size_t MaxClusters = 10000;
    static constexpr size_t MaxJets = 10000;

    static constexpr int nDims=7;
    static constexpr int NHIST = nDims*nDims;
    static constexpr int gridSize = 7;
    static constexpr double vzcut = 60; 
    static constexpr double oppnum = 7;
    static constexpr double oppden = 8;
    static constexpr double oppcut = oppnum*M_PI/oppden;
    static constexpr double tcut = 4;
    static constexpr double tlowcut = 0;
    static constexpr double thighcut = 4;
    static constexpr double radius = 93; 
    static constexpr double etacut = 1.1;
    static constexpr double etamin = -etacut;
    static constexpr double etamax = etacut;
    static constexpr double cluster_pt_cut = 10;

    static constexpr int nCalibBins = 7; // Uncalib, JES calibrated, JES calibrated and JER corrected, JES+JER+reweighted, JES+JER+photon smeared, JERhigh, JERlow
    static constexpr int nPtBins = 5;
    // Bin 0 (13-15 GeV) and the last bin (35-100 GeV) are migration-only buffers, not
    // reported physics results: bin 0 lets true/reco migration across the analysis's low-pT
    // edge be modeled instead of truncated (see the refolding-closure-test discussion that
    // motivated it - low-pT trigger/reco/purity performance was checked to hold down to
    // 13 GeV before adding this), matching the high-pT overflow bin's existing role. The
    // nPtBinsUsed reported bins are ptBins[firstUsedPtBin..firstUsedPtBin+nPtBinsUsed].
    static constexpr int nPtBinsUsed = 3;
    static constexpr int firstUsedPtBin = 1;
    static constexpr int nUnfoldXjBins = 16;
    // findabcdBin's `bin` selector: index 0 is the nominal ABCD grid; indices 1-5 each
    // override exactly one of the four cut boundaries below (isoBins/isoBinsHigh/
    // bdtGoodLow/bdtBadLow) to a Data/MC systag reprocessing - see systagAbcdBinArr in
    // unfolder.cc's fill_matrix() for the systag-name -> bin mapping, and CLAUDE.md /
    // drawing/draw_systematics.C's purityMembers for how the resulting five systags are
    // combined into one "Purity" systematic.
    static constexpr int nIsoBdtBins = 6;
    static constexpr int nBdtBins = 4;
    static constexpr int nxjBins = 3;
    static constexpr int nJetR = 7;
    static constexpr int n3jetBins = 2;
    static constexpr int nHadronBins = 3;
    static constexpr int nEmfracBins = 3;
    static constexpr int singleptlow = 15;
    static constexpr int singlepthigh = 35;

    // Definitive systag reprocessing list: nominal + every named systematic-variation
    // reprocessing of the same input tree (see unfolder.h's constructor comment for what
    // each one changes). This is the SINGLE place to add/remove a systag - every
    // consumer below reads it, so a new entry here automatically propagates:
    //   - macros/unfold_allsys.C's default systags list (unfolder.cc then fills that
    //     systag's histograms AND its insitu_tree for every jet radius - see
    //     unfolder.h's per-(systag,radius) insitu_tree construction)
    //   - insitu/run_grid.sh's systag loop (the in-situ JES scan then runs for it)
    //   - drawing/draw_systematics.C's `systematics` source list and asymmetric/symmetric
    //     classification (see asymmetricSystagPairs below)
    // A systag that isn't a "nominal +/- something" reprocessing at all (a different
    // generator sample like herwig, or a different unfolding setting like niterLow/
    // niterHigh/priorSensitivity) isn't a systags entry - those have no insitu_tree
    // equivalent and are wired up independently (draw_nonclosure.C,
    // draw_iteration_halfclosure.C, draw_prior_sensitivity.C).
    static const vector<string> systags;

    // Two-point (high/low) systematic pairs within `systags` - both names in each pair
    // must appear in `systags`. Per this project's asymmetric-systematic ground rule
    // (see CLAUDE.md / drawing/draw_systematics.C's sign-split combination), each pair
    // contributes its own signed per-bin value to whichever total (up/down) matches its
    // sign; every `systags` entry NOT listed here (e.g. threejet, narrowBDT, narrowISO)
    // is a symmetric single-sided source instead, contributing its full magnitude to
    // both totals. Adding a new "*_high"/"*_low" pair to `systags` also needs an entry
    // here to be treated as asymmetric - it is NOT auto-detected from the name, so a
    // typo'd or renamed pair fails loud (missing from every total) rather than silently
    // misclassifying.
    static const vector<pair<string,string>> asymmetricSystagPairs;
    static constexpr double ptBins[nPtBins+1] = {13,15,20,25,35,100};
    // Plain duplicate of ptBins[firstUsedPtBin..firstUsedPtBin+nPtBinsUsed] (15,20,25,35) -
    // kept as its own literal array (not pointer arithmetic into ptBins) since a few
    // variable-binning TH1/TH2 constructors need just the reported bins' edges, and ROOT's
    // Cling JIT (macros that #include unfolder.h directly, rather than only linking the
    // compiled .so, parse/compile it there) has had rough edges with pointer offsets into
    // static constexpr array members in default member initializers.
    static constexpr double ptBinsUsed[nPtBinsUsed+1] = {15,20,25,35};
    static constexpr double unfoldXjBins[nUnfoldXjBins+1] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.5,1.7,2.0};
    // ABCD grid boundaries, one column per nIsoBdtBins bin index (0=nominal, 1=narrowBDT,
    // 2=narrowISO, 3=narrowBDTbkg, 4=narrowISObkg, 5=wideISObkg - see nIsoBdtBins comment
    // above). Each non-nominal bin overrides exactly ONE of the four boundaries below and
    // leaves the other three at nominal, following PPG12's ABCD sideband-systematic
    // treatment (sPHENIX isolated-photon analysis note, Sec. 5.3): narrowBDT/narrowISO
    // shift the SIGNAL-side cut (bdtGoodLow/isoBins) that separates region A from B/C;
    // narrowBDTbkg/narrowISObkg/wideISObkg shift the BACKGROUND-side cut (bdtBadLow/
    // isoBinsHigh) that separates the B/C/D sideband definition from the excluded middle
    // ground - PPG12's analog of this boundary is a two-sided variation (tighter/looser
    // gap), unlike the one-sided BDT edges, so it gets two bins (tight=4, loose=5); this
    // project treats each as its own independent symmetrized single-sided source rather
    // than sign-splitting them against each other as an asymmetric pair (see
    // ana::asymmetricSystagPairs's comment). See drawing/draw_systematics.C's
    // purityMembers for how all five are then combined into one "Purity" systematic.
    static constexpr double isoBins[nIsoBdtBins]      = {2,   2,   1.5, 2,   2,   2  };
    static constexpr double isoBinsHigh[nIsoBdtBins]  = {4,   4,   4,   4,   3,   5  };
    static constexpr double bdtGoodHigh[nIsoBdtBins]  = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    static constexpr double bdtGoodLow[nIsoBdtBins]   = {0.8, 0.7, 0.8, 0.8, 0.8, 0.8}; 
    static constexpr double bdtBadHigh[nIsoBdtBins]   = {0.6, 0.6, 0.6, 0.6, 0.6, 0.6};
    static constexpr double bdtBadLow[nIsoBdtBins]    = {0.2, 0.2, 0.2, 0.1, 0.2, 0.2};
    static constexpr double bdtBins[nBdtBins+1] = {0.4,0.7,0.8,0.9,1.0};
    static constexpr double xjBins[nxjBins+1] = {0,0.3,0.7,2.0};
    static constexpr double hadronBins[nHadronBins][2] = {{20,25},{35,45},{50,60}};
    static constexpr double emfracBins[nEmfracBins+1] = {0,0.5,0.8,1};
    static constexpr double JetRs[nJetR] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    static constexpr const char * rnames[nJetR] = {"R02", "R03", "R04", "R05", "R06", "R07", "R08"};
    static constexpr double drcut[nJetR] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    static constexpr double jet_pt_cut[nJetR] = {3,3,3,3,3,3,3};
    //static constexpr double jet_calib_pt_cut[nJetR] = {3,3,3,3,3,3,3};
    static constexpr double jet_calib_pt_cut[nJetR] = {5,5,5,5,5,5,5};

    // In-situ JES calibration, one value per jet radius (insitu/grid_insitu.C's
    // purity-corrected best-fit p_a, full 7-radius scan - replaces the old single
    // R=0.4-derived 0.9446 that used to be reused across every radius in
    // unfolder.cc's fill_matrix()). jesTotalErrLow/High are that same scan's own
    // total (stat+syst combined in quadrature) asymmetric uncertainty on p_a, used
    // directly as the jes_high/jes_low systematic shift - replacing the old flat,
    // unsourced +-0.03 placeholder. jes_high subtracts jesTotalErrLow (a smaller p_a
    // means a bigger 1/p_a correction, i.e. Data's jet pT scaled UP more relative to
    // nominal); jes_low adds jesTotalErrHigh (a larger p_a, smaller correction) -
    // the same sign convention the old symmetric +-0.03 followed.
    // These three arrays are generated, not hand-edited: insitu/run_grid.sh's full
    // systag sweep ends by running insitu/draw_jes_summary.C, which rewrites these exact
    // literals in place from its own freshly computed per-radius pa/totalLow/totalHigh -
    // see that macro's updateAnaHeader()/completeness guard. A manual edit here survives
    // only until the next full run_grid.sh sweep overwrites it; rebuild via src/make.sh
    // afterward for the new numbers to reach the compiled library.
    static constexpr double jesNominal[nJetR]      = {0.9162, 0.9266, 0.9358, 0.9310, 0.9305, 0.9279, 0.9482};
    static constexpr double jesTotalErrLow[nJetR]  = {0.0320, 0.0307, 0.0237, 0.0164, 0.0121, 0.0120, 0.0205};
    static constexpr double jesTotalErrHigh[nJetR] = {0.0261, 0.0324, 0.0300, 0.0218, 0.0180, 0.0194, 0.0253};

  private:
};

#endif
