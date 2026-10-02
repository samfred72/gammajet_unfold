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

    // EM-calorimeter scale/resolution systematic inputs - values and prescription taken from
    // the PPG12 isolated-photon analysis note (sPH-ppg-2024-012, Sec. 5.1-5.2), so this
    // analysis's EM systematics are the same as PPG12's. (PPG12's additional residual
    // energy-scale non-linearity term is deliberately not included - this analysis keeps
    // its own emscale_high/low + EMRhigh/low systag structure, see ana::systags.)
    //
    // emscaleShift: +-1.48% cluster-pT scale shift (Calorimeter Calibration Working Group),
    // applied by unfolder.cc as emscale_high/emscale_low.
    static constexpr float emscaleShift = 0.0148;
    // emResolutionSigma: fractional width of the extra Gaussian smearing applied to MC
    // cluster pT so MC's resolution matches Data's. sigma_extra(E) = sqrt(max(0,
    // sigma_data(E)^2 - sigma_MC(E)^2)), with sigma(E)/E = sqrt(p0^2/E + p1^2/E^2 + p2^2),
    // E in GeV, evaluated at the TRUTH cluster pT (PPG12 convention - the smearing is
    // then drawn as N(0, sigma_extra*E_truth) and added to the reco pT). Zero for E below
    // ~12.8 GeV (MC already wider than Data), ~1% at 15 GeV, 2% at 20 GeV, 2.4% at 36 GeV
    // for the nominal parameters. emrVariant selects the Data-resolution parameter set:
    // emrNominal (0.15,0.05,0.05), emrHigh = wider Data (0.13,0.08,0.08, ~6% roughly flat),
    // emrLow = no extra smearing at all (sigma_extra = 0). MC parameters are always
    // (0.185,0,0.040), fit to unsmeared pythia8 tight+iso prompt-photon response.
    static constexpr int emrNominal = 0;
    static constexpr int emrHigh    = 1;
    static constexpr int emrLow     = -1;
    static float    emResolutionSigma(float truthPt, int emrVariant = emrNominal);

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
    // Calorimeter acceptance |eta| < 1.1: jets are required within |eta_jet| < etacut - R.
    static constexpr double etacut = 1.1;
    // Photon acceptance |eta^gamma| < photonEtaCut, applied to reco AND truth photons
    // (unfolder::check_pair serves both, so it also defines the particle-level fiducial
    // region). 0.7, as in PPG12 (sPH-CONF-JET-2025-02, |eta^gamma| < 0.7): at |eta| > 0.7
    // the R = 0.4 isolation cone leaves the EMCal acceptance and the purity becomes
    // eta-dependent (PPG18 review issue 10, decided Sep 28 2026 - see claude_checks/purity/).
    static constexpr double photonEtaCut = 0.7;
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
    // threejet systematic: veto the event if the third jet (treeuser::thirdjet_pt, in the
    // same pT definition as the recoil jet) is above this. Reco only, any eta and any dR
    // from the leading jet; no truth-level veto, so the unfolded observable is unchanged.
    static constexpr double thirdJetPtCut = 5;

    // In-situ JES calibration (insitu/grid_insitu.C's purity-corrected best-fit p_a), one
    // value per jet radius. Data jet pT is divided by p_a in unfolder.cc's fill_matrix().
    //
    // Each systematic variation is propagated ONCE, coherently (PPG18 review issue 5,
    // Sep 28 2026): every systag that has its own in-situ scan (JER, EM scale, EMR,
    // threejet, the ABCD-boundary variations) corrects Data with the p_a ITS OWN scan
    // extracted - jesBySystag below - so the variation's effect on the jet energy scale
    // is carried inside that variation's own systematic. jesStatErrLow/High are only the
    // nominal fit's statistical (Delta chi2 = 1) uncertainty; jes_high/jes_low shift the
    // nominal p_a by exactly that, so the JES systematic no longer re-counts the other
    // variations. (Previously jes_high/jes_low used a "total" uncertainty that already
    // contained every other systag's in-situ shift in quadrature, while those systags were
    // also unfolded with the NOMINAL p_a and entered the total as their own sources - so
    // each one was counted twice, once directly and once inside JES.)
    // jes_high subtracts jesStatErrLow (a smaller p_a means a bigger 1/p_a correction,
    // i.e. Data's jet pT scaled UP more relative to nominal); jes_low adds jesStatErrHigh.
    //
    // Everything from here to "END jesBySystag" is generated, not hand-edited:
    // insitu/run_grid.sh's full systag sweep ends by running insitu/draw_jes_summary.C,
    // which rewrites these literals from the scan outputs (see its updateAnaHeader() and
    // completeness guard). Rebuild via src/make.sh afterward for the numbers to reach the
    // compiled library. jesSystagNames must list exactly ana::systags, in order;
    // jesForSystag() below looks entries up by name and warns on a mismatch.
    static constexpr double jesNominal[nJetR]      = {0.9034, 0.9158, 0.9265, 0.9259, 0.9292, 0.9295, 0.9519};
    static constexpr double jesStatErrLow[nJetR]   = {0.0050, 0.0114, 0.0058, 0.0069, 0.0077, 0.0080, 0.0087};
    static constexpr double jesStatErrHigh[nJetR]  = {0.0105, 0.0070, 0.0091, 0.0086, 0.0073, 0.0080, 0.0086};
    // BEGIN jesBySystag
    static constexpr int nJesSystags = 15;
    static constexpr const char * jesSystagNames[nJesSystags] = {"nominal", "JERhigh", "JERlow", "emscale_high", "emscale_low", "jes_high", "jes_low", "threejet", "narrowBDT", "narrowISO", "EMRhigh", "EMRlow", "narrowBDTbkg", "narrowISObkg", "wideISObkg"};
    static constexpr double jesBySystag[nJesSystags][nJetR] = {
      {0.9034, 0.9158, 0.9265, 0.9259, 0.9292, 0.9295, 0.9519}, // nominal
      {0.8860, 0.9033, 0.9217, 0.9194, 0.9266, 0.9265, 0.9489}, // JERhigh
      {0.9145, 0.9244, 0.9322, 0.9308, 0.9336, 0.9327, 0.9539}, // JERlow
      {0.9091, 0.9180, 0.9321, 0.9310, 0.9343, 0.9344, 0.9566}, // emscale_high
      {0.9034, 0.9038, 0.9226, 0.9215, 0.9243, 0.9245, 0.9478}, // emscale_low
      {0.9034, 0.9158, 0.9265, 0.9258, 0.9294, 0.9291, 0.9514}, // jes_high
      {0.9034, 0.9158, 0.9265, 0.9266, 0.9295, 0.9298, 0.9519}, // jes_low
      {0.8000, 0.8686, 0.8964, 0.9157, 0.9213, 0.9231, 0.9568}, // threejet
      {0.9034, 0.9180, 0.9322, 0.9297, 0.9272, 0.9289, 0.9506}, // narrowBDT
      {0.8946, 0.9033, 0.9200, 0.9182, 0.9218, 0.9256, 0.9424}, // narrowISO
      {0.9145, 0.9262, 0.9419, 0.9396, 0.9440, 0.9430, 0.9649}, // EMRhigh
      {0.9034, 0.9158, 0.9265, 0.9258, 0.9272, 0.9287, 0.9514}, // EMRlow
      {0.9091, 0.9180, 0.9303, 0.9310, 0.9348, 0.9361, 0.9536}, // narrowBDTbkg
      {0.9034, 0.9098, 0.9265, 0.9252, 0.9272, 0.9280, 0.9510}, // narrowISObkg
      {0.9077, 0.9158, 0.9265, 0.9266, 0.9295, 0.9297, 0.9512}, // wideISObkg
    };
    // END jesBySystag
    // p_a for a given systag and radius from jesBySystag; falls back to jesNominal[ir]
    // (with a warning) if the systag has no entry.
    static double   jesForSystag(const string & systag, int ir);

  private:
};

#endif
