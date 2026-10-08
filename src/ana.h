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
    // Repository root ($GAMMAJET_UNFOLD); exits if unset.
    static const char * dir();
    // dir() + "/" + rel, valid for the whole program (usable wherever a string literal was).
    static const char * path(const std::string & rel);

    ana();
    ~ana() = default;

    static Bool_t   PassEtaCut(float eta, float vz);
    static Double_t GetShiftedEta(float _vz, float _eta);
    // Purity of paired photons, per jet radius ir (pairing differs by radius). Default R = 0.4.
    static float    getPurity(float low, float high, string systag = "nominal", int ir = 2);
    static float    getPurity(float val, string systag = "nominal", int ir = 2);
    static float    getPurityErrorLow(float low, float high, string systag = "nominal", int ir = 2);
    static float    getPurityErrorHigh(float low, float high, string systag = "nominal", int ir = 2);
    // Region-C purity from the same puritymaker.C fit ("combined_C" graph).
    static float    getPurityC(float low, float high, string systag = "nominal", int ir = 2);
    static float    getPurityCErrorLow(float low, float high, string systag = "nominal", int ir = 2);
    static float    getPurityCErrorHigh(float low, float high, string systag = "nominal", int ir = 2);
    // hists/purity_<systag>.root, one ana::rnames[ir] subdirectory per radius (written by puritymaker.C).
    static string   purityFilename(const string & systag);
    static Int_t findPtBin(double value);
    static Int_t findUnfoldXjBin(double value);
    // Keep a single overload: extra overloads made Cling's resolution non-deterministic and once
    // silently passed the BDT score as an int.
    static Int_t findabcdBin(double iso, double bdt, int bin = 0);
    static Int_t findxjBin(double value);
    static Int_t findBdtBin(double value);
    static Int_t findUnfoldBin(double xj, double pt);
    static Int_t findHadronBin(double value);
    static Int_t findEmfracBin(double value);

    // EM scale/resolution systematics as in PPG12 (sPH-ppg-2024-012, Sec. 5.1-5.2).
    // emscaleShift: +-1.48% cluster-pT scale (emscale_high/low).
    static constexpr float emscaleShift = 0.0148;
    // emResolutionSigma: extra Gaussian MC cluster smearing, sigma_extra = sqrt(max(0, sigma_data^2 -
    // sigma_MC^2)) with sigma/E = sqrt(p0^2/E + p1^2/E^2 + p2^2), at the truth pT (PPG12). Data parameters:
    // nominal (0.15,0.05,0.05), EMRhigh (0.13,0.08,0.08), EMRlow none. MC (0.185,0,0.040).
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
    // Data event timing cut, as the multiJet skim (|t_lead + 2| < 6 ns, |t_lead - t_sub| < 3 ns) with the
    // photon cluster as the leading object: |t_cluster - timingClusterCenter| < timingClusterHalfWidth and,
    // when a jet is stored, |t_cluster - t_jet| < timingDeltaMax (ns). MC has no timing.
    static constexpr double timingClusterCenter = -2;
    static constexpr double timingClusterHalfWidth = 6;
    static constexpr double timingDeltaMax = 3;
    static constexpr double radius = 93; 
    // Jets: |eta_jet| < etacut - R.
    static constexpr double etacut = 1.1;
    // Photons: |eta^gamma| < 0.7 for reco and truth (also the fiducial definition), as in PPG12;
    // beyond it the isolation cone leaves the EMCal.
    static constexpr double photonEtaCut = 0.7;
    static constexpr double etamin = -etacut;
    static constexpr double etamax = etacut;
    static constexpr double cluster_pt_cut = 10;

    static constexpr int nCalibBins = 7; // calib stages: uncalib, JES, JES+JER, +reweight, +photon smear, JERhigh, JERlow
    static constexpr int nPtBins = 5;
    // First (13-15 GeV) and last (35-100 GeV) pT bins are migration buffers, not reported.
    static constexpr int nPtBinsUsed = 3;
    static constexpr int firstUsedPtBin = 1;
    static constexpr int nUnfoldXjBins = 16;
    // findabcdBin's bin selector: 0 nominal; 1-5 each move one ABCD boundary (systagAbcdBinArr in
    // unfolder.cc; combined into one Purity systematic in draw_systematics.C).
    static constexpr int nIsoBdtBins = 6;
    static constexpr int nBdtBins = 4;
    static constexpr int nxjBins = 3;
    static constexpr int nJetR = 7;
    static constexpr int n3jetBins = 2;
    static constexpr int nHadronBins = 3;
    static constexpr int nEmfracBins = 3;
    static constexpr int singleptlow = 15;
    static constexpr int singlepthigh = 35;

    // The systag list - the single place to add one. Read by unfold_allsys.C (histograms and in-situ
    // trees), insitu/run_grid.sh (keep its copy in sync) and draw_systematics.C. Generator/unfolding
    // variations (herwig, niter, prior) are not systags.
    static const vector<string> systags;

    // Two-point (high/low) pairs: each contributes its signed per-bin value to the up or down total
    // (CLAUDE.md). Other systags are symmetric. Pairs must be listed here; names are not parsed.
    static const vector<pair<string,string>> asymmetricSystagPairs;
    static constexpr double ptBins[nPtBins+1] = {13,15,20,25,35,100};
    // ptBins[firstUsedPtBin..] as its own array (Cling had trouble with pointer offsets into static
    // constexpr members).
    static constexpr double ptBinsUsed[nPtBinsUsed+1] = {15,20,25,35};
    static constexpr double unfoldXjBins[nUnfoldXjBins+1] = {0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.5,1.7,2.0};
    // ABCD boundaries per variation (0 nominal, 1 narrowBDT, 2 narrowISO, 3 narrowBDTbkg, 4 narrowISObkg,
    // 5 wideISObkg); each moves one boundary, as in PPG12 Sec. 5.3: signal side (bdtGoodLow/isoBins)
    // or background side (bdtBadLow/isoBinsHigh). Each is a symmetric source.
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
    static constexpr double jet_calib_pt_cut[nJetR] = {5,5,5,5,5,5,5};
    // threejet: veto if the third jet's pT (recoil-jet pT definition) is above this. Reco only, so the
    // unfolded observable is unchanged.
    static constexpr double thirdJetPtCut = 5;

    // In-situ JES p_a per radius (grid_insitu.C purity-corrected fit); Data jet pT is divided by it.
    // Each systag with its own scan uses its own p_a (jesBySystag), so its JES effect stays inside that
    // systematic. jes_high/jes_low shift nominal p_a by its statistical error only (jes_high uses
    // jesStatErrLow: smaller p_a, larger correction).
    // Generated by insitu/draw_jes_summary.C (after run_grid.sh) down to END jesBySystag - do not edit
    // by hand; rebuild after it runs. jesSystagNames must match ana::systags.
    static constexpr double jesNominal[nJetR]      = {0.9034, 0.9033, 0.9158, 0.9116, 0.9172, 0.9151, 0.9289};
    static constexpr double jesStatErrLow[nJetR]   = {0.0056, 0.0041, 0.0083, 0.0070, 0.0065, 0.0072, 0.0079};
    static constexpr double jesStatErrHigh[nJetR]  = {0.0112, 0.0103, 0.0067, 0.0041, 0.0063, 0.0074, 0.0085};
    // BEGIN jesBySystag
    static constexpr int nJesSystags = 15;
    static constexpr const char * jesSystagNames[nJesSystags] = {"nominal", "JERhigh", "JERlow", "emscale_high", "emscale_low", "jes_high", "jes_low", "threejet", "narrowBDT", "narrowISO", "EMRhigh", "EMRlow", "narrowBDTbkg", "narrowISObkg", "wideISObkg"};
    static constexpr double jesBySystag[nJesSystags][nJetR] = {
      {0.9034, 0.9033, 0.9158, 0.9116, 0.9172, 0.9151, 0.9289}, // nominal
      {0.8946, 0.9001, 0.9072, 0.9060, 0.9131, 0.9114, 0.9276}, // JERhigh
      {0.9145, 0.9180, 0.9203, 0.9141, 0.9197, 0.9169, 0.9304}, // JERlow
      {0.9097, 0.9099, 0.9203, 0.9141, 0.9217, 0.9194, 0.9322}, // emscale_high
      {0.9006, 0.9033, 0.9101, 0.9064, 0.9124, 0.9092, 0.9250}, // emscale_low
      {0.9034, 0.9033, 0.9142, 0.9115, 0.9172, 0.9143, 0.9286}, // jes_high
      {0.9034, 0.9038, 0.9158, 0.9116, 0.9176, 0.9152, 0.9290}, // jes_low
      {0.8946, 0.9033, 0.9149, 0.9208, 0.9171, 0.9141, 0.9377}, // threejet
      {0.9013, 0.9093, 0.9229, 0.9194, 0.9188, 0.9152, 0.9300}, // narrowBDT
      {0.8973, 0.9026, 0.9101, 0.9059, 0.9139, 0.9150, 0.9235}, // narrowISO
      {0.9145, 0.9184, 0.9265, 0.9295, 0.9305, 0.9285, 0.9441}, // EMRhigh
      {0.9034, 0.9033, 0.9142, 0.9116, 0.9172, 0.9142, 0.9285}, // EMRlow
      {0.9097, 0.9096, 0.9187, 0.9186, 0.9218, 0.9225, 0.9322}, // narrowBDTbkg
      {0.9034, 0.9033, 0.9142, 0.9114, 0.9172, 0.9140, 0.9285}, // narrowISObkg
      {0.9034, 0.9038, 0.9158, 0.9116, 0.9181, 0.9151, 0.9288}, // wideISObkg
    };
    // END jesBySystag
    // p_a for a systag and radius; nominal (with a warning) if the systag has no entry.
    static double   jesForSystag(const string & systag, int ir);

  private:
};

#endif
