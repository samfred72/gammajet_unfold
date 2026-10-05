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
    // p_a for a systag and radius; nominal (with a warning) if the systag has no entry.
    static double   jesForSystag(const string & systag, int ir);

  private:
};

#endif
