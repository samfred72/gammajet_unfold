#ifndef INSITU_UTILITY_H
#define INSITU_UTILITY_H

#include "ana.h"
#include <vector>
#include <string>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"

// Helpers shared by the insitu/*.C JES scans.
struct DataEvent { float pho_pt, jet_pt; int ptbin; };
// Multijet balance event: leading pT, recoil jets' pT/phi, weight (1 for Data).
struct MultijetEvent { float lead, sl, slphi, ssl, sslphi, w; };

class insitu_utility {
  public:
    // JES scan window [scanLow, scanHigh) in scanN steps of 1e-4. unfolder.cc keeps in-situ events
    // down to scanLow (ispairedInsitu floor), so lowering it needs a rerun of the unfolding.
    static constexpr float scanLow = 0.80;
    static constexpr float scanHigh = 1.00;
    static constexpr int scanN = 2000;

    // Data events from unfolder.cc's insitutree for one ABCD region and jet radius. restrictToUsed
    // keeps ana::ptBinsUsed (ptbin re-indexed from 0); false keeps every ana::ptBins bin (the unfolded
    // scans need the full vector for pT migration). The MC "weight" branch is not read (always 1 in Data).
    static vector<DataEvent> cacheDataEvents(const char * filename, int abcdSelect, int ir, bool restrictToUsed = true);

    // Low-x_J floor for radius ir and pT-bin lower edge ptLow: the x_J-bin upper edge above
    // jet_calib_pt_cut[ir]/ptLow, exactly as unfolder::check_pair. Every scan applies it at each trial pa.
    static double lowXjFloor(int ir, double ptLow);

    // insitutree file name as written by unfolder.cc; sim = "" for Data.
    static string insituFilename(const char * insitu_dir, const char * trigger,
        const char * sim, const string & systag);

    // Points where chi2 first exceeds minchisq + 1 on each side of the minimum.
    static void findError(TGraph * g, int ibest, float minchisq, float & errLow, float & errHigh);

    // Mean x_J vs photon pT over ana::ptBinsUsed.
    static TGraphErrors * meanGraph(const float mean[], const float err[], const char * name);

    // Ratio of two mean arrays, errors added in quadrature.
    static TGraphErrors * ratioGraph(const float meanNum[], const float errNum[],
        const float meanDen[], const float errDen[], const char * name);

    // Same label block as drawer::drawAll, without constructing a drawer (which opens many files).
    static void drawSPhenixLabel(vector<string> samples, vector<string> features,
        float drawx, float drawy, int fontsize, float csize);

    // ----- Shared by the grid_insitu*.C scans -----

    // Sum fine bins [startBin, nBinsForFit) in groups of groupSize (shape method, pT bin 2).
    static vector<double> coarsenSum(const vector<double> & fine, int startBin, int nBinsForFit, int groupSize);
    static vector<double> coarsenQuadrature(const vector<double> & fineErr, int startBin, int nBinsForFit, int groupSize);

    // Weighted mean x_J and error per used pT bin over (file, weight) MC samples - the reference the
    // Data scan is compared to. Same filters as cacheDataEvents.
    static void referenceMeans(const vector<pair<string,double>> & samples, int abcdSelect, int ir,
        float refMean[], float refMeanErr[], const float lowXj[]);

    // Reference x_J shape (bin fractions) and per-bin errors per used pT bin.
    static void referenceShape(const vector<pair<string,double>> & samples, int abcdSelect, int ir,
        vector<vector<double>> & refFrac, vector<vector<double>> & refFracErr, const float lowXj[]);

    // Region-A mean x_J at trial scale pa (display where the fit uses another criterion).
    static void computeRegionAMeans(const vector<DataEvent> & dataA, float pa, float mean[], float err[], const float lowXj[]);

    // Purity-corrected mean x_J at trial pa (A and C scaled alike; purity fixed).
    static void computeCorrectedMeans(const vector<DataEvent> & dataA, const vector<DataEvent> & dataC,
        float pa, const float purity[], const float purityC[], float mean[], float err[], const float lowXj[]);

    // Data x_J histogram per pT bin at trial pa. nBins = ana::nPtBinsUsed, or ana::nPtBins for the
    // unfolded scans.
    static vector<TH1D*> buildXjByPtBin(const vector<DataEvent> & data, float pa, int nBins,
        const char * prefix, const float lowXj[]);

    // MC reference x_J histogram per used pT bin (weights as referenceMeans).
    static vector<TH1D*> buildMCXjByPtBin(const vector<pair<string,double>> & samples, int abcdSelect, int ir,
        const char * prefix, const float lowXj[]);

    // Two-purity correction per pT bin (unfold_utility::purityCorrect). The unfolded scans pass
    // zero purity errors.
    static vector<TH1D*> purityCorrectByPtBin(const vector<TH1D*> & hA, const vector<TH1D*> & hC, int nBins,
        const float purity[], const float purityErrLow[], const float purityErrHigh[],
        const float purityC[], const float purityCErrLow[], const float purityCErrHigh[],
        const char * prefix);

    // ----- Multijet balance (grid_insitu.C combined mode) -----
    // Reads multijet/analysis.cc's per-radius trees from multijet_analysis_<sim>.root, so the multijet
    // selection and MC weighting live only there:
    //   Data: ttree_data_r<10R>; MC: ttree_<Jet5..Jet30>_r<10R>_<RECO|HIGH|LOW>
    // B = pT,lead / |pT,sub + pT,subsub|. A constant JES cancels in B, so it constrains only the slope pb
    // of f(pT) = pa + pb*pT. The Data trees keep events down to analysis.cc's jesFloor, so the pT cuts
    // below are applied after each jet is divided by f.
    static constexpr int nMultijetPtBins = 4;
    static constexpr double multijetPtBins[nMultijetPtBins+1] = {20, 25, 30, 35, 50}; // = analysis.cc pTBins
    static constexpr float multijetBalanceLow = 0.4, multijetBalanceHigh = 2.65; // = analysis.cc hxj range
    // = analysis.cc cuts: leading jet, jets 2 and 3 each, recoil |pT2 + pT3|
    static constexpr float multijetLeadCut = 20, multijetRecoilJetCut = 7, multijetRecoilCut = 14;
    static constexpr int nMultijetMCSamples = 5;
    static constexpr const char * multijetMCSamples[nMultijetMCSamples] = {"Jet5", "Jet8", "Jet12", "Jet20", "Jet30"};

    static int findMultijetPtBin(double leadPt);
    static string multijetAnalysisFilename(const char * multijet_dir, const char * sim);
    // JER variant for a systag: HIGH (JERhigh), LOW (JERlow), else RECO.
    static string multijetSysName(const string & systag);
    // Empty vector, with a warning, if the file or tree is missing.
    static vector<MultijetEvent> cacheMultijetEvents(const char * filename, const string & treename);
    // B with each jet divided by f at its own pT, and the corrected leading pT; -1 if the corrected
    // jets fail the pT cuts.
    // applyCuts=false for the MC trees: analysis.cc already selected them, and in its two-truth-jet
    // events jets 2 and 3 are stored scaled to the smeared recoil.
    static double multijetBalance(const MultijetEvent & ev, double pa, double pb, double & leadCorr, bool applyCuts = true);
    // Weighted mean B and Kish error per leading-pT bin, for B in [multijetBalanceLow, multijetBalanceHigh).
    static void multijetMeans(const vector<MultijetEvent> & events, double pa, double pb,
        float mean[], float err[], bool applyCuts = true);
};

#endif // INSITU_UTILITY_H
