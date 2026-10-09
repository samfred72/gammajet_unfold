#ifndef PURITY_UTILITY_H
#define PURITY_UTILITY_H

#include "ana.h"
#include <string>
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
using namespace std;

// Photon purity from the ABCD regions, shared by macros/puritymaker.C (from the unfolder's
// hclusterpt_abcd histograms) and insitu/grid_insitu.C (from the insitu trees at each trial p_a).
class purity_utility {
  public:
    // Leakage-corrected ABCD solve per ana::ptBins bin (h: Data counts A-D, f: MC leakage fractions
    // f^X = N_sig^X/N_sig^A), with a 10^4-toy bootstrap from a fresh default-seeded TRandom3, so equal
    // inputs give identical results. Returns P_A (median, 16th/84th percentiles; "combined"); P_C in
    // *graphCOut ("combined_C"). write: also write the bootstrap histograms, graphs and erf fit into
    // the current directory (puritymaker.C's file layout); otherwise they are deleted.
    static TGraphAsymmErrors * combine(TH1D * h[], TH1D * f[], TGraphAsymmErrors ** graphCOut = nullptr, bool write = true);

    // MC truth-matched leakage fractions for one systag and radius, as puritymaker.C. f must hold 4
    // entries; f[0] is 1.
    static void leakageFractions(const string & systag, int ir, TH1D * f[]);

    // Data ABCD counts vs cluster pT (ana::ptBins) from an insitu tree at jet scale pa, with
    // unfolder.cc's pairing: jet pT = jet_pt/pa above ana::jet_calib_pt_cut and the x_J floor of
    // unfolder::check_pair (the tree already has the photon, eta and dphi cuts), plus the threejet
    // veto when the tree defers it (thirdjet_pt >= 0). Reproduces hclusterpt_abcd<ir>_<j> exactly at
    // the pa the unfolding used (claude_checks/insitu_selfconsistent/validate_tree_purity.C).
    static void dataCountsFromTree(const string & filename, int ir, float pa, TH1D * h[]);
};

#endif // PURITY_UTILITY_H
