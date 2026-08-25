#ifndef REWEIGHT_UTILITY_H
#define REWEIGHT_UTILITY_H

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include <iostream>

// Loads the Data/MC v_z and cluster-pT weights written by
// reweight/make_vz_pt_reweight.C and evaluates them for a given event. Usage:
//   Reweighter rw;
//   float w = rw.GetWeight(tu.vz, tu.cluster_pt);
//
// v_z uses hVzWeight directly (one lookup per bin, no fit - stats are fine across the
// whole |v_z| < 60 cm range). Cluster pT uses hPtWeight's first bin (10-11 GeV) as a flat
// value below 11 GeV, and the fPtWeight "expo" fit (fit over the well-populated 11-25 GeV
// range, see make_vz_pt_reweight.C) for pT >= 11 - the raw bins above ~25 GeV are too
// low-statistics to use directly.
class Reweighter {
  public:
    Reweighter(const char * filename = "/home/samson72/sphnx/gammajet_unfold/reweight/vz_pt_reweight.root") {
      f = TFile::Open(filename, "read");
      if (!f || f->IsZombie()) {
        std::cout << "Reweighter: could not open " << filename << std::endl;
        return;
      }
      hVzWeight = (TH1D*)f->Get("hVzWeight");
      hPtWeight = (TH1D*)f->Get("hPtWeight");
      fPtWeight = (TF1*)f->Get("fPtWeight");

      ptLowWeight = hPtWeight->GetBinContent(1);
      ptLowEdge = hPtWeight->GetXaxis()->GetBinUpEdge(1); // 11 GeV
    }

    double GetVzWeight(double vz) {
      int b = hVzWeight->GetXaxis()->FindFixBin(vz);
      if (b < 1 || b > hVzWeight->GetNbinsX()) return 1; // outside the |v_z| < 60 cm input range
      return hVzWeight->GetBinContent(b);
    }

    double GetPtWeight(double pt) {
      if (pt < ptLowEdge) return ptLowWeight;
      return fPtWeight->Eval(pt);
    }

    double GetWeight(double vz, double pt) {
      return GetVzWeight(vz) * GetPtWeight(pt);
    }

    TFile * f = nullptr;
    TH1D * hVzWeight = nullptr;
    TH1D * hPtWeight = nullptr;
    TF1 * fPtWeight = nullptr;
    double ptLowWeight = 1;
    double ptLowEdge = 11;
};

#endif // REWEIGHT_UTILITY_H
