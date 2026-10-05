#ifndef REWEIGHT_UTILITY_H
#define REWEIGHT_UTILITY_H

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include <iostream>

// Data/MC weights from reweight/make_vz_pt_reweight.C: v_z from hVzWeight bins; cluster pT from
// hPtWeight's first bin below 11 GeV and the fPtWeight expo fit (11-25 GeV) above.
class Reweighter {
  public:
    Reweighter(const char * filename = ana::path("reweight/vz_pt_reweight.root")) {
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
