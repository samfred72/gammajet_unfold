#ifndef PHO_OBJECT_H
#define PHO_OBJECT_H

#include "object.h"
#include <cmath>
#include <TMath.h>
#include <vector>

class pho_object : public object {
  public:
    pho_object();
    pho_object(float pt_, float e_, float eta_, float phi_,
        float i3, float i4, float t_,
        float bdt_, int showershape_)
      : object(pt_, e_, eta_, phi_, t_),
      iso3(i3), iso4(i4), bdt(bdt_), showershape(showershape_)
  {}

    ~pho_object();
    static int get_showershape(float showershapes[], float pt);
    void print(bool s = 0);

    float iso3 = 0;
    float iso4 = 0;
    float bdt = 0;
    int showershape = 0;
};

//--------------------------------------
// Helper functions to create vectors
//--------------------------------------


inline std::vector<pho_object> make_clusters(
    const std::vector<float>& pt,
    const std::vector<float>& e,
    const std::vector<float>& eta,
    const std::vector<float>& phi,
    const std::vector<std::vector<float>>& showershapes = {},
    const std::vector<float>& t = {},
    const std::vector<std::vector<float>>& bdt = {}) 
{
  std::vector<pho_object> v;
  for (size_t i = 0; i < pt.size(); i++) {
    if (!t.empty()) {
      float iso3 = showershapes.at(8).at(i) + showershapes.at(9).at(i) + showershapes.at(10).at(i);
      float iso4 = showershapes.at(11).at(i) + showershapes.at(12).at(i) + showershapes.at(13).at(i);

      bool e1133     = showershapes.at(5).at(i) < 0.98;
      bool wetacogx  = showershapes.at(3).at(i) < 0.6;
      bool et1       = 0.60 < showershapes.at(0).at(i) && showershapes.at(0).at(i) < 1;
      bool e3235     = 0.80 < showershapes.at(6).at(i) && showershapes.at(6).at(i) < 1;
      bool wetacogxt = 0.00 < showershapes.at(3).at(i) && showershapes.at(3).at(i) < 0.15 * pt[i];
      bool wphicogxt = 0.00 < showershapes.at(3).at(i) && showershapes.at(3).at(i) < 0.15 * pt[i];
      bool e1133t    = 0.40 < showershapes.at(5).at(i) && e1133;
      bool et1t      = 0.90 < showershapes.at(0).at(i) && showershapes.at(0).at(i) < 1;
      bool e3235t    = 0.92 < showershapes.at(6).at(i) && showershapes.at(6).at(i) < 1;

      bool isloose = e1133 && wetacogx && et1 && e3235 &&
        (wetacogxt + wphicogxt + e1133t + et1t + e3235t <= 3);
      bool istight = e1133 && wetacogx && et1 && e3235 &&
        wetacogxt && wphicogxt && e1133t && et1t && e3235t;

      int showershape = isloose * 1 + istight * 2; // 0: nothing, 1: loose, 2: tight

      v.emplace_back(pt[i], e[i], eta[i], phi[i], iso3, iso4,
          t[i] * 17.6, bdt.empty() ? 0 : bdt[9][i], showershape);
    } else {
      v.emplace_back(pt[i], e[i], eta[i], phi[i], 0, 0, 0, 0, 0);
    }
  }
  return v;
}

#endif // PHO_OBJECT_H
