#ifndef JET_OBJECT_H
#define JET_OBJECT_H

#include "object.h"
#include <cmath>
#include <TMath.h>
#include <vector>

class jet_object : public object {
  public:
    jet_object();
    jet_object(float pt_, float e_, float eta_, float phi_,
        float efrac, float ifrac, float ofrac, float t_)
      : object(pt_, e_, eta_, phi_, t_),
      emfrac(efrac),
      ihfrac(ifrac),
      ohfrac(ofrac)
  {}

    ~jet_object();
    void print(bool s = 0);

    float emfrac = 0;
    float ihfrac = 0;
    float ohfrac = 0;
};

//--------------------------------------
// Helper functions to create vectors
//--------------------------------------

inline std::vector<jet_object> make_jets(
    const std::vector<float>& pt,
    const std::vector<float>& e,
    const std::vector<float>& eta,
    const std::vector<float>& phi,
    const std::vector<float>& emfrac = {},
    const std::vector<float>& ihfrac = {},
    const std::vector<float>& ohfrac = {},
    const std::vector<float>& t = {}) 
{
  std::vector<jet_object> v;
  for (size_t i = 0; i < pt.size(); i++) {
    if (!t.empty()) {
      v.emplace_back(pt[i], e[i], eta[i], phi[i],
          emfrac.empty() ? 0 : emfrac[i],
          ihfrac.empty() ? 0 : ihfrac[i],
          ohfrac.empty() ? 0 : ohfrac[i],
          t[i] * 17.6);
    } else {
      v.emplace_back(pt[i], e[i], eta[i], phi[i], 0, 0, 0, 0);
    }
  }
  return v;
}

#endif // JET_OBJECT_H
