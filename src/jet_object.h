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

#endif // JET_OBJECT_H
