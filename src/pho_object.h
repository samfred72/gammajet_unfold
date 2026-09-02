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

#endif // PHO_OBJECT_H
