#ifndef OBJECT_H
#define OBJECT_H

#include <iostream>
#include <cmath>
#include <TMath.h>
#include <vector>

class object {
  public:
    object();
    object(float pt_, float e_, float eta_, float phi_, float t_)
      : e(e_),
      p(pt_ * TMath::CosH(eta_)),
      et(e_ / TMath::CosH(eta_)),
      pt(pt_),
      eta(eta_),
      phi(phi_),
      t(t_)
  {}

    virtual ~object(); // Virtual destructor for safe polymorphism

    float deltaR(const object& obj);
    float deltaPhi(const object& obj);

    float e = 0;
    float p = 0;
    float et = 0;
    float pt = 0;
    float eta = 0;
    float phi = 0;
    float t = 0;
};
#endif // OBJECT_H
