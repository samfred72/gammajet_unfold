#include "object.h"

object::object() {}
object::~object() {}

float object::deltaR(const object& obj) {
  float dphi = deltaPhi(obj);
  float deta = eta - obj.eta;
  return std::sqrt(deta * deta + dphi * dphi);
}

float object::deltaPhi(const object& obj) {
  float dphi = std::fabs(phi - obj.phi);
  if (dphi > M_PI) dphi = 2 * M_PI - dphi;
  return dphi;
}
