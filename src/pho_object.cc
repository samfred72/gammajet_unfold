#include "pho_object.h"

pho_object::pho_object() {}
pho_object::~pho_object() {}

int pho_object::get_showershape(float showershapes[], float pt) {
  bool e1133     = showershapes[5] < 0.98;
  bool wetacogx  = showershapes[3] < 0.6;
  bool et1       = 0.60 < showershapes[0] && showershapes[0] < 1;
  bool e3235     = 0.80 < showershapes[6] && showershapes[6] < 1;
  bool wetacogxt = 0.00 < showershapes[3] && showershapes[3] < 0.15 * pt;
  bool wphicogxt = 0.00 < showershapes[3] && showershapes[3] < 0.15 * pt;
  bool e1133t    = 0.40 < showershapes[5] && e1133;
  bool et1t      = 0.90 < showershapes[0] && showershapes[0] < 1;
  bool e3235t    = 0.92 < showershapes[6] && showershapes[6] < 1;

  bool isloose = e1133 && wetacogx && et1 && e3235 &&
    (wetacogxt + wphicogxt + e1133t + et1t + e3235t <= 3);
  bool istight = e1133 && wetacogx && et1 && e3235 &&
    wetacogxt && wphicogxt && e1133t && et1t && e3235t;

  int showershape = isloose * 1 + istight * 2; // 0: nothing, 1: loose, 2: tight
                                               
  return showershape;
}

void pho_object::print(bool s) {
  if (s) {
    std::cout << pt << " " << eta << " " << phi << " " << t << " " << iso4 << " " << bdt << std::endl;
  }
  else {
    std::cout << "Photon vals: "                << std::endl;
    std::cout << "pt:          " << pt          << std::endl; 
    std::cout << "eta:         " << eta         << std::endl; 
    std::cout << "phi:         " << phi         << std::endl; 
    std::cout << "e:           " << e           << std::endl; 
    std::cout << "time:        " << t           << std::endl; 
    std::cout << "iso3:        " << iso3        << std::endl; 
    std::cout << "iso4:        " << iso4        << std::endl; 
    std::cout << "bdt score:   " << bdt         << std::endl; 
    std::cout << "showershape: " << showershape << std::endl; 
  }
}
