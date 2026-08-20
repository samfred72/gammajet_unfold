#include "jet_object.h"

jet_object::jet_object() {}
jet_object::~jet_object() {}

void jet_object::print(bool s) {
  if (s) {
    std::cout << pt << " " << eta << " " << phi << " " << t << " " << emfrac << std::endl;
  }
  else {
    std::cout << "Jet vals:"               << std::endl;
    std::cout << "pt:          " << pt     << std::endl; 
    std::cout << "eta:         " << eta    << std::endl; 
    std::cout << "phi:         " << phi    << std::endl; 
    std::cout << "e:           " << e      << std::endl; 
    std::cout << "time:        " << t      << std::endl; 
    std::cout << "efrac:       " << emfrac << std::endl; 
  }
}
