// Writes a GAMMAJET_JES_TABLE file holding the compiled ana.h p_a: every systag row of jesBySystag.
// run_iter.sh uses it as the default start table. Run without GAMMAJET_JES_TABLE set.
#include "../src/ana.h"
R__LOAD_LIBRARY(libgammajet_unfold.so);

void write_table(string outname) {
  FILE * f = fopen(outname.c_str(), "w");
  fprintf(f, "# in-situ p_a per systag (R=0.2 ... 0.8), from ana.h jesBySystag; format: see ana.cc jesTableOverride\n");
  for (int i = 0; i < ana::nJesSystags; i++) {
    fprintf(f, "%s", ana::jesSystagNames[i]);
    for (int ir = 0; ir < ana::nJetR; ir++) fprintf(f, " %.4f", ana::jesBySystag[i][ir]);
    fprintf(f, "\n");
  }
  fclose(f);
  printf("wrote %s\n", outname.c_str());
}
