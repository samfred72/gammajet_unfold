R__LOAD_LIBRARY(libgammajet_unfold.so);
#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
void abcdcomp(){
  TFile*fd=TFile::Open("/home/samson72/sphnx/gammajet_unfold/hists/Data_nominal_unfolding.root");
  TFile*fj=TFile::Open("/home/samson72/sphnx/gammajet_unfold/hists/Jet12_long_pythia_nominal_unfolding.root");
  drawer d("pythia","nominal");
  for(int ir: {0,2}){
    TH1D*hd[4],*hj[4],*hjt[4],*hpt[4];
    for(int k=0;k<4;k++){hd[k]=(TH1D*)fd->Get(Form("hclusterpt_abcd%d_%d",ir,k)); hj[k]=(TH1D*)fj->Get(Form("hclusterpt_abcd%d_%d",ir,k)); hjt[k]=(TH1D*)fj->Get(Form("hclusterpt_abcd_truthmatched%d_%d",ir,k)); hpt[k]=d.get(Form("hclusterpt_abcd_truthmatched%d_%d",ir,k),1);}
    printf("\nR=%.1f\n%-7s | %-38s | %-38s | %s\n",ana::JetRs[ir],"pT","Data A:B:C:D (norm to D)","Jet12 A:B:C:D (norm to D)","Jet12 bkg-only (untagged): C/D  B/D  A*D/(B*C)");
    for(int b=2;b<=4;b++){
      double D=hd[3]->GetBinContent(b), J=hj[3]->GetBinContent(b);
      double bk[4]; for(int k=0;k<4;k++) bk[k]=hj[k]->GetBinContent(b)-hjt[k]->GetBinContent(b);
      printf("%3.0f-%-3.0f | %6.2f %5.2f %5.2f 1 (A/C=%.2f, B/D=%.2f) | %6.2f %5.2f %5.2f 1 (A/C=%.2f, B/D=%.2f) | %.2f %.2f %.2f\n",
        hd[0]->GetBinLowEdge(b),hd[0]->GetXaxis()->GetBinUpEdge(b),
        hd[0]->GetBinContent(b)/D,hd[1]->GetBinContent(b)/D,hd[2]->GetBinContent(b)/D,hd[0]->GetBinContent(b)/hd[2]->GetBinContent(b),hd[1]->GetBinContent(b)/D,
        hj[0]->GetBinContent(b)/J,hj[1]->GetBinContent(b)/J,hj[2]->GetBinContent(b)/J,hj[0]->GetBinContent(b)/hj[2]->GetBinContent(b),hj[1]->GetBinContent(b)/J,
        bk[2]/bk[3],bk[1]/bk[3],bk[0]*bk[3]/(bk[1]*bk[2]));
      printf("        leakage (Photon MC truth-matched) b=%.3f c=%.3f d=%.3f | Jet12 tagged leakage b=%.3f c=%.3f d=%.3f\n",
        hpt[1]->GetBinContent(b)/hpt[0]->GetBinContent(b),hpt[2]->GetBinContent(b)/hpt[0]->GetBinContent(b),hpt[3]->GetBinContent(b)/hpt[0]->GetBinContent(b),
        hjt[1]->GetBinContent(b)/hjt[0]->GetBinContent(b),hjt[2]->GetBinContent(b)/hjt[0]->GetBinContent(b),hjt[3]->GetBinContent(b)/hjt[0]->GetBinContent(b));
    }
  }
}
