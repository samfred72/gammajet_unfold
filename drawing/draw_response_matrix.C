#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
// The original gammajet project builds its OWN, differently-laid-out drawer/ana classes
// into /home/samson72/root/lib/libgammajet.so, sitting on the same library search path
// as this project's libgammajet_unfold.so. Without forcing which one loads first, ROOT's
// implicit symbol autoload can bind drawer/ana calls to the wrong (mismatched-layout)
// library and segfault - explicit load removes the ambiguity.
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Plots the raw response matrix used as unfolding input: hxjresponse%i, written by
// src/unfolder.cc:469-470 as RooUnfoldResponse::Hresponse() for the Photon5/10/20
// combination. Per RooUnfoldResponse.h:66, axes are (x,y)=(measured,truth), flattened
// over pT x xJ (with under/overflow) bins on both axes - see draw_covariance.C's
// matrixIndex() for the exact flat-index convention, reused here only for the pT-bin
// gridlines below. This macro does no unfolding/closure/systematics - just a direct look
// at the migration matrix itself.

const int ir = 2; // nominal jet radius index (R=0.4), matching every other macro in this directory

void draw_response_matrix(string systag = "nominal") {
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);

  drawer d("pythia", systag);
  string pdfPath = Form("/home/samson72/sphnx/gammajet_unfold/pdfs/draw_response_matrix_%s.pdf", systag.c_str());

  TH2D * hResponse = d.get2d(Form("hxjresponse%i", ir), 1);

  // NaN/Inf poisons ROOT's axis auto-ranging and silently produces a blank canvas.
  for (int bx = 0; bx <= hResponse->GetNbinsX()+1; bx++) {
    for (int by = 0; by <= hResponse->GetNbinsY()+1; by++) {
      double v = hResponse->GetBinContent(bx, by);
      if (!std::isfinite(v)) hResponse->SetBinContent(bx, by, 0);
    }
  }

  TCanvas * c = new TCanvas("c", "", 900, 800);
  c->cd();
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.12);
  gPad->SetTicks(1, 1);

  hResponse->SetTitle("");
  hResponse->GetXaxis()->SetTitle("Measured bin (flattened p_{T}^{#gamma} #times x_{J#gamma})");
  hResponse->GetYaxis()->SetTitle("Truth bin (flattened p_{T}^{#gamma} #times x_{J#gamma})");
  hResponse->GetYaxis()->SetTitleOffset(1.4);
  hResponse->GetZaxis()->SetTitle("Entries");
  hResponse->Draw("colz");

  // pT-bin boundary gridlines: the flattened axes stack nUnfoldXjBins+2 (under/overflow
  // included) xJ bins per pT bin, ana::nPtBins pT bins total - same flat-index convention
  // as draw_covariance.C's matrixIndex(). Purely visual, to show the block-diagonal
  // structure (migration is expected to stay within a pT block, only crossing pT
  // boundaries near the block edges).
  int nFlatPerPt = ana::nUnfoldXjBins + 2;
  int nFlatTotal = nFlatPerPt * ana::nPtBins;
  for (int ipt = 1; ipt < ana::nPtBins; ipt++) {
    double edge = hResponse->GetXaxis()->GetBinLowEdge(ipt*nFlatPerPt + 1);
    TLine * lx = new TLine(edge, hResponse->GetYaxis()->GetXmin(), edge, hResponse->GetYaxis()->GetXmax());
    TLine * ly = new TLine(hResponse->GetXaxis()->GetXmin(), edge, hResponse->GetXaxis()->GetXmax(), edge);
    lx->SetLineStyle(3);
    ly->SetLineStyle(3);
    lx->Draw("same");
    ly->Draw("same");
  }

  d.drawAll({"Pythia8 #gamma+jet MC"}, {Form("Jet R=%.1f", ana::JetRs[ir])}, .18, .95, 14, 800);
  c->SaveAs(pdfPath.c_str());
  cout << "Done. Wrote " << pdfPath << endl;
}
