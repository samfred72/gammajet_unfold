#include "/home/samson72/sphnx/gammajet_unfold/src/ana.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/drawer.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/pho_object.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/jet_object.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/insitu_utility.h"
#include "/home/samson72/sphnx/gammajet_unfold/src/unfold_utility.h"
// Explicit load; run interpreted, never with ACLiC "+" (sibling libgammajet.so collision).
R__LOAD_LIBRARY(libgammajet_unfold.so);

// PPG18 review round 2, issue 3: why is the R = 0.2 in-situ scan pinned at scanLow = 0.900?
//
// The production scan (insitu/grid_insitu.C) cannot look below 0.900: its data input
// (insitu/inputs/Data_*_insitu.root) only keeps events with raw xJ >= 0.90 x floor
// (unfolder.cc ispairedInsitu, floorScale = insitu_utility::scanLow). Here the data A/C
// events are rebuilt straight from trees/gammajet_Data.root with NO xJ pre-cut, so the
// same chi2 can be evaluated over pa = 0.70-1.10. Everything else is the production scan:
//   - data jet pT = jet_pt_calib[ir] (raw, no in-situ factor), same ABCD, same check_pair
//     cuts minus the floor (photon |eta|, jet |eta| < 1.1-R, dphi > 7pi/8)
//   - reference = insitu_utility::referenceMeans over Photon5+10+20 region A (unchanged)
//   - purities ana::getPurity/getPurityC per R, purityCorrectCoeffs, chi2 exactly as
//     grid_insitu.C's purity-corrected chi2
// Cross-check: with raw x >= 0.90 x floor, the rebuilt A/C counts must equal
// insitu_utility::cacheDataEvents on Data_nominal_insitu.root.
//
// Two scan variants per R:
//   "production"  floor applied to the SCALED xJ at every pa (as grid_insitu.C)
//   "fixed pop."  floor applied once at pa = 1, then the same events are scaled - the
//                 mean then scales exactly as 1/pa, so this isolates what the moving
//                 floor does to the fit
// Output: chi2 curves (pdf), minima and per-pT-bin preferred pa (log).

namespace {
  struct Ev { float pho_pt, jet_pt; int ptbin; };
}

void scan_extended_range()
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  const int irs[3] = {0, 1, 2};
  const int nR = 3, nPt = ana::nPtBinsUsed;
  const float paLo = 0.70, paHi = 1.10; const int nPa = 400;
  map<int,double> photon_scale = {{5,146359.3},{10,6944.675},{20,130.4461}}; // = grid_insitu.C
  const char * inDir = "/home/samson72/sphnx/gammajet_unfold/insitu/inputs";

  // ---------------- rebuild data A/C events from the tree ----------------
  vector<Ev> evA[nR], evC[nR];
  TFile * fin = TFile::Open("/home/samson72/sphnx/gammajet_unfold/trees/gammajet_Data.root", "read");
  TTree * t = (TTree*)fin->Get("towerntup");
  float vz, cluster_pt, cluster_e, cluster_eta, cluster_phi, cluster_time;
  float cluster_showershape[12], cluster_bdt_scores[11];
  float jet_pt_calib[7], jet_e[7], jet_eta[7], jet_phi[7], jet_emfrac[7];
  t->SetBranchStatus("*", 0);
  for (const char * b : {"vz","cluster_pt","cluster_e","cluster_eta","cluster_phi","cluster_time",
        "cluster_showershape","cluster_bdt_scores","jet_pt_calib","jet_e","jet_eta","jet_phi","jet_emfrac"})
    t->SetBranchStatus(b, 1);
  t->SetBranchAddress("vz", &vz);
  t->SetBranchAddress("cluster_pt", &cluster_pt);
  t->SetBranchAddress("cluster_e", &cluster_e);
  t->SetBranchAddress("cluster_eta", &cluster_eta);
  t->SetBranchAddress("cluster_phi", &cluster_phi);
  t->SetBranchAddress("cluster_time", &cluster_time);
  t->SetBranchAddress("cluster_showershape", cluster_showershape);
  t->SetBranchAddress("cluster_bdt_scores", cluster_bdt_scores);
  t->SetBranchAddress("jet_pt_calib", jet_pt_calib);
  t->SetBranchAddress("jet_e", jet_e);
  t->SetBranchAddress("jet_eta", jet_eta);
  t->SetBranchAddress("jet_phi", jet_phi);
  t->SetBranchAddress("jet_emfrac", jet_emfrac);
  for (Long64_t e = 0; e < t->GetEntries(); e++) {
    t->GetEntry(e);
    if (fabs(vz) > ana::vzcut) continue;
    pho_object pho(cluster_pt, cluster_e, cluster_eta, cluster_phi, cluster_showershape[10], cluster_showershape[11],
        cluster_time, cluster_bdt_scores[9], pho_object::get_showershape(cluster_showershape, cluster_pt));
    int iabcd = ana::findabcdBin(pho.iso4, pho.bdt, 0);
    if (iabcd != 0 && iabcd != 2) continue;
    if (fabs(pho.eta) > ana::etacut) continue;
    int ptbin = ana::findPtBin(pho.pt);
    if (ptbin < ana::firstUsedPtBin || ptbin >= ana::firstUsedPtBin + nPt) continue;
    for (int r = 0; r < nR; r++) {
      int ir = irs[r];
      if (jet_pt_calib[ir] <= 0) continue;
      jet_object jet(jet_pt_calib[ir], jet_e[ir], jet_eta[ir], jet_phi[ir], jet_emfrac[ir], 0, 0, 0);
      if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) continue;
      if (jet.deltaPhi(pho) < ana::oppcut) continue;
      Ev ev{pho.pt, jet_pt_calib[ir], ptbin - ana::firstUsedPtBin};
      (iabcd == 0 ? evA[r] : evC[r]).push_back(ev);
    }
  }

  // ---------------- scan ----------------
  drawer d("pythia", "nominal");
  string outdir = "/home/samson72/sphnx/gammajet_unfold/claude_checks/insitu_scan_edge/pdfs";
  gSystem->mkdir(outdir.c_str(), true);
  string pdf = outdir + "/scan_extended_range.pdf";
  TCanvas * c = new TCanvas("c", "", 1500, 550);
  c->SaveAs((pdf+"[").c_str());
  c->Divide(3, 1);
  TGraph * gProd[nR], * gFix[nR], * gBin[nR][nPt];

  for (int r = 0; r < nR; r++) {
    int ir = irs[r];
    float lowXj[nPt], purity[nPt], purityC[nPt], refMean[nPt], refMeanErr[nPt];
    for (int p = 0; p < nPt; p++) {
      lowXj[p] = insitu_utility::lowXjFloor(ir, ana::ptBinsUsed[p]);
      purity[p]  = ana::getPurity(ana::ptBinsUsed[p], ana::ptBinsUsed[p+1], "nominal", ir);
      purityC[p] = ana::getPurityC(ana::ptBinsUsed[p], ana::ptBinsUsed[p+1], "nominal", ir);
    }
    insitu_utility::referenceMeans({
        {insitu_utility::insituFilename(inDir, "Photon5",  "pythia", "nominal"), photon_scale[5]},
        {insitu_utility::insituFilename(inDir, "Photon10", "pythia", "nominal"), photon_scale[10]},
        {insitu_utility::insituFilename(inDir, "Photon20", "pythia", "nominal"), photon_scale[20]},
      }, 0, ir, refMean, refMeanErr, lowXj);

    // cross-check against the production input
    string dataFile = insitu_utility::insituFilename(inDir, "Data", "", "nominal");
    vector<DataEvent> prodA = insitu_utility::cacheDataEvents(dataFile.c_str(), 0, ir);
    vector<DataEvent> prodC = insitu_utility::cacheDataEvents(dataFile.c_str(), 2, ir);
    long nA09 = 0, nC09 = 0;
    for (auto & ev : evA[r]) if (ev.jet_pt/ev.pho_pt >= insitu_utility::scanLow*lowXj[ev.ptbin]) nA09++;
    for (auto & ev : evC[r]) if (ev.jet_pt/ev.pho_pt >= insitu_utility::scanLow*lowXj[ev.ptbin]) nC09++;
    printf("\n===== R = %.1f =====\ncross-check (raw xJ >= %.2f x floor): rebuilt A %ld / C %ld  vs  production input A %zu / C %zu\n",
           ana::JetRs[ir], insitu_utility::scanLow, nA09, nC09, prodA.size(), prodC.size());
    printf("reference <xJ> (gamma+jet MC, A): %.4f %.4f %.4f;  purities P_A %.3f %.3f %.3f  P_C %.3f %.3f %.3f\n",
           refMean[0], refMean[1], refMean[2], purity[0], purity[1], purity[2], purityC[0], purityC[1], purityC[2]);

    gProd[r] = new TGraph(); gFix[r] = new TGraph();
    for (int p = 0; p < nPt; p++) gBin[r][p] = new TGraph();
    for (int ia = 0; ia <= nPa; ia++) {
      float pa = paLo + ia*(paHi-paLo)/nPa;
      double chiProd = 0, chiFix = 0;
      for (int variant = 0; variant < 2; variant++) {
        vector<double> sA(nPt,0), sA2(nPt,0), sC(nPt,0), sC2(nPt,0); vector<int> nAc(nPt,0), nCc(nPt,0);
        for (int reg = 0; reg < 2; reg++) {
          for (auto & ev : (reg ? evC[r] : evA[r])) {
            float x = (ev.jet_pt/pa)/ev.pho_pt;
            float xFloor = variant == 0 ? x : ev.jet_pt/ev.pho_pt; // production: scaled; fixed: pa = 1
            if (xFloor < lowXj[ev.ptbin]) continue;
            if (reg) { sC[ev.ptbin] += x; sC2[ev.ptbin] += x*x; nCc[ev.ptbin]++; }
            else     { sA[ev.ptbin] += x; sA2[ev.ptbin] += x*x; nAc[ev.ptbin]++; }
          }
        }
        for (int p = 0; p < nPt; p++) {
          if (refMean[p] <= 0 || nAc[p] == 0 || nCc[p] == 0) continue;
          float cA, cC; unfold_utility::purityCorrectCoeffs(purity[p], purityC[p], nAc[p], nCc[p], cA, cC);
          double Ncorr = cA*nAc[p] - cC*nCc[p];
          if (Ncorr <= 0) continue;
          double mean = (cA*sA[p] - cC*sC[p])/Ncorr;
          double var = (cA*sA2[p] - cC*sC2[p])/Ncorr - mean*mean;
          double err = sqrt(std::max(var, 0.)/Ncorr);
          double diff = 1 - mean/refMean[p];
          double errt = sqrt(err*err/(refMean[p]*refMean[p]) + mean*mean*refMeanErr[p]*refMeanErr[p]/pow(refMean[p],4));
          double term = errt > 0 ? diff*diff/(errt*errt) : 0;
          if (variant == 0) { chiProd += term; gBin[r][p]->SetPoint(gBin[r][p]->GetN(), pa, term); }
          else chiFix += term;
        }
      }
      gProd[r]->SetPoint(ia, pa, chiProd); gFix[r]->SetPoint(ia, pa, chiFix);
    }
    auto minOf = [](TGraph * g, double lo, double hi, double & xmin, double & ymin) {
      ymin = 1e30; xmin = 0;
      for (int i = 0; i < g->GetN(); i++) if (g->GetX()[i] >= lo - 1e-6 && g->GetX()[i] <= hi + 1e-6 && g->GetY()[i] < ymin) { ymin = g->GetY()[i]; xmin = g->GetX()[i]; }
    };
    double xm, ym;
    minOf(gProd[r], insitu_utility::scanLow, insitu_utility::scanHigh, xm, ym);
    printf("production chi2, restricted to [%.2f, %.2f]: min at pa = %.3f (chi2 %.2f)\n", insitu_utility::scanLow, insitu_utility::scanHigh, xm, ym);
    minOf(gProd[r], paLo, paHi, xm, ym);
    printf("production chi2, full [%.2f, %.2f]:        min at pa = %.3f (chi2 %.2f)\n", paLo, paHi, xm, ym);
    minOf(gFix[r], paLo, paHi, xm, ym);
    printf("fixed-population chi2, full range:          min at pa = %.3f (chi2 %.2f)\n", xm, ym);
    for (int p = 0; p < nPt; p++) {
      minOf(gBin[r][p], paLo, paHi, xm, ym);
      printf("   pT %2.0f-%2.0f GeV alone (production): min at pa = %.3f (chi2 %.2f)\n", ana::ptBinsUsed[p], ana::ptBinsUsed[p+1], xm, ym);
    }

    // draw
    c->cd(r+1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.13); gPad->SetRightMargin(0.03); gPad->SetTopMargin(0.05);
    double ymax = 0; for (int i = 0; i < gProd[r]->GetN(); i++) ymax = std::max({ymax, gProd[r]->GetY()[i], gFix[r]->GetY()[i]});
    TH1D * fr = new TH1D(Form("frS%d", r), ";p_{a};#chi^{2} (purity corrected)", 1, paLo, paHi);
    fr->SetMinimum(0); fr->SetMaximum(std::min(ymax, 80.0)*1.2); fr->GetYaxis()->SetTitleSize(0.05); fr->GetXaxis()->SetTitleSize(0.05); fr->Draw();
    TBox * box = new TBox(insitu_utility::scanLow, 0, insitu_utility::scanHigh, fr->GetMaximum()); box->SetFillColor(kGray); box->SetFillStyle(1001); box->Draw();
    fr->Draw("axis same");
    gProd[r]->SetLineColor(kBlack); gProd[r]->SetLineWidth(2); gProd[r]->Draw("l same");
    gFix[r]->SetLineColor(kRed+1); gFix[r]->SetLineWidth(2); gFix[r]->SetLineStyle(2); gFix[r]->Draw("l same");
    int bc[3] = {kAzure+7, kGreen+2, kOrange+1};
    for (int p = 0; p < nPt; p++) { gBin[r][p]->SetLineColor(bc[p]); gBin[r][p]->SetLineWidth(1); gBin[r][p]->Draw("l same"); }
    TLegend * l = new TLegend(0.45, 0.55, 0.97, 0.78); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.035);
    l->AddEntry(gProd[r], "production (floor on scaled x_{J})", "l");
    l->AddEntry(gFix[r], "fixed population (floor at p_{a}=1)", "l");
    for (int p = 0; p < nPt; p++) l->AddEntry(gBin[r][p], Form("%.0f-%.0f GeV term", ana::ptBinsUsed[p], ana::ptBinsUsed[p+1]), "l");
    l->AddEntry(box, "production scan range", "f");
    l->Draw();
    d.drawAll({"p+p Run24 Data vs Pythia8 #gamma+jet"}, {Form("Jet R=%.1f", ana::JetRs[ir]), "rebuilt from tree, no x_{J} pre-cut"}, .17, .9, 12, gPad->GetWh()*0.8);
  }
  c->SaveAs(pdf.c_str());
  c->SaveAs((pdf+"]").c_str());
}
