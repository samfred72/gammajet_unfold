#include "../../src/ana.h"
#include "../../src/drawer.h"
#include "../../src/pho_object.h"
#include "../../src/jet_object.h"
// Explicit load: the sibling gammajet project's libgammajet.so has same-named classes. Run
// interpreted (root -b -q draw_timing_distributions.C), never with ACLiC "+".
R__LOAD_LIBRARY(libgammajet_unfold.so);

// Timing distributions of region-A photon-jet pairs in data, one page each:
//   t_cluster, t_jet, t_cluster - t_MBD, t_jet - t_MBD, t_cluster - t_jet,
// then each versus pT (2D, with the mean per pT bin overlaid): cluster quantities versus
// photon pT, jet quantities versus corrected jet pT, and t_cluster - t_jet versus both.
// Pairing follows unfolder::check_pair: |v_z| < ana::vzcut, photon |eta| < ana::photonEtaCut,
// jet |eta| < ana::etacut - R, jet pT = jet_pt_calib/jesNominal > jet_calib_pt_cut,
// Delta phi > ana::oppcut, and the xJ floor.
// The trees hold only objects inside the treemaking windows (data only, CaloAna.cc):
// photon cluster 0 < t_MBD - t_cluster < 4 ns, jet -2 < t_MBD - t_jet < 5 ns. In the
// object-minus-MBD convention used here those are -4 < t_cl - t_MBD < 0 and
// -5 < t_jet - t_MBD < 2; they are drawn dashed. Jets with no tower above 0.1 GeV have
// t_jet = -999 and fail the jet window, so none are stored.

namespace {
  // NaN/Inf poisons ROOT's auto-ranging (CLAUDE.md) - zero such bins before drawing.
  void sanitize(TH1 * h) {
    for (int b = 0; b <= h->GetNcells(); b++) {
      if (!std::isfinite(h->GetBinContent(b)) || !std::isfinite(h->GetBinError(b))) {
        h->SetBinContent(b, 0); h->SetBinError(b, 0);
      }
    }
  }
  // Visible x range: the central (1 - 2*tail) of the entries, padded by 10%.
  void autoRange(TH1D * h, double tail = 0.001) {
    double q[2], p[2] = {tail, 1 - tail};
    h->GetQuantiles(2, q, p);
    double pad = 0.1*(q[1] - q[0]);
    h->GetXaxis()->SetRangeUser(q[0] - pad, q[1] + pad);
  }
}

void draw_timing_distributions(int ir = 2)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TH1::SetDefaultSumw2();

  TFile * fin = TFile::Open(ana::path("trees/gammajet_Data.root"), "read");
  TTree * t = (TTree*)fin->Get("towerntup");

  float vz, mbd_time, cluster_pt, cluster_e, cluster_eta, cluster_phi, cluster_time;
  float cluster_showershape[12], cluster_bdt_scores[11];
  float jet_pt_calib[7], jet_e[7], jet_eta[7], jet_phi[7], jet_emfrac[7], jet_time[7];
  t->SetBranchStatus("*", 0);
  for (const char * b : {"vz","mbd_time","cluster_pt","cluster_e","cluster_eta","cluster_phi",
        "cluster_time","cluster_showershape","cluster_bdt_scores","jet_pt_calib","jet_e",
        "jet_eta","jet_phi","jet_emfrac","jet_time"}) t->SetBranchStatus(b, 1);
  t->SetBranchAddress("vz", &vz);
  t->SetBranchAddress("mbd_time", &mbd_time);
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
  t->SetBranchAddress("jet_time", jet_time);

  const double phoLow = 15, phoHigh = 35; // reported photon-pT range (ana::ptBins adds unfolding buffer bins)
  // Fine bins over a wide range; the drawn range is set from the quantiles.
  const int nb = 800; const double lo = -40, hi = 40;
  enum { kCl, kJet, kClMbd, kJetMbd, kClJet, nVar };
  const char * xtitle[nVar] = {"t_{cluster} [ns]", "t_{jet} [ns]", "t_{cluster} - t_{MBD} [ns]",
                               "t_{jet} - t_{MBD} [ns]", "t_{cluster} - t_{jet} [ns]"};
  // Treemaking windows in this convention (none on the raw times or on cluster - jet).
  const double winLo[nVar] = {0, 0, -ana::thighcut, -5, 0};
  const double winHi[nVar] = {0, 0, -ana::tlowcut,   2, 0};
  const bool hasWin[nVar]  = {false, false, true, true, false};
  TH1D * h[nVar];
  for (int v = 0; v < nVar; v++) h[v] = new TH1D(Form("h%d", v), Form(";%s;pairs / 0.1 ns", xtitle[v]), nb, lo, hi);
  // 2D versus pT: x = photon pT (1 GeV bins) or corrected jet pT (variable bins), y = time.
  const int nPho = 20;
  const int nJpt = 9;
  const double jptEdges[nJpt+1] = {5, 8, 11, 14, 17, 20, 24, 28, 35, 50};
  struct Vs { int var; bool jetPt; TH2D * h2; };
  vector<Vs> vs = {{kCl, false, nullptr}, {kClMbd, false, nullptr}, {kJet, true, nullptr},
                   {kJetMbd, true, nullptr}, {kClJet, false, nullptr}, {kClJet, true, nullptr}};
  for (size_t i = 0; i < vs.size(); i++) {
    const char * xt = vs[i].jetPt ? "corrected jet p_{T} [GeV]" : "p_{T}^{#gamma} [GeV]";
    vs[i].h2 = vs[i].jetPt
      ? new TH2D(Form("h2_%zu", i), Form(";%s;%s", xt, xtitle[vs[i].var]), nJpt, jptEdges, nb, lo, hi)
      : new TH2D(Form("h2_%zu", i), Form(";%s;%s", xt, xtitle[vs[i].var]), nPho, phoLow, phoHigh, nb, lo, hi);
  }

  long nPairs = 0;
  // Multijet-style timing cuts applied to the photon-jet pair (cluster in the role of the leading jet):
  // multijet/analysis.cc:433-436 |t_lead| < 6 ns, |t_lead - t_sub| < 3 ns; DijetTreeMaker.cc:773-778
  // skim |t_lead + 2| < 6 ns. Counted on the stored (already windowed) pairs.
  long nDt3 = 0, nJet6 = 0, nJet6off = 0, nCl6 = 0, nCl6off = 0, nMjAna = 0, nMjSkim = 0;
  Long64_t n = t->GetEntries();
  for (Long64_t e = 0; e < n; e++) {
    t->GetEntry(e);
    if (fabs(vz) > ana::vzcut) continue;
    if (cluster_pt < phoLow || cluster_pt >= phoHigh) continue;
    if (jet_pt_calib[ir] <= 0) continue; // no in-window jet stored for this radius

    pho_object pho(cluster_pt, cluster_e, cluster_eta, cluster_phi,
        cluster_showershape[10], cluster_showershape[11], cluster_time,
        cluster_bdt_scores[9], pho_object::get_showershape(cluster_showershape, cluster_pt));
    if (ana::findabcdBin(pho.iso4, pho.bdt, 0) != 0) continue; // region A
    if (fabs(pho.eta) > ana::photonEtaCut) continue;

    float jetPt = jet_pt_calib[ir] / ana::jesNominal[ir];
    if (jetPt <= ana::jet_calib_pt_cut[ir]) continue;
    jet_object jet(jetPt, jet_e[ir], jet_eta[ir], jet_phi[ir], jet_emfrac[ir], 0, 0, jet_time[ir]);
    if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) continue;
    if (jet.deltaPhi(pho) < ana::oppcut) continue;
    int ptbin = ana::findPtBin(pho.pt);
    float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(ana::jet_calib_pt_cut[ir]/ana::ptBins[ptbin])+1];
    if (jetPt / pho.pt < lowbin) continue;

    double val[nVar];
    val[kCl] = cluster_time;
    val[kJet] = jet_time[ir];
    val[kClMbd] = cluster_time - mbd_time;
    val[kJetMbd] = jet_time[ir] - mbd_time;
    val[kClJet] = cluster_time - jet_time[ir];
    for (int v = 0; v < nVar; v++) h[v]->Fill(val[v]);
    bool dt3 = fabs(val[kClJet]) < 3;
    nDt3 += dt3;
    nJet6 += fabs(val[kJet]) < 6;  nJet6off += fabs(val[kJet] + 2) < 6;
    nCl6 += fabs(val[kCl]) < 6;    nCl6off += fabs(val[kCl] + 2) < 6;
    nMjAna += dt3 && fabs(val[kCl]) < 6;
    nMjSkim += dt3 && fabs(val[kCl] + 2) < 6;
    for (auto & x : vs) x.h2->Fill(x.jetPt ? jetPt : pho.pt, val[x.var]);
    nPairs++;
  }
  printf("R=%.1f region-A pairs: %ld\n", ana::JetRs[ir], nPairs);
  auto frac = [&](long k) { return 100.0*k/std::max(nPairs, 1L); };
  printf("  multijet-style cuts on these pairs (fraction kept):\n");
  printf("    |t_cluster - t_jet| < 3 ns                    %.2f%%\n", frac(nDt3));
  printf("    |t_jet| < 6 ns / |t_jet + 2| < 6 ns          %.2f%% / %.2f%%\n", frac(nJet6), frac(nJet6off));
  printf("    |t_cluster| < 6 ns / |t_cluster + 2| < 6 ns  %.2f%% / %.2f%%\n", frac(nCl6), frac(nCl6off));
  printf("    analysis-style (|t_cl| < 6, |dt| < 3)        %.2f%%\n", frac(nMjAna));
  printf("    skim-style (|t_cl + 2| < 6, |dt| < 3)        %.2f%%\n", frac(nMjSkim));
  for (int v = 0; v < nVar; v++)
    printf("  %-28s mean %7.3f  RMS %6.3f  underflow %.0f  overflow %.0f\n", xtitle[v], h[v]->GetMean(), h[v]->GetRMS(),
           h[v]->GetBinContent(0), h[v]->GetBinContent(nb+1));

  drawer d("pythia", "nominal"); // ctor args only pick files for other helpers; labels come from drawAll
  string phoLabel = Form("%.0f GeV < p_{T}^{#gamma} < %.0f GeV", phoLow, phoHigh);
  string rLabel = Form("Jet R=%.1f, p_{T}^{jet} > %.0f GeV", ana::JetRs[ir], ana::jet_calib_pt_cut[ir]);
  string outdir = ana::path("claude_checks/jet_timing/pdfs");
  gSystem->mkdir(outdir.c_str(), true);
  string pdf = outdir + Form("/draw_timing_distributions_r%02d.pdf", (int)(ana::JetRs[ir]*10+0.5));
  TCanvas * c = new TCanvas("c", "", 1200, 550);
  for (int v = 0; v < nVar; v++) {
    c->Clear();
    c->Divide(2, 1);
    sanitize(h[v]);
    h[v]->SetLineColor(kBlack); h[v]->SetLineWidth(2);
    autoRange(h[v]);
    for (int ip = 1; ip <= 2; ip++) {
      TPad * p = (TPad*)c->cd(ip);
      p->SetLeftMargin(0.15); p->SetBottomMargin(0.13); p->SetRightMargin(0.04); p->SetTopMargin(0.05);
      p->SetLogy(ip == 2);
      TH1D * hd = (TH1D*)h[v]->Clone(Form("%s_p%d", h[v]->GetName(), ip));
      hd->SetMinimum(ip == 2 ? 0.5 : 0);
      hd->SetMaximum(ip == 2 ? 50*hd->GetMaximum() : 1.6*hd->GetMaximum());
      hd->Draw("hist");
      if (hasWin[v]) {
        for (double x : {winLo[v], winHi[v]}) {
          TLine * l = new TLine(x, hd->GetMinimum(), x, ip == 2 ? hd->GetMaximum()/50 : hd->GetMaximum()/1.6);
          l->SetLineStyle(2); l->SetLineColor(kRed+1); l->SetLineWidth(2); l->Draw();
        }
      }
      vector<string> feats = {phoLabel, rLabel, "Region A, paired",
                              Form("mean %.2f ns, RMS %.2f ns", h[v]->GetMean(), h[v]->GetRMS())};
      if (hasWin[v]) feats.push_back("dashed: treemaking window");
      d.drawAll({"p+p Run24 Data"}, feats, .19, .9, 14, p->GetWh()*0.8);
    }
    c->SaveAs((pdf + (v == 0 ? "(" : "")).c_str());
  }

  // 2D pages: counts (colz) with the mean per pT bin, linear and log z.
  for (size_t i = 0; i < vs.size(); i++) {
    TH2D * h2 = vs[i].h2;
    sanitize(h2);
    h2->RebinY(2); // 0.2 ns (before setting the range: rebinning resets it)
    // y range from the matching 1D distribution's quantiles, with headroom for the labels
    double q[2], pr[2] = {0.001, 0.999};
    h[vs[i].var]->GetQuantiles(2, q, pr);
    double span = q[1] - q[0];
    h2->GetYaxis()->SetRangeUser(q[0] - 0.1*span, q[1] + 1.0*span);
    TProfile * prof = h2->ProfileX(Form("%s_prof", h2->GetName()));
    prof->SetLineColor(kRed+1); prof->SetMarkerColor(kRed+1); prof->SetMarkerStyle(20); prof->SetMarkerSize(0.8); prof->SetLineWidth(2);
    c->Clear();
    c->Divide(2, 1);
    for (int ip = 1; ip <= 2; ip++) {
      TPad * p = (TPad*)c->cd(ip);
      p->SetLeftMargin(0.15); p->SetBottomMargin(0.13); p->SetRightMargin(0.14); p->SetTopMargin(0.05);
      p->SetLogz(ip == 2);
      h2->Draw("colz");
      prof->Draw("same");
      vector<string> feats = {phoLabel, rLabel, "Region A, paired", "red: mean per p_{T} bin"};
      if (hasWin[vs[i].var]) feats.push_back("window cuts the tails (treemaking)");
      d.drawAll({"p+p Run24 Data"}, feats, .19, .9, 14, p->GetWh()*0.8);
    }
    c->SaveAs((pdf + (i == vs.size()-1 ? ")" : "")).c_str());
    printf("  mean %-28s vs %-10s:", xtitle[vs[i].var], vs[i].jetPt ? "jet pT" : "photon pT");
    for (int b = 1; b <= prof->GetNbinsX(); b++)
      if (prof->GetBinEntries(b) > 0) printf(" [%.0f-%.0f] %.2f", prof->GetXaxis()->GetBinLowEdge(b), prof->GetXaxis()->GetBinUpEdge(b), prof->GetBinContent(b));
    printf("\n");
  }
  printf("Wrote %s\n", pdf.c_str());
}
