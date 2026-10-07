// Recoil = |pT2 + pT3| (vector sum of the subleading and subsubleading jets), two ways, each selected on
// its OWN smeared recoil: smeared recoil >= recoilMin (no per-jet cuts on jets 2 and 3), leading jet's
// smeared pT >= 20, |eta| < 1.1-R for the three jets, dphi12 > 3pi/4, dphi13 > pi/2. Events with exactly
// two truth jets >= truthMin.
//   A (current): jets ranked by their individually smeared pT (jet_pt_smear_reco); smeared jets 2+3 added.
//   B (proposed): jets ranked by calib pT (leading jet keeps its own smeared pT); jets 2+3 added unsmeared,
//      then the sum smeared once with the nominal JER template width at the pT of the truth jet matched to
//      the recoil (the one not matched to the leading jet).
// Unweighted (single sample).
#include <vector>
#include <algorithm>

static float dphi(float a, float b) { float d = std::fabs(a-b); return d > M_PI ? 2*M_PI-d : d; }

void recoil_smear_compare_recoilcut(const char * file = "/sphenix/tg/tg01/jets/samfred/multiJet_full_hadded/multijet_pythia_Jet20.root",
                          const char * templ = "/sphenix/user/samfred/projects/gammajet/treemaking/macros/jerband_smearing_templates.root",
                          const char * outpdf = "recoil_smear_compare_recoilcut.pdf", float truthMin = 7, float recoilMin = 15) {
  gStyle->SetOptStat(0);
  TFile ft(templ);
  TH1D * width = (TH1D*)ft.Get("h_jer_smear_r04_pileup_EMfracJES_nominal");
  width->SetDirectory(0);
  TRandom3 rnd(12345);

  TFile f(file);
  TTree * t = (TTree*)f.Get("ttree");
  const int nR = 7; int radii[nR] = {2,3,4,5,6,7,8};
  const int nPt = 6; const double ptb[nPt+1] = {20,25,30,35,40,50,60};
  std::vector<float> *cal[nR], *sm[nR], *eta[nR], *etad[nR], *phi[nR], *tpt[nR], *teta[nR], *tphi[nR];
  float z;
  t->SetBranchStatus("*", 0);
  for (int i = 0; i < nR; i++) {
    cal[i] = sm[i] = eta[i] = etad[i] = phi[i] = tpt[i] = teta[i] = tphi[i] = nullptr;
    auto addr = [&](const char * n, std::vector<float> ** v) { t->SetBranchStatus(Form(n, radii[i]), 1); t->SetBranchAddress(Form(n, radii[i]), v); };
    addr("jet_pt_calib_%d", &cal[i]); addr("jet_pt_smear_reco_%d", &sm[i]);
    addr("jet_eta_%d", &eta[i]); addr("jet_eta_det_%d", &etad[i]); addr("jet_phi_%d", &phi[i]);
    addr("truth_jet_pt_%d", &tpt[i]); addr("truth_jet_eta_%d", &teta[i]); addr("truth_jet_phi_%d", &tphi[i]);
  }
  t->SetBranchStatus("mbd_vertex_z", 1); t->SetBranchAddress("mbd_vertex_z", &z);

  // [radius][method A/B][subset: all selected, 2-truth-jet events][pT bin]
  TH1D * hRec[nR][2][2][nPt], * hXj[nR][2][2][nPt];
  for (int i = 0; i < nR; i++) for (int m = 0; m < 2; m++) for (int s = 0; s < 2; s++) for (int k = 0; k < nPt; k++) {
    hRec[i][m][s][k] = new TH1D(Form("rec_%d_%d_%d_%d",i,m,s,k), ";recoil |p_{T,2}+p_{T,3}| (GeV);events", 50, 0, 50);
    hXj[i][m][s][k]  = new TH1D(Form("xj_%d_%d_%d_%d",i,m,s,k), ";x_{j};events", 45, 0.4, 2.65);
  }

  // passes the selection with these jets (i0 leading) and these pT values
  auto select = [&](int i, int i0, int i1, int i2, float p0, float p1, float p2) {
    float R = radii[i]/10.f;
    if (p0 < 20 || p1 < 7 || p2 < 7 || p1 >= 30 || p2 >= 30) return false;
    for (int j : {i0, i1, i2}) if (std::fabs(etad[i]->at(j)) > 1.1-R) return false;
    return dphi(phi[i]->at(i0), phi[i]->at(i1)) >= 3*M_PI/4 && dphi(phi[i]->at(i0), phi[i]->at(i2)) >= M_PI/2;
  };
  auto geom = [&](int i, int i0, int i1, int i2, float p0) {
    float R = radii[i]/10.f;
    if (p0 < 20) return false;
    for (int j : {i0, i1, i2}) if (std::fabs(etad[i]->at(j)) > 1.1-R) return false;
    return dphi(phi[i]->at(i0), phi[i]->at(i1)) >= 3*M_PI/4 && dphi(phi[i]->at(i0), phi[i]->at(i2)) >= M_PI/2;
  };
  auto vsum = [](float a, float pa, float b, float pb) {
    float x = a*std::cos(pa) + b*std::cos(pb), y = a*std::sin(pa) + b*std::sin(pb); return std::sqrt(x*x + y*y);
  };
  auto top3 = [](const std::vector<float> & p) {
    std::vector<int> idx(p.size()); for (size_t k = 0; k < p.size(); k++) idx[k] = k;
    std::partial_sort(idx.begin(), idx.begin()+3, idx.end(), [&](int a, int b){ return p[a] > p[b]; });
    return idx;
  };
  auto bin = [&](float lead) { for (int k = 0; k < nPt; k++) if (lead >= ptb[k] && lead < ptb[k+1]) return k; return -1; };

  Long64_t n = t->GetEntries();
  for (Long64_t e = 0; e < n; e++) {
    t->GetEntry(e);
    if (std::fabs(z) > 60) continue;
    for (int i = 0; i < nR; i++) {
      if (cal[i]->size() < 3) continue;
      // truth jets above truthMin
      std::vector<int> tj;
      for (size_t k = 0; k < tpt[i]->size(); k++) if (tpt[i]->at(k) >= truthMin) tj.push_back(k);
      bool twoTruth = (tj.size() == 2);

      if (!twoTruth) continue;
      int k;
      // A: smeared ranking, smeared jets added, cut on that recoil
      auto a = top3(*sm[i]);
      float a0 = sm[i]->at(a[0]);
      float recA = vsum(sm[i]->at(a[1]), phi[i]->at(a[1]), sm[i]->at(a[2]), phi[i]->at(a[2]));
      if (geom(i, a[0], a[1], a[2], a0) && recA >= recoilMin && (k = bin(a0)) >= 0) {
        hRec[i][0][1][k]->Fill(recA); hXj[i][0][1][k]->Fill(a0/recA);
      }
      // B: calib ranking, unsmeared jets added, the sum smeared once, cut on that recoil
      auto b = top3(*cal[i]);
      float b0 = sm[i]->at(b[0]);
      auto dR = [&](int r, int tk) { float de = eta[i]->at(r) - teta[i]->at(tk), dp = dphi(phi[i]->at(r), tphi[i]->at(tk)); return std::sqrt(de*de + dp*dp); };
      int tLead = dR(b[0], tj[0]) <= dR(b[0], tj[1]) ? tj[0] : tj[1];
      int tRec = (tLead == tj[0]) ? tj[1] : tj[0];
      float ptTruth = tpt[i]->at(tRec);
      float recCal = vsum(cal[i]->at(b[1]), phi[i]->at(b[1]), cal[i]->at(b[2]), phi[i]->at(b[2]));
      float recB = rnd.Gaus(recCal, ptTruth * width->Interpolate(std::min(std::max(ptTruth, 5.01f), 79.9f)));
      if (geom(i, b[0], b[1], b[2], b0) && recB >= recoilMin && (k = bin(b0)) >= 0) {
        hRec[i][1][1][k]->Fill(recB); hXj[i][1][1][k]->Fill(b0/recB);
        hRec[i][0][0][k]->Fill(recCal); hXj[i][0][0][k]->Fill(b0/recCal); // B's events, unsmeared recoil
      }
    }
  }

  printf("Each method cut on its own smeared recoil >= %.0f GeV (two truth jets >= %.0f GeV): events A / B, <recoil>, <xj>\n", recoilMin, truthMin);
  printf("columns: B's events unsmeared | A smear-then-add | B add-then-smear\n");
  for (int i = 0; i < nR; i++) {
    printf("R=0.%d", radii[i]);
    for (int k = 0; k < 5; k++)
      printf(" | %2.0f-%2.0f: N=%6.0f/%6.0f  rec %5.2f %5.2f %5.2f  xj %.3f %.3f %.3f", ptb[k], ptb[k+1], hRec[i][0][1][k]->GetEntries(), hRec[i][1][1][k]->GetEntries(),
             hRec[i][0][0][k]->GetMean(), hRec[i][0][1][k]->GetMean(), hRec[i][1][1][k]->GetMean(),
             hXj[i][0][0][k]->GetMean(), hXj[i][0][1][k]->GetMean(), hXj[i][1][1][k]->GetMean());
    printf("\n");
  }

  // pages: per radius, recoil spectra (two-truth-jet events) then x_j (all selected events), A vs B, absolute counts
  TCanvas * c = new TCanvas("c", "", 1500, 900);
  c->Print(Form("%s[", outpdf));
  for (int i = 0; i < nR; i++) {
    for (int page = 0; page < 2; page++) {
      c->Clear(); c->Divide(3, 2);
      int s = 1;
      for (int k = 0; k < nPt; k++) {
        TVirtualPad * p = c->cd(k+1);
        p->SetLeftMargin(0.14); p->SetRightMargin(0.03); p->SetTopMargin(0.08); p->SetBottomMargin(0.13);
        TH1D * ha = page == 0 ? hRec[i][0][s][k] : hXj[i][0][s][k];
        TH1D * hb = page == 0 ? hRec[i][1][s][k] : hXj[i][1][s][k];
        ha->SetLineColor(kBlue+1); ha->SetLineWidth(2);
        hb->SetLineColor(kRed+1); hb->SetLineWidth(2);
        ha->SetMaximum(1.35*std::max(ha->GetMaximum(), hb->GetMaximum())); ha->SetMinimum(0);
        ha->GetXaxis()->SetTitleSize(0.05); ha->GetYaxis()->SetTitleSize(0.05); ha->GetYaxis()->SetTitleOffset(1.4);
        TH1D * h0 = page == 0 ? hRec[i][0][0][k] : hXj[i][0][0][k];
        h0->SetLineColor(kGray+1); h0->SetLineStyle(2); h0->SetLineWidth(2);
        ha->SetMaximum(1.35*std::max({ha->GetMaximum(), hb->GetMaximum(), h0->GetMaximum()}));
        ha->Draw("hist"); hb->Draw("hist same"); h0->Draw("hist same");
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.045);
        tx.DrawLatex(0.16, 0.94, Form("Jet20, R=0.%d, %.0f < p_{T}^{lead} < %.0f GeV, %s", radii[i], ptb[k], ptb[k+1],
                                     Form("smeared recoil #geq %.0f GeV", recoilMin)));
        TLegend * l = new TLegend(0.45, 0.68, 0.97, 0.88); l->SetBorderSize(0); l->SetFillStyle(0); l->SetTextSize(0.04);
        l->AddEntry(ha, Form("A: smear jets, then add (mean %.2f)", ha->GetMean()), "l");
        l->AddEntry(hb, Form("B: add, then smear (mean %.2f)", hb->GetMean()), "l");
        l->AddEntry(h0, Form("B events, unsmeared (mean %.2f)", h0->GetMean()), "l");
        l->Draw();
      }
      c->Print(outpdf);
    }
  }
  c->Print(Form("%s]", outpdf));
}
