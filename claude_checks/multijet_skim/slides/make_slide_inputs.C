// Inputs for the multijet slides, from the new pythia trees:
//  1. eta-phi event displays (R=0.4, Jet20) of events the current analysis selects (individually smeared
//     jets 2 and 3 >= 7 GeV, leading >= 20, eta/dphi cuts) that have exactly two truth jets >= 7 GeV;
//  2. per radius, the fraction of selected events with exactly two truth jets >= 7 GeV, per sample and
//     for the cross-section-weighted Jet12+20+30 (stitched on the leading truth jet).
#include <vector>
#include <algorithm>

static float dphi(float a, float b) { float d = std::fabs(a-b); return d > M_PI ? 2*M_PI-d : d; }

// displaysOnly: only Jet20, stop after the displays (no fraction table)
void make_slide_inputs(const char * dir = "/sphenix/tg/tg01/jets/samfred/multiJet_full_hadded", int nDisplays = 6, bool displaysOnly = false) {
  gStyle->SetOptStat(0);
  const int nR = 7; int radii[nR] = {2,3,4,5,6,7,8};
  const char * samples[3] = {"Jet12", "Jet20", "Jet30"};
  const double weight[3] = {1.4903e+06/2.5298e+03, 6.2623e+04/2.5298e+03, 1};
  // stitching on the leading truth jet with Jet12 from 0 (Jet5/Jet8 not used)
  const double lo[3][nR] = {{0,0,0,0,0,0,0}, {20,21,21,27,29,32,34}, {30,31,32,38,41,45,47}};
  const double hi[3][nR] = {{20,21,21,27,29,32,34}, {30,31,32,38,41,45,47}, {1e9,1e9,1e9,1e9,1e9,1e9,1e9}};

  double nSel[3][nR] = {}, nTwo[3][nR] = {}, wSel[nR] = {}, wTwo[nR] = {};
  int nDrawn = 0;
  TCanvas * c = new TCanvas("c", "", 620, 1050);

  for (int s = 0; s < 3; s++) {
    if (displaysOnly && s != 1) continue;
    TFile f(Form("%s/multijet_pythia_%s.root", dir, samples[s]));
    TTree * t = (TTree*)f.Get("ttree");
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

    Long64_t n = t->GetEntries();
    for (Long64_t e = 0; e < n; e++) {
      if (displaysOnly && nDrawn >= nDisplays) break;
      t->GetEntry(e);
      if (std::fabs(z) > 60) continue;
      for (int i = 0; i < nR; i++) {
        auto & p = *sm[i]; float R = radii[i]/10.f;
        if (p.size() < 3 || tpt[i]->empty()) continue;
        float leadTruth = *std::max_element(tpt[i]->begin(), tpt[i]->end());
        if (leadTruth < lo[s][i] || leadTruth >= hi[s][i]) continue;
        std::vector<int> idx(p.size()); for (size_t k = 0; k < p.size(); k++) idx[k] = k;
        std::partial_sort(idx.begin(), idx.begin()+3, idx.end(), [&](int a, int b){ return p[a] > p[b]; });
        int j0 = idx[0], j1 = idx[1], j2 = idx[2];
        if (p[j0] < 20 || p[j1] < 7 || p[j2] < 7 || p[j1] >= 30 || p[j2] >= 30) continue;
        if (std::fabs(etad[i]->at(j0)) > 1.1-R || std::fabs(etad[i]->at(j1)) > 1.1-R || std::fabs(etad[i]->at(j2)) > 1.1-R) continue;
        if (dphi(phi[i]->at(j0), phi[i]->at(j1)) < 3*M_PI/4 || dphi(phi[i]->at(j0), phi[i]->at(j2)) < M_PI/2) continue;
        int nTruth = 0; for (float pt : *tpt[i]) if (pt >= 7) nTruth++;
        nSel[s][i]++; wSel[i] += weight[s];
        if (nTruth != 2) continue;
        nTwo[s][i]++; wTwo[i] += weight[s];

        if (s == 1 && radii[i] == 4 && nDrawn < nDisplays) {
          c->Clear(); c->SetLeftMargin(0.15); c->SetRightMargin(0.05); c->SetTopMargin(0.07); c->SetBottomMargin(0.08);
          TH2D * fr = new TH2D(Form("fr%d", nDrawn), ";#eta;#phi", 10, -1.5, 1.5, 10, -M_PI, M_PI);
          fr->GetXaxis()->SetTitleSize(0.06); fr->GetYaxis()->SetTitleSize(0.06); fr->GetYaxis()->SetTitleOffset(0.9);
          fr->GetXaxis()->SetLabelSize(0.045); fr->GetYaxis()->SetLabelSize(0.045); fr->GetXaxis()->SetTitleOffset(0.6);
          fr->Draw();
          TLine * acc1 = new TLine(-(1.1-R), -M_PI, -(1.1-R), M_PI), * acc2 = new TLine(1.1-R, -M_PI, 1.1-R, M_PI);
          for (auto l : {acc1, acc2}) { l->SetLineStyle(3); l->SetLineColor(kGray+1); l->Draw(); }
          const int roleColor[3] = {kBlack, kBlue+1, kRed+1};
          const char * roleName[3] = {"leading", "jet 2", "jet 3"};
          int roles[3] = {j0, j1, j2};
          for (size_t k = 0; k < cal[i]->size(); k++) {
            int role = -1; for (int r = 0; r < 3; r++) if ((int)k == roles[r]) role = r;
            TEllipse * el = new TEllipse(eta[i]->at(k), phi[i]->at(k), R, R);
            el->SetFillStyle(0); el->SetLineWidth(role >= 0 ? 3 : 1);
            el->SetLineColor(role >= 0 ? roleColor[role] : kGray+1); el->Draw();
            float ly = phi[i]->at(k) + R + 0.08;
            if (ly > M_PI - 0.55) ly = phi[i]->at(k) - R - 0.55; // keep the label inside the frame
            TLatex * lab = new TLatex(std::min(eta[i]->at(k) - 0.45, 0.25), ly,
              role >= 0 ? Form("#splitline{%s}{%.1f #rightarrow %.1f GeV}", roleName[role], cal[i]->at(k), sm[i]->at(k))
                        : Form("%.1f GeV", cal[i]->at(k)));
            lab->SetTextSize(0.045); lab->SetTextColor(role >= 0 ? roleColor[role] : kGray+2); lab->Draw();
          }
          for (size_t k = 0; k < tpt[i]->size(); k++) {
            if (tpt[i]->at(k) < 3) continue;
            TMarker * m = new TMarker(teta[i]->at(k), tphi[i]->at(k), 29);
            m->SetMarkerSize(tpt[i]->at(k) >= 7 ? 2.6 : 1.4); m->SetMarkerColor(kGreen+2); m->Draw();
            TLatex * lab = new TLatex(std::min(teta[i]->at(k) + 0.05, 0.7), tphi[i]->at(k) - 0.32, Form("truth %.1f", tpt[i]->at(k)));
            lab->SetTextSize(0.042); lab->SetTextColor(kGreen+3); lab->Draw();
          }
          float rec = std::sqrt(std::pow(sm[i]->at(j1)*std::cos(phi[i]->at(j1)) + sm[i]->at(j2)*std::cos(phi[i]->at(j2)), 2) +
                                std::pow(sm[i]->at(j1)*std::sin(phi[i]->at(j1)) + sm[i]->at(j2)*std::sin(phi[i]->at(j2)), 2));
          TLatex head; head.SetNDC(); head.SetTextSize(0.045);
          head.DrawLatex(0.15, 0.95, Form("recoil = %.1f GeV,  x_{j} = %.2f", rec, sm[i]->at(j0)/rec));
          c->SaveAs(Form("event_display_%d.pdf", nDrawn));
          c->SaveAs(Form("event_display_%d.png", nDrawn));
          nDrawn++;
        }
      }
    }
  }
  printf("Fraction of selected events with exactly two truth jets >= 7 GeV (current analysis selection)\n");
  printf("R     Jet12          Jet20          Jet30          xsec-weighted MC\n");
  for (int i = 0; i < nR; i++) {
    printf("0.%d", radii[i]);
    for (int s = 0; s < 3; s++) printf("   %5.1f%% (%6.0f)", nSel[s][i] ? 100*nTwo[s][i]/nSel[s][i] : 0., nSel[s][i]);
    printf("   %5.1f%%\n", wSel[i] ? 100*wTwo[i]/wSel[i] : 0.);
  }
}
