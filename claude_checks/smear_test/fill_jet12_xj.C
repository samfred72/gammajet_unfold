// Plain ROOT. Jet12 (Pythia) multijet xJ = pT,1 / |pT,2 + pT,3| with three low-pT extrapolations of the
// JER smearing width, in two smearing configurations. Output jet12_xj_R<RR>.root, drawn by draw_jet12_xj.C.
//
// Width functions (nominal only), as in fill_jet12_extrap.C:
//   tmpl   h_jer_smear_r04_pileup_EMfracJES_nominal->Interpolate(pt) (the treemaker's; flat below 5.19 GeV)
//   lin15  template above 15 GeV, its tangent at 15 GeV below
//   lin5   template above 5.19 GeV, its first segment (slope at ~5 GeV) continued below
// plus "calib" (no smearing) as a reference.
//
// Configurations (follow multijet/analysis.cc, d14878d, line refs to that file):
//   perjet   every jet smeared on its own (--recoil-smear none): pT = calib + z calib w(calib), with the
//            treemaker's per-jet deviate z (jet_smear_z_R) and the width at the jet's own pT_calib - the
//            "_reco" flavour analysis.cc reads by default (no --truth-smear). Jets ranked by smeared pT.
//   truth2   the nominal configuration (--recoil-smear truth2, truth jets >= 7 GeV): in events where jets 2
//            and 3 (ranked by pT_calib) belong to one truth jet (analysis.cc:456-496, ported verbatim), they
//            are added unsmeared and the sum is smeared once, pT,23 -> sum + z_recoil pT_truth w(pT_truth),
//            with the stored per-event recoil_smear_z_R; the leading jet keeps its per-jet smear and
//            ranking is by pT_calib (analysis.cc:518-538). Other events as perjet. The width function is
//            swapped in the recoil smear too (the truth jet is >= 7 GeV, so lin5 only differs from the
//            template there through the per-jet smears).
//
// Selection (analysis.cc:430-575): |z_vtx| <= 60 cm, >= 3 jets, leading truth jet >= 14 GeV (Jet12 full
// efficiency at R=0.4; no upper stitching cut, so this is Jet12 alone over its whole efficient range),
// leading reco jet 20-35 GeV (35 = the Jet12 cap), jets 2 and 3 >= 7 GeV each (before the recoil smear),
// recoil >= 14 GeV, |eta_det| < 1.1 - R for all three, dphi12 > 3pi/4, dphi13 > pi/2. Weighted by the
// nominal reweighting, zvtx ratio x ratio_func_JERreco(lead) from ../../multijet/aux/ratio<R>_pythia.root.
// xJ binned as analysis.cc (45 bins, 0.4-2.65) in leading-pT bins 20-25, 25-30, 30-35 GeV.
//
// Storage floor (see fill_jet12_extrap.C): lin5 jets missing from the tree all lie below 6.4 GeV, under the
// 7 GeV jet 2/3 cut, so they cannot enter or be ranked into a selected event; the xJ here is complete.

float dPhiAbs(float a, float b) { float d = std::fabs(a - b); return d > M_PI ? 2*M_PI - d : d; }
float vecSum(float pt1, float phi1, float pt2, float phi2, bool wantPhi) {
  const float x = pt1*std::cos(phi1) + pt2*std::cos(phi2), y = pt1*std::sin(phi1) + pt2*std::sin(phi2);
  return wantPhi ? std::atan2(y, x) : std::sqrt(x*x + y*y);
}
std::vector<int> top3(const std::vector<float> &pt) {
  std::vector<int> idx(pt.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(), [&](int a, int b) { return pt[a] > pt[b]; });
  idx.resize(3, -1);
  return idx;
}

void fill_jet12_xj(int radius = 4,
                   const char *tree = "/projects/sali6653/PPG18_analysis/trees/multijet_pythia_Jet12.root",
                   const char *auxdir = "/projects/sali6653/PPG18_analysis/multijet/aux",
                   Long64_t maxEvents = -1) {
  TH1::SetDefaultSumw2();
  const float R = radius / 10.0f;
  const double truthThresholdJet12[9] = {0, 0, 12, 13, 14, 19, 22, 24, 25}; // analysis.cc truthThreshold, by radius
  const double truthLeadMin = truthThresholdJet12[radius];
  const float leadMin = 20, leadMax = 35, jetMin = 7, recoilMin = 14, truthJetMin = 7;
  const float SLDPHI = 3*M_PI/4, SSLDPHI = M_PI/2;
  const std::vector<float> ptBins = {20, 25, 30, 35};
  const int nPt = 3; // ptBins.size() - 1

  // widths
  TFile fjer("jerband_smearing_templates.root");
  TH1D *hNom = (TH1D*)fjer.Get("h_jer_smear_r04_pileup_EMfracJES_nominal");
  hNom->SetDirectory(0);
  fjer.Close();
  const double p15 = 15, w15 = hNom->Interpolate(p15);
  const double s15 = (hNom->Interpolate(p15 + 0.01) - hNom->Interpolate(p15 - 0.01)) / 0.02;
  const double c1 = hNom->GetBinCenter(1), w5 = hNom->GetBinContent(1);
  const double s5 = (hNom->GetBinContent(2) - hNom->GetBinContent(1)) / (hNom->GetBinCenter(2) - c1);
  const int nW = 4;
  const char *wNames[nW] = {"calib", "tmpl", "lin15", "lin5"};
  auto width = [&](int w, double pt) -> double {
    if (w == 1) return hNom->Interpolate(pt);
    if (w == 2) return pt >= p15 ? hNom->Interpolate(pt) : w15 + s15*(pt - p15);
    if (w == 3) return pt >= c1 ? hNom->Interpolate(pt) : w5 + s5*(pt - c1);
    return 0;
  };

  // reweighting (analysis.cc:243-256, 665-667)
  TFile *fr = TFile::Open(Form("%s/ratio%d_pythia.root", auxdir, radius));
  TF1 *fit = (TF1*)fr->Get("ratio_func_JERreco");
  TH1D *zr = (TH1D*)fr->Get(Form("hratio_zvtx%d_pythia", radius));
  if (!fit || !zr) { printf("missing reweighting objects in %s\n", fr->GetName()); return; }

  const int nCfg = 2;
  const char *cfgNames[nCfg] = {"perjet", "truth2"};
  TH1D *hxj[nCfg][nW][nPt];
  for (int c = 0; c < nCfg; c++) for (int w = 0; w < nW; w++) for (int k = 0; k < nPt; k++)
    hxj[c][w][k] = new TH1D(Form("hxj_%s_%s_%.0f", cfgNames[c], wNames[w], ptBins[k]), ";x_{J};Weighted events", 45, 0.4, 2.65);

  TFile *f = TFile::Open(tree);
  TTree *t = (TTree*)f->Get("ttree");
  std::vector<float> *calib = nullptr, *eta = nullptr, *etadet = nullptr, *phi = nullptr, *z = nullptr, *stored = nullptr,
                     *tpt = nullptr, *teta = nullptr, *tphi = nullptr, *jtruth = nullptr;
  float vz = 0, recoilZ = 0;
  std::vector<TString> active = {TString::Format("jet_pt_calib_%d", radius), TString::Format("jet_eta_%d", radius),
      TString::Format("jet_eta_det_%d", radius), TString::Format("jet_phi_%d", radius), TString::Format("jet_smear_z_%d", radius),
      TString::Format("jet_pt_smear_reco_%d", radius), TString::Format("truth_jet_pt_%d", radius),
      TString::Format("truth_jet_eta_%d", radius), TString::Format("truth_jet_phi_%d", radius),
      TString::Format("jet_truth_pt_%d", radius), TString::Format("recoil_smear_z_%d", radius), TString("mbd_vertex_z")};
  t->SetBranchStatus("*", 0);
  for (const TString &b : active) t->SetBranchStatus(b, 1);
  t->SetBranchAddress(active[0], &calib);  t->SetBranchAddress(active[1], &eta);
  t->SetBranchAddress(active[2], &etadet); t->SetBranchAddress(active[3], &phi);
  t->SetBranchAddress(active[4], &z);      t->SetBranchAddress(active[5], &stored);
  t->SetBranchAddress(active[6], &tpt);    t->SetBranchAddress(active[7], &teta);
  t->SetBranchAddress(active[8], &tphi);   t->SetBranchAddress(active[9], &jtruth);
  t->SetBranchAddress(active[10], &recoilZ); t->SetBranchAddress(active[11], &vz);

  const Long64_t n = (maxEvents > 0) ? std::min(maxEvents, t->GetEntries()) : t->GetEntries();
  Long64_t nJets = 0, nMismatch = 0, nRecoilEvt = 0, nTruthOK = 0;
  double nSel[nCfg][nW] = {{0}};
  for (Long64_t e = 0; e < n; e++) {
    t->GetEntry(e);
    if (std::fabs(vz) > 60) continue;
    if (calib->size() <= 2 || tpt->empty()) continue;
    if (*std::max_element(tpt->begin(), tpt->end()) < truthLeadMin) continue;
    nTruthOK++;

    // per-jet smears for every width function; validate tmpl against the stored _reco smear
    const size_t nj = calib->size();
    std::vector<float> ptw[nW];
    for (int w = 0; w < nW; w++) ptw[w].resize(nj);
    for (size_t j = 0; j < nj; j++) {
      const double pc = calib->at(j);
      ptw[0][j] = pc;
      for (int w = 1; w < nW; w++) ptw[w][j] = pc + z->at(j) * pc * width(w, pc);
      nJets++;
      if (std::fabs(ptw[1][j] - stored->at(j)) > 1e-3 * std::max(1.0, std::fabs((double)ptw[1][j]))) nMismatch++;
    }

    // truth2 recoil matching (analysis.cc:453-496), on pT_calib ranking
    std::vector<int> truthJets;
    for (size_t k = 0; k < tpt->size(); k++) if (tpt->at(k) >= truthJetMin) truthJets.push_back(k);
    int recoilTruth = -1;
    const std::vector<int> c = top3(*calib);
    if (truthJets.size() >= 2 && c[2] != -1) {
      const float rMatch = 0.75*R, rIso = 1.5*R;
      auto dRj = [&](int jet, int k) { float de = eta->at(jet) - teta->at(k), dp = dPhiAbs(phi->at(jet), tphi->at(k)); return std::sqrt(de*de + dp*dp); };
      auto nearest = [&](int jet, int exclude) {
        int best = -1; float bestDR = rMatch;
        for (int k : truthJets) if (k != exclude && dRj(jet, k) < bestDR) { bestDR = dRj(jet, k); best = k; }
        return best;
      };
      int tL = -1;
      const float leadMatch = jtruth->at(c[0]);
      for (int k : truthJets) if (leadMatch > 0 && std::fabs(tpt->at(k) - leadMatch) < 1e-3) tL = k;
      if (tL < 0) tL = nearest(c[0], -1);
      const int tR = (tL >= 0) ? nearest(c[1], tL) : -1;
      bool ok = (tL >= 0 && tR >= 0);
      if (ok) {
        const int t3 = nearest(c[2], -1);
        if (dRj(c[2], tL) < rMatch || (t3 >= 0 && t3 != tR && dRj(c[2], tR) >= rMatch)) ok = false;
        for (int k : truthJets) if (k != tL && k != tR && (dRj(c[1], k) < rIso || dRj(c[2], k) < rIso)) ok = false;
      }
      if (ok) {
        double px = 0, py = 0, pz = 0;
        for (int k : {c[1], c[2]}) {
          const double pt = calib->at(k);
          px += pt*std::cos(phi->at(k)); py += pt*std::sin(phi->at(k)); pz += pt*std::sinh(eta->at(k));
        }
        const double ptSum = std::hypot(px, py);
        if (ptSum > 0) {
          const double de = std::asinh(pz/ptSum) - teta->at(tR), dp = dPhiAbs(std::atan2(py, px), tphi->at(tR));
          if (std::sqrt(de*de + dp*dp) < rMatch) recoilTruth = tR;
        }
      }
    }
    if (recoilTruth >= 0) nRecoilEvt++;

    const double wz = zr->GetBinContent(zr->FindBin(vz));
    for (int cfg = 0; cfg < nCfg; cfg++) {
      const bool recoilSmear = (cfg == 1 && recoilTruth >= 0);
      for (int w = 0; w < nW; w++) {
        const std::vector<float> &cur = ptw[w];
        const std::vector<int> idx = recoilSmear ? c : top3(cur);
        if (idx[2] == -1) continue;
        const float lead = cur[idx[0]];
        float pt2 = recoilSmear ? calib->at(idx[1]) : cur[idx[1]];
        float pt3 = recoilSmear ? calib->at(idx[2]) : cur[idx[2]];
        const float cut2 = pt2, cut3 = pt3;
        if (recoilSmear && w > 0) {
          const float sum = vecSum(pt2, phi->at(idx[1]), pt3, phi->at(idx[2]), false);
          const float ptT = tpt->at(recoilTruth);
          const float ptRef = (w == 1) ? std::min(std::max(ptT, 5.01f), 79.9f) : std::min(ptT, 79.9f); // analysis.cc:536 clamp for the template
          const float smeared = sum + recoilZ * ptT * width(w, ptRef);
          if (sum <= 0 || smeared <= 0) continue;
          pt2 *= smeared / sum; pt3 *= smeared / sum;
        }
        const float recoil = vecSum(pt2, phi->at(idx[1]), pt3, phi->at(idx[2]), false);
        if (lead < leadMin || lead > leadMax || cut2 < jetMin || cut3 < jetMin || recoil < recoilMin) continue;
        if (std::fabs(etadet->at(idx[0])) > 1.1 - R || std::fabs(etadet->at(idx[1])) > 1.1 - R || std::fabs(etadet->at(idx[2])) > 1.1 - R) continue;
        if (dPhiAbs(phi->at(idx[0]), phi->at(idx[2])) < SSLDPHI || dPhiAbs(phi->at(idx[0]), phi->at(idx[1])) < SLDPHI) continue;
        const double wt = wz * fit->Eval(lead);
        for (int k = 0; k < nPt; k++)
          if (lead > ptBins[k] && lead < ptBins[k+1]) hxj[cfg][w][k]->Fill(lead / recoil, wt);
        nSel[cfg][w]++;
      }
    }
  }
  printf("events %lld, truth-window %lld, truth2 recoil-matched %lld; jets %lld, tmpl vs jet_pt_smear_reco mismatches %lld\n",
         n, nTruthOK, nRecoilEvt, nJets, nMismatch);
  for (int cfg = 0; cfg < nCfg; cfg++) {
    printf("%-7s selected (unweighted):", cfgNames[cfg]);
    for (int w = 0; w < nW; w++) printf("  %s %.0f", wNames[w], nSel[cfg][w]);
    printf("\n");
    for (int k = 0; k < nPt; k++) {
      printf("   %.0f-%.0f GeV  <xJ>:", ptBins[k], ptBins[k+1]);
      for (int w = 0; w < nW; w++) printf("  %s %.4f+-%.4f", wNames[w], hxj[cfg][w][k]->GetMean(), hxj[cfg][w][k]->GetMeanError());
      printf("\n");
    }
  }

  TFile fout(Form("jet12_xj_R%02d.root", radius), "RECREATE");
  for (int cfg = 0; cfg < nCfg; cfg++) for (int w = 0; w < nW; w++) for (int k = 0; k < nPt; k++) hxj[cfg][w][k]->Write();
  TParameter<double>("truthLeadMin", truthLeadMin).Write();
  fout.Close();
}
