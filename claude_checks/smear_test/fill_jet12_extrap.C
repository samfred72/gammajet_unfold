// Plain ROOT (no libgammajet_unfold): runs on SDCC, where the trees are. Output
// jet12_extrap_R<RR>_<flav>.root is read by draw_jet12_extrap.C.
//
// Jet12 (Pythia) leading / subleading / subsubleading jet pT for four jet-pT definitions:
//   calib   unsmeared pt_calib
//   tmpl    the treemaker's nominal JER smear: w(pt) = h_jer_smear_r04_pileup_EMfracJES_nominal->Interpolate(pt),
//           flat below the first bin centre (5.1875 GeV)
//   lin15   template above 15 GeV; below, its tangent at 15 GeV continued down: w(15) + w'(15) (pt - 15)
//   lin5    template above the first bin centre c1 = 5.1875 GeV; below, the template's first segment
//           (the line through bins 1 and 2, i.e. its slope at ~5 GeV) continued down. There is no slope
//           "at 5 GeV" to take directly, since Interpolate is flat below c1.
// Smearing as in the treemaker (multiJet DijetTreeMaker.cc:1739-1760, smear_pt() at :2373):
// smeared = pt_calib + z * pt_ref * w(pt_ref), with the stored per-jet z (jet_smear_z_R) and pt_ref = matched
// truth pT (jet_truth_pt_R) or pt_calib if unmatched (useTruthRef), or pt_calib (!useTruthRef). The recomputed
// template smear is checked against the stored jet_pt_smear_{truth,reco}_R.
//
// Selection: |z_vtx| < 60 cm and leading truth jet pT above the Jet12 full-efficiency threshold, as
// multijet/analysis.cc; jets with |eta_det| < 1.1 - R, ranked by pT separately in each definition.
//
// Storage floor: the treemaker stores a jet only if raw pT >= 4 GeV or pt_calib or one of its six template
// smears is >= 4 GeV (DijetTreeMaker.cc:1767); the widest of those is the sys-up template. Where a linear
// width exceeds the sys-up width (lin5 below ~4.5 GeV), jets with pt_calib < 4 that it would push above
// 4 GeV can be missing from the tree. Estimated here from the stored ones: a stored jet with pt_calib = c < 4,
// raw pT < 4 and no truth match (pt_ref = c) passed with probability P_up = 1 - Phi(z_up),
// z_up = (4 - c)/(c w_up(c)), so it stands for 1/P_up jets, of which a fraction Phi(z_up) - Phi(z_L)
// (z_L the same with the linear width) would pass the linear smear but are missing. h_missing_<L> holds
// that expected count, placed at the median smeared pT of the missing band (all ranks together - the
// missing jets are not ranked into events). The 15/5 GeV event skim (any flavour) can likewise drop events
// that only a wider linear smear would have kept; not estimated.

double tangentSlope(TH1D *h, double x) { return (h->Interpolate(x + 0.01) - h->Interpolate(x - 0.01)) / 0.02; }

void fill_jet12_extrap(int radius = 4, bool useTruthRef = true,
                       const char *tree = "/sphenix/tg/tg01/jets/samfred/multiJet_full_hadded/multijet_pythia_Jet12.root",
                       Long64_t maxEvents = -1) {
  TH1::SetDefaultSumw2();
  const float R = radius / 10.0f;
  const double truthThresholdJet12[8] = {0, 0, 12, 13, 14, 19, 22, 24}; // multijet/analysis.cc, by radius
  const double truthLeadMin = (radius == 10) ? 25 : truthThresholdJet12[radius];
  const double vzMax = 60;
  const double storeCut = 4; // DijetTreeMaker m_store_pt_cut and pt_cut
  const char *flav = useTruthRef ? "truth" : "reco";

  TFile fjer("jerband_smearing_templates.root");
  TH1D *hNom = (TH1D*)fjer.Get("h_jer_smear_r04_pileup_EMfracJES_nominal");
  TH1D *hUp  = (TH1D*)fjer.Get("h_jer_smear_r04_pileup_EMfracJES_sysup");
  hNom->SetDirectory(0); hUp->SetDirectory(0);
  fjer.Close();

  const double p15 = 15, w15 = hNom->Interpolate(p15), s15 = tangentSlope(hNom, p15);
  const double c1 = hNom->GetBinCenter(1), w5 = hNom->GetBinContent(1);
  const double s5 = (hNom->GetBinContent(2) - hNom->GetBinContent(1)) / (hNom->GetBinCenter(2) - c1);
  printf("lin15: w(15) = %.4f, slope = %.5f /GeV -> w(0) = %.4f, w(5) = %.4f, w(10) = %.4f\n",
         w15, s15, w15 - s15*p15, w15 + s15*(5 - p15), w15 + s15*(10 - p15));
  printf("lin5 : w(%.4f) = %.4f, slope = %.5f /GeV -> w(0) = %.4f, w(2) = %.4f, w(4) = %.4f\n",
         c1, w5, s5, w5 - s5*c1, w5 + s5*(2 - c1), w5 + s5*(4 - c1));
  printf("template: w(10) = %.4f, w(5) = %.4f; sys up w(<5.19) = %.4f\n", hNom->Interpolate(10), hNom->Interpolate(5), hUp->Interpolate(3));

  const int nScen = 4, nRank = 3;
  const char *scenNames[nScen] = {"calib", "tmpl", "lin15", "lin5"};
  const char *rankNames[nRank] = {"leading", "subleading", "subsubleading"};
  auto width = [&](int sc, double pt) -> double {
    if (sc == 1) return hNom->Interpolate(pt);
    if (sc == 2) return pt >= p15 ? hNom->Interpolate(pt) : w15 + s15*(pt - p15);
    if (sc == 3) return pt >= c1 ? hNom->Interpolate(pt) : w5 + s5*(pt - c1);
    return 0;
  };

  TH1D *h[nScen][nRank];
  for (int sc = 0; sc < nScen; sc++)
    for (int r = 0; r < nRank; r++)
      h[sc][r] = new TH1D(Form("h_%s_%s", scenNames[sc], rankNames[r]), ";p_{T}^{jet} [GeV];Counts", 100, 0, 50);
  TH1D *hMiss[nScen] = {nullptr};
  for (int sc = 2; sc < nScen; sc++)
    hMiss[sc] = new TH1D(Form("h_missing_%s", scenNames[sc]), ";p_{T}^{jet} [GeV];Expected missing jets", 100, 0, 50);

  TFile *f = TFile::Open(tree);
  TTree *t = (TTree*)f->Get("ttree");
  std::vector<float> *calib = nullptr, *raw = nullptr, *truthpt = nullptr, *z = nullptr, *etadet = nullptr, *stored = nullptr, *tjet = nullptr;
  float vz = 0;
  std::vector<TString> active = {TString::Format("jet_pt_calib_%d", radius), TString::Format("jet_pt_%d", radius),
                                 TString::Format("jet_truth_pt_%d", radius), TString::Format("jet_smear_z_%d", radius),
                                 TString::Format("jet_eta_det_%d", radius), TString::Format("jet_pt_smear_%s_%d", flav, radius),
                                 TString::Format("truth_jet_pt_%d", radius), TString("mbd_vertex_z")};
  t->SetBranchStatus("*", 0);
  for (const TString &b : active) t->SetBranchStatus(b, 1);
  t->SetBranchAddress(Form("jet_pt_calib_%d", radius), &calib);
  t->SetBranchAddress(Form("jet_pt_%d", radius), &raw);
  t->SetBranchAddress(Form("jet_truth_pt_%d", radius), &truthpt);
  t->SetBranchAddress(Form("jet_smear_z_%d", radius), &z);
  t->SetBranchAddress(Form("jet_eta_det_%d", radius), &etadet);
  t->SetBranchAddress(Form("jet_pt_smear_%s_%d", flav, radius), &stored);
  t->SetBranchAddress(Form("truth_jet_pt_%d", radius), &tjet);
  t->SetBranchAddress("mbd_vertex_z", &vz);

  const Long64_t n = (maxEvents > 0) ? std::min(maxEvents, t->GetEntries()) : t->GetEntries();
  Long64_t nPass = 0, nJetsChecked = 0, nMismatch = 0;
  std::vector<float> pts[nScen];
  // missing-jet estimate broken down by pt_calib (1 GeV slices): stored jets used, expected missing
  double missStored[nScen][4] = {{0}}, missSum[nScen][4] = {{0}};
  for (Long64_t e = 0; e < n; e++) {
    t->GetEntry(e);
    if (std::fabs(vz) > vzMax) continue;
    if (tjet->empty() || *std::max_element(tjet->begin(), tjet->end()) < truthLeadMin) continue;
    nPass++;
    for (int sc = 0; sc < nScen; sc++) pts[sc].clear();
    for (size_t j = 0; j < calib->size(); j++) {
      const double pc = calib->at(j);
      const bool matched = useTruthRef && truthpt->at(j) > 0;
      const double ref = matched ? truthpt->at(j) : pc;
      const double zz = z->at(j);
      nJetsChecked++;
      const double tm = pc + zz * ref * width(1, ref);
      if (std::fabs(tm - stored->at(j)) > 1e-3 * std::max(1.0, std::fabs(tm))) nMismatch++;

      if (std::fabs(etadet->at(j)) > 1.1 - R) continue;
      pts[0].push_back(pc);
      for (int sc = 1; sc < nScen; sc++) pts[sc].push_back(pc + zz * ref * width(sc, ref));

      // missing-jet estimate (see header)
      if (pc > 0 && pc < storeCut && raw->at(j) < storeCut && !matched) {
        const double zUp = (storeCut - pc) / (pc * hUp->Interpolate(pc));
        const double pUp = 1 - ROOT::Math::normal_cdf(zUp);
        if (pUp <= 0) continue;
        for (int sc = 2; sc < nScen; sc++) {
          const double zL = (storeCut - pc) / (pc * width(sc, pc));
          if (zL >= zUp) continue;
          const double fL = ROOT::Math::normal_cdf(zL), fU = ROOT::Math::normal_cdf(zUp);
          const double zMid = ROOT::Math::normal_quantile(0.5 * (fL + fU), 1);
          hMiss[sc]->Fill(pc + zMid * pc * width(sc, pc), (fU - fL) / pUp);
          missStored[sc][(int)pc]++; missSum[sc][(int)pc] += (fU - fL) / pUp;
        }
      }
    }
    for (int sc = 0; sc < nScen; sc++) {
      std::sort(pts[sc].begin(), pts[sc].end(), std::greater<float>());
      for (int r = 0; r < nRank && r < (int)pts[sc].size(); r++) h[sc][r]->Fill(pts[sc][r]);
    }
  }
  printf("events %lld, passing %lld; jets checked %lld, template-recompute mismatches vs jet_pt_smear_%s_%d: %lld\n",
         n, nPass, nJetsChecked, flav, radius, nMismatch);
  for (int sc = 2; sc < nScen; sc++) {
    const double sub = h[sc][2]->Integral(h[sc][2]->FindBin(4.001), h[sc][2]->FindBin(5.999));
    printf("%-6s expected missing jets: %.0f total, %.0f at 4-6 GeV (vs %.0f subsubleading jets at 4-6 GeV)\n",
           scenNames[sc], hMiss[sc]->Integral(), hMiss[sc]->Integral(hMiss[sc]->FindBin(4.001), hMiss[sc]->FindBin(5.999)), sub);
    for (int k = 0; k < 4; k++)
      printf("   pt_calib %d-%d GeV: %6.0f stored jets -> %8.0f expected missing\n", k, k + 1, missStored[sc][k], missSum[sc][k]);
  }

  TFile fout(Form("jet12_extrap_R%02d_%s.root", radius, flav), "RECREATE");
  for (int sc = 0; sc < nScen; sc++) for (int r = 0; r < nRank; r++) h[sc][r]->Write();
  for (int sc = 2; sc < nScen; sc++) hMiss[sc]->Write();
  TParameter<double>("truthLeadMin", truthLeadMin).Write();
  fout.Close();
}
