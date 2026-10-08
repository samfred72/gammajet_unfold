// Plain ROOT (no libgammajet_unfold): runs on SDCC, where the trees are. Output
// jet12_smear_R<RR>_<flav>.root is read by draw_jet12_smear.C.

// Jet12 (Pythia) leading / subleading / subsubleading jet pT with the treemaker's JER smearing
// vs. a version whose width is linearly extrapolated below 10 GeV.
//
// Treemaker smearing (multiJet DijetTreeMaker.cc:1739-1760, smear_pt() at :2373; same scheme as
// gammajet CaloAna.cc): smeared = pt_calib + z * pt_ref * w(pt_ref), one standard-normal z per jet
// shared by nominal/up/down, w = h_jer_smear_r04_pileup_EMfracJES_*->Interpolate(pt_ref). The
// template starts at 5 GeV and Interpolate holds w flat below 5.19 GeV (see draw_smear_function.C).
// The tree stores pt_calib, z (jet_smear_z_R) and the matched truth pT (jet_truth_pt_R, -1 if
// unmatched), so each jet is re-smeared here with the same z and only w changed:
//
//   template  w(pt) = Interpolate(pt)                                (as in the treemaker)
//   linear    w(pt) = Interpolate(pt) for pt >= 10 GeV,
//             w(10) + w'(10) * (pt - 10) below, w'(10) the template's slope at 10 GeV
//
// separately for nominal, sys up and sys down. useTruthRef picks the "_truth" flavour (pt_ref = matched
// truth pT, else pt_calib; what CaloAna uses for MC) or "_reco" (pt_ref = pt_calib). The recomputed
// template smear is checked against the stored jet_pt_smear_{truth,reco}_R branch.
//
// Selection: |z_vtx| < 60 cm and leading truth jet pT above the Jet12 full-efficiency threshold, as
// multijet/analysis.cc; jets with |eta_det| < 1.1 - R, ranked by pT separately in each scenario.
// Storage floor: the treemaker only stores a jet if raw pT, pt_calib or one of its six template smears
// reaches 4 GeV (DijetTreeMaker.cc:1767). With the shared z, a linear smear can only pass where a
// template smear did not if its width is larger than the template sys-up width: the linear nominal and
// down never are, and linear up only below pT_ref ~ 0.9 GeV (not reachable to 4 GeV). So no stored-jet
// loss affects these spectra.

void fill_jet12_smear(int radius = 4, bool useTruthRef = true,
                      const char *tree = "/sphenix/tg/tg01/jets/samfred/multiJet_full_hadded/multijet_pythia_Jet12.root",
                      Long64_t maxEvents = -1) {
  TH1::SetDefaultSumw2();

  const float R = radius / 10.0f;
  const double truthThresholdJet12[8] = {0, 0, 12, 13, 14, 19, 22, 24}; // multijet/analysis.cc, by radius
  const double truthLeadMin = (radius == 10) ? 25 : truthThresholdJet12[radius];
  const double vzMax = 60;
  const double linPivot = 10; // GeV: below this the width is linearly extrapolated
  const char *flav = useTruthRef ? "truth" : "reco";

  // width functions
  TFile fjer("jerband_smearing_templates.root");
  const char *jerNames[3] = {"nominal", "sysup", "sysdown"};
  TH1D *hW[3];
  double w10[3], slope[3];
  for (int s = 0; s < 3; s++) {
    hW[s] = (TH1D*)fjer.Get(Form("h_jer_smear_r04_pileup_EMfracJES_%s", jerNames[s]));
    hW[s]->SetDirectory(0);
    w10[s] = hW[s]->Interpolate(linPivot);
    slope[s] = (hW[s]->Interpolate(linPivot + 0.01) - hW[s]->Interpolate(linPivot - 0.01)) / 0.02;
    printf("%-8s w(10) = %.4f  slope = %.5f /GeV  -> w_lin(0) = %.4f, w_lin(5) = %.4f (template %.4f)\n",
           jerNames[s], w10[s], slope[s], w10[s] - slope[s]*linPivot, w10[s] + slope[s]*(5 - linPivot),
           hW[s]->Interpolate(5));
  }
  fjer.Close();
  auto wTmpl = [&](int s, double pt) { return hW[s]->Interpolate(pt); };
  auto wLin  = [&](int s, double pt) { return pt >= linPivot ? hW[s]->Interpolate(pt) : w10[s] + slope[s]*(pt - linPivot); };

  // scenarios: 0 unsmeared, 1-3 template nom/up/down, 4-6 linear nom/up/down
  const int nScen = 7, nRank = 3;
  const char *scenNames[nScen] = {"calib", "tmpl_nom", "tmpl_up", "tmpl_down", "lin_nom", "lin_up", "lin_down"};
  const char *rankNames[nRank] = {"leading", "subleading", "subsubleading"};
  TH1D *h[nScen][nRank];
  for (int sc = 0; sc < nScen; sc++)
    for (int r = 0; r < nRank; r++)
      h[sc][r] = new TH1D(Form("h_%s_%s", scenNames[sc], rankNames[r]), ";p_{T}^{jet} [GeV];Counts", 100, 0, 50);

  TFile *f = TFile::Open(tree);
  TTree *t = (TTree*)f->Get("ttree");
  std::vector<float> *calib = nullptr, *truthpt = nullptr, *z = nullptr, *etadet = nullptr, *stored = nullptr, *tjet = nullptr;
  float vz = 0;
  t->SetBranchStatus("*", 0);
  std::vector<TString> active = {TString::Format("jet_pt_calib_%d", radius), TString::Format("jet_truth_pt_%d", radius),
                                 TString::Format("jet_smear_z_%d", radius), TString::Format("jet_eta_det_%d", radius),
                                 TString::Format("jet_pt_smear_%s_%d", flav, radius), TString::Format("truth_jet_pt_%d", radius),
                                 TString("mbd_vertex_z")};
  for (const TString &b : active)
    t->SetBranchStatus(b, 1);
  t->SetBranchAddress(Form("jet_pt_calib_%d", radius), &calib);
  t->SetBranchAddress(Form("jet_truth_pt_%d", radius), &truthpt);
  t->SetBranchAddress(Form("jet_smear_z_%d", radius), &z);
  t->SetBranchAddress(Form("jet_eta_det_%d", radius), &etadet);
  t->SetBranchAddress(Form("jet_pt_smear_%s_%d", flav, radius), &stored);
  t->SetBranchAddress(Form("truth_jet_pt_%d", radius), &tjet);
  t->SetBranchAddress("mbd_vertex_z", &vz);

  const Long64_t n = (maxEvents > 0) ? std::min(maxEvents, t->GetEntries()) : t->GetEntries();
  Long64_t nPass = 0, nJetsChecked = 0, nMismatch = 0;
  std::vector<float> pts[nScen];
  for (Long64_t e = 0; e < n; e++) {
    t->GetEntry(e);
    if (std::fabs(vz) > vzMax) continue;
    if (tjet->empty() || *std::max_element(tjet->begin(), tjet->end()) < truthLeadMin) continue;
    nPass++;
    for (int sc = 0; sc < nScen; sc++) pts[sc].clear();
    for (size_t j = 0; j < calib->size(); j++) {
      const double pc = calib->at(j);
      const double ref = (useTruthRef && truthpt->at(j) > 0) ? truthpt->at(j) : pc;
      const double zz = z->at(j);
      // validation against the treemaker's own smear (all jets, before the eta cut)
      const double tm = pc + zz * ref * wTmpl(0, ref);
      nJetsChecked++;
      if (std::fabs(tm - stored->at(j)) > 1e-3 * std::max(1.0, std::fabs(tm))) nMismatch++;

      if (std::fabs(etadet->at(j)) > 1.1 - R) continue;
      pts[0].push_back(pc);
      for (int s = 0; s < 3; s++) {
        pts[1+s].push_back(pc + zz * ref * wTmpl(s, ref));
        pts[4+s].push_back(pc + zz * ref * wLin(s, ref));
      }
    }
    for (int sc = 0; sc < nScen; sc++) {
      std::sort(pts[sc].begin(), pts[sc].end(), std::greater<float>());
      for (int r = 0; r < nRank && r < (int)pts[sc].size(); r++) h[sc][r]->Fill(pts[sc][r]);
    }
  }
  printf("events %lld, passing %lld; jets checked %lld, template-recompute mismatches vs jet_pt_smear_%s_%d: %lld\n",
         n, nPass, nJetsChecked, flav, radius, nMismatch);

  TFile fout(Form("jet12_smear_R%02d_%s.root", radius, flav), "RECREATE");
  for (int sc = 0; sc < nScen; sc++) for (int r = 0; r < nRank; r++) h[sc][r]->Write();
  TParameter<double>("truthLeadMin", truthLeadMin).Write();
  TParameter<double>("nPass", nPass).Write();
  TParameter<double>("nMismatch", nMismatch).Write();
  fout.Close();

}
