// Run interpreted: root -l -b -q tcal_runindep.C (output: tcal_runindep.log)
//
// Can a run-INDEPENDENT tower-time correction recover the gain from Dading Chen's per-run
// calibration? delta(calo, E) = median(t_corrected - t_standard) per calorimeter and tower-
// energy bin, derived on the even-indexed condor files (10 runs) and applied to the standard
// time on the odd-indexed files (10 other runs). Compares window in-fractions vs pT and EM
// fraction for: standard, std + delta (run-independent), corrected (full per-run).
// Same tower selection as draw_timingana_tcal.C: E > 0.5 GeV, EMCal ZS excluded.
void tcal_runindep() {
  const char *dir = "/home/samson72/sphnx/gammajet_unfold/claude_checks/jet_timing/timingana/tcal";
  std::vector<TString> files;
  TSystemDirectory sd("d", dir); TList *fl = sd.GetListOfFiles(); fl->Sort();
  for (auto o : *fl) { TString n = o->GetName(); if (n.EndsWith(".root")) files.push_back(TString(dir) + "/" + n); }
  TChain Ttrain("T"), Ttest("T");
  for (size_t i = 0; i < files.size(); i++) (i % 2 ? Ttest : Ttrain).Add(files[i]);
  printf("train files %d, test files %d\n", Ttrain.GetNtrees(), Ttest.GetNtrees());

  const int nte = 20; double teb[nte + 1];
  for (int i = 0; i <= nte; i++) teb[i] = 0.5 * std::pow(40.0, i / (double) nte);
  TH2D *hd[3];
  for (int c = 0; c < 3; c++) hd[c] = new TH2D(Form("hd%d", c), "", nte, teb, 1600, -40, 40);

  struct Br { float vz, mbd, pt, eta, ef; std::vector<int> *calo = 0, *st = 0; std::vector<float> *e = 0, *t = 0, *ts = 0, *ee = 0, *et = 0, *ets = 0; std::vector<bool> *eg = 0; };
  auto attach = [](TChain &T, Br &b) {
    T.SetBranchAddress("vz", &b.vz); T.SetBranchAddress("mbd_time", &b.mbd); T.SetBranchAddress("jet_pt_calib", &b.pt);
    T.SetBranchAddress("jet_eta", &b.eta); T.SetBranchAddress("jet_emfrac", &b.ef);
    T.SetBranchAddress("tw_calo", &b.calo); T.SetBranchAddress("tw_e", &b.e); T.SetBranchAddress("tw_time", &b.t); T.SetBranchAddress("tw_time_std", &b.ts);
    T.SetBranchAddress("em_e", &b.ee); T.SetBranchAddress("em_time", &b.et); T.SetBranchAddress("em_time_std", &b.ets);
    T.SetBranchAddress("em_isgood", &b.eg); T.SetBranchAddress("em_status", &b.st);
  };
  // loop over selected towers: f(calo, E, t_corr, t_std)
  auto loop = [](TChain &T, Br &b, std::function<void(int, float, float, float)> tower, std::function<void(float, float, float)> jetEnd) {
    for (Long64_t i = 0; i < T.GetEntries(); i++) {
      T.GetEntry(i);
      if (!std::isfinite(b.mbd) || !std::isfinite(b.vz) || std::fabs(b.vz) > 60 || std::fabs(b.eta) > 0.7) continue;
      float jpt = b.pt / 0.9265; if (jpt < 5) continue;
      auto add = [&](int c, float e, float tc, float ts) { if (e > 0.5 && std::isfinite(tc) && std::isfinite(ts)) tower(c, e, tc, ts); };
      for (size_t j = 0; j < b.ee->size(); j++) if (b.eg->at(j) && !((b.st->at(j) >> 5) & 1)) add(0, b.ee->at(j), b.et->at(j), b.ets->at(j));
      for (size_t j = 0; j < b.e->size(); j++) if (b.calo->at(j) > 0) add(b.calo->at(j), b.e->at(j), b.t->at(j), b.ts->at(j));
      jetEnd(jpt, b.ef, b.mbd);
    }
  };

  // pass 1: derive delta on training runs
  Br b1; attach(Ttrain, b1);
  loop(Ttrain, b1, [&](int c, float e, float tc, float ts) { hd[c]->Fill(e, tc - ts); }, [](float, float, float) {});
  double delta[3][nte + 2];
  const char *cn[3] = {"EMCal", "iHCal", "oHCal"};
  printf("\ndelta(calo, E) = median(t_corr - t_std) [ns], training runs\n");
  for (int c = 0; c < 3; c++) {
    double fallback = 0;
    { TH1D *p = hd[c]->ProjectionY("all", 1, nte); double q = 0.5, m = 0; if (p->GetEntries() > 0) p->GetQuantiles(1, &m, &q); fallback = m; delete p; }
    printf("  %-6s", cn[c]);
    for (int bx = 0; bx <= nte + 1; bx++) {
      TH1D *p = hd[c]->ProjectionY(Form("p%d%d", c, bx), bx, bx);
      double q = 0.5, m = fallback;
      if (p->GetEntries() >= 30) p->GetQuantiles(1, &m, &q);
      delta[c][bx] = m; delete p;
      if (bx >= 1 && bx <= nte && bx % 4 == 1) printf("  E=%.1f: %.2f", hd[c]->GetXaxis()->GetBinCenter(bx), m);
    }
    printf("\n");
  }

  // pass 2: evaluate on test runs
  const int NV = 3; const char *vn[NV] = {"standard", "std+delta", "corrected"};
  struct J { float pt, ef, d[NV]; }; std::vector<J> jets;
  double tw[NV], te[NV];
  Br b2; attach(Ttest, b2);
  for (int v = 0; v < NV; v++) tw[v] = te[v] = 0;
  loop(Ttest, b2,
       [&](int c, float e, float tc, float ts) {
         int bx = std::min(std::max(hd[c]->GetXaxis()->FindBin(e), 1), nte);
         float t[NV] = {ts, (float) (ts + delta[c][bx]), tc};
         for (int v = 0; v < NV; v++) { tw[v] += e * t[v]; te[v] += e; }
       },
       [&](float jpt, float ef, float mbd) {
         if (te[0] > 0) { J j; j.pt = jpt; j.ef = ef; for (int v = 0; v < NV; v++) j.d[v] = mbd - tw[v] / te[v]; jets.push_back(j); }
         for (int v = 0; v < NV; v++) tw[v] = te[v] = 0;
       });
  printf("\ntest jets %zu\n", jets.size());
  double med[NV];
  for (int v = 0; v < NV; v++) {
    std::vector<double> x; for (auto &j : jets) x.push_back(j.d[v]); std::sort(x.begin(), x.end());
    med[v] = x[x.size() / 2];
    printf("  %-10s median %6.2f  half 16-84 %.2f ns\n", vn[v], med[v], 0.5 * (x[size_t(.84 * x.size())] - x[size_t(.16 * x.size())]));
  }
  double ptb[] = {5, 7, 9, 11, 13, 15, 18, 22, 28, 40}; double efb[] = {0, 0.2, 0.4, 0.6, 0.8, 1.0001};
  for (double hw : {2.0, 3.0}) {
    printf("\nwindow |t_MBD - t_jet - median| < %.0f ns, in-fraction\n   pT:       ", hw);
    for (int b = 0; b < 9; b++) printf(" %5.0f", ptb[b]);
    printf("   | EMfrac:"); for (int b = 0; b < 5; b++) printf(" %4.1f", efb[b]); printf("\n");
    for (int v = 0; v < NV; v++) {
      printf("  %-10s", vn[v]);
      for (int b = 0; b < 9; b++) { int n = 0, k = 0; for (auto &j : jets) if (j.pt >= ptb[b] && j.pt < ptb[b + 1]) { n++; if (std::fabs(j.d[v] - med[v]) < hw) k++; } printf(" %5.3f", n ? double(k) / n : -1); }
      printf("   |       ");
      for (int b = 0; b < 5; b++) { int n = 0, k = 0; for (auto &j : jets) if (j.ef >= efb[b] && j.ef < efb[b + 1]) { n++; if (std::fabs(j.d[v] - med[v]) < hw) k++; } printf(" %4.2f", n ? double(k) / n : -1); }
      printf("\n");
    }
  }
}
