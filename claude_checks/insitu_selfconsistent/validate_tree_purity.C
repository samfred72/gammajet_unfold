// Validation for a self-consistent in-situ scan: can the in-situ trees (insitu/inputs/Data_<systag>_insitu.root)
// replace the unfold pass as the source of the Data ABCD counts that puritymaker.C uses?
// For each systag and radius, at the p_a the unfold pass used (the pass-2 table of the Oct 8 run,
// logs/full_pipeline_20261008_162751/jes_used_by_pass2.txt, hardcoded below):
//   1. rebuild hclusterpt_abcd<ir>_<j> from the tree with unfolder.cc's pairing (jet pT = raw/p_a > 5 GeV,
//      x_J floor at floorScale 1; eta/dphi/photon cuts are already applied when the tree is filled),
//      and compare bin by bin with hists/Data_<systag>_unfolding.root;
//   2. run puritymaker.C's combine_hists on the rebuilt counts and compare with hists/purity_<systag>.root.
#include "../../macros/puritymaker.C"

// Pass-2 table (jes_used_by_pass2.txt). jes_high/low follow unfolder.cc: nominal -/+ stat.
const map<string, vector<float>> paTable = {
  {"nominal",      {0.9034, 0.9033, 0.9158, 0.9116, 0.9172, 0.9151, 0.9289}},
  {"JERhigh",      {0.8946, 0.9001, 0.9072, 0.9060, 0.9131, 0.9114, 0.9276}},
  {"JERlow",       {0.9145, 0.9180, 0.9203, 0.9141, 0.9197, 0.9169, 0.9304}},
  {"emscale_high", {0.9097, 0.9099, 0.9203, 0.9141, 0.9217, 0.9194, 0.9322}},
  {"emscale_low",  {0.9006, 0.9033, 0.9101, 0.9064, 0.9124, 0.9092, 0.9250}},
  {"threejet",     {0.8946, 0.9033, 0.9149, 0.9211, 0.9170, 0.9142, 0.9377}},
  {"narrowBDT",    {0.9013, 0.9093, 0.9229, 0.9194, 0.9188, 0.9152, 0.9300}},
  {"narrowISO",    {0.8973, 0.9026, 0.9101, 0.9059, 0.9139, 0.9150, 0.9235}},
  {"EMRhigh",      {0.9145, 0.9184, 0.9265, 0.9295, 0.9305, 0.9285, 0.9441}},
  {"EMRlow",       {0.9034, 0.9033, 0.9142, 0.9116, 0.9172, 0.9142, 0.9285}},
  {"narrowBDTbkg", {0.9097, 0.9096, 0.9187, 0.9186, 0.9218, 0.9225, 0.9322}},
  {"narrowISObkg", {0.9034, 0.9033, 0.9142, 0.9114, 0.9172, 0.9140, 0.9285}},
  {"wideISObkg",   {0.9034, 0.9038, 0.9158, 0.9116, 0.9181, 0.9151, 0.9288}},
};
const float nomPa[7]   = {0.9034, 0.9033, 0.9158, 0.9116, 0.9172, 0.9151, 0.9289};
const float statLo[7]  = {0.0056, 0.0041, 0.0083, 0.0070, 0.0065, 0.0072, 0.0079};
const float statHi[7]  = {0.0112, 0.0103, 0.0067, 0.0041, 0.0063, 0.0074, 0.0085};

float paFor(const string & systag, int ir) {
  if (systag == "jes_high") return nomPa[ir] - statLo[ir];
  if (systag == "jes_low")  return nomPa[ir] + statHi[ir];
  return paTable.at(systag)[ir];
}

void validate_tree_purity() {
  TH1::AddDirectory(kFALSE);
  vector<string> systags = {"nominal","JERhigh","JERlow","emscale_high","emscale_low","jes_high","jes_low",
                            "threejet","narrowBDT","narrowISO","EMRhigh","EMRlow","narrowBDTbkg","narrowISObkg","wideISObkg"};
  TFile * fscratch = TFile::Open("validate_tree_purity_scratch.root", "RECREATE");
  long totBins = 0, badBins = 0; double maxCountDiff = 0, maxPurDiff = 0, maxPurErrDiff = 0;
  for (const string & systag : systags) {
    TFile * fu = TFile::Open(Form("%s/hists/Data_%s_unfolding.root", ana::dir(), systag.c_str()));
    TFile * ft = TFile::Open(Form("%s/insitu/inputs/Data_%s_insitu.root", ana::dir(), systag.c_str()));
    TFile * fp = TFile::Open(ana::purityFilename(systag).c_str());
    TTree * t = (TTree*)ft->Get("insitutree");
    float phoPt, jetPt, w; int abcd, irT;
    t->SetBranchAddress("pho_pt", &phoPt); t->SetBranchAddress("jet_pt", &jetPt);
    t->SetBranchAddress("abcd", &abcd); t->SetBranchAddress("weight", &w); t->SetBranchAddress("ir", &irT);
    vector<array<TH1D*,4>> hT(ana::nJetR);
    for (int ir = 0; ir < ana::nJetR; ir++) for (int j = 0; j < 4; j++)
      hT[ir][j] = new TH1D(Form("t_%s_%d_%d", systag.c_str(), ir, j), "", ana::nPtBins, ana::ptBins);
    for (Long64_t e = 0; e < t->GetEntries(); e++) {
      t->GetEntry(e);
      float pa = paFor(systag, irT);
      float recoJetPt = jetPt / pa;                          // unfolder.cc: jet_pt_calib / p_a (float)
      if (recoJetPt <= ana::jet_calib_pt_cut[irT]) continue; // ispaired pre-filter
      int ptbin = ana::findPtBin(phoPt);
      if (ptbin == -1) continue;
      float lowval = ana::jet_calib_pt_cut[irT]/ana::ptBins[ptbin];
      float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(lowval)+1];
      if (recoJetPt/phoPt < lowbin) continue;                // check_pair floor, floorScale 1
      hT[irT][abcd]->Fill(phoPt, w);
    }
    int sysBad = 0; double sysMaxPur = 0;
    for (int ir = 0; ir < ana::nJetR; ir++) {
      TH1D * h[4]; TH1D * fpx[4];
      for (int j = 0; j < 4; j++) {
        TH1D * hu = (TH1D*)fu->Get(Form("hclusterpt_abcd%d_%d", ir, j));
        for (int b = 1; b <= ana::nPtBins; b++) {
          totBins++;
          double d = fabs(hu->GetBinContent(b) - hT[ir][j]->GetBinContent(b));
          maxCountDiff = max(maxCountDiff, d);
          if (d > 0) { badBins++; sysBad++;
            printf("  MISMATCH %s R%d region %d ptbin %d: unfold %.0f tree %.0f\n", systag.c_str(), ir, j, b,
                   hu->GetBinContent(b), hT[ir][j]->GetBinContent(b)); }
        }
        h[j] = hT[ir][j];
      }
      // Leakage fractions exactly as puritymaker.C (MC truth-matched, same drawer).
      drawer d("pythia", systag);
      TH1D * hp[4];
      for (int j = 0; j < 4; j++) hp[j] = d.get(Form("hclusterpt_abcd_truthmatched%i_%i", ir, j), 1);
      for (int j = 0; j < 4; j++) { fpx[j] = (TH1D*)hp[j]->Clone(Form("fp%i", j)); fpx[j]->Divide(hp[j], hp[0]); }
      fscratch->cd(); fscratch->mkdir(Form("%s_%s", systag.c_str(), ana::rnames[ir]))->cd();
      TGraphAsymmErrors * gT = combine_hists(h, fpx);
      TGraphAsymmErrors * gP = (TGraphAsymmErrors*)fp->Get(Form("%s/combined", ana::rnames[ir]));
      for (int i = 0; i < gP->GetN(); i++) {
        double dv = fabs(gT->GetPointY(i) - gP->GetPointY(i));
        double de = max(fabs(gT->GetErrorYlow(i) - gP->GetErrorYlow(i)), fabs(gT->GetErrorYhigh(i) - gP->GetErrorYhigh(i)));
        maxPurDiff = max(maxPurDiff, dv); maxPurErrDiff = max(maxPurErrDiff, de); sysMaxPur = max(sysMaxPur, dv);
      }
    }
    printf("%-13s entries %7lld  count mismatches %3d  max |dP_A| %.2e\n", systag.c_str(), t->GetEntries(), sysBad, sysMaxPur);
    fu->Close(); ft->Close(); fp->Close();
  }
  printf("\nTOTAL: %ld/%ld (systag, R, region, pT) count bins differ (max |diff| %.3g); max |dP_A| %.2e, max |d err| %.2e\n",
         badBins, totBins, maxCountDiff, maxPurDiff, maxPurErrDiff);
  fscratch->Close();
}
