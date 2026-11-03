// Compares two unfolding outputs (every TH1, recursively: bin contents and errors) and two in-situ
// trees (every float leaf of every tree, entry by entry). Prints the number of differing objects.
int nH = 0, nHdiff = 0;
double maxRel = 0;

void compareDir(TDirectory * a, TDirectory * b) {
  for (TObject * k : *a->GetListOfKeys()) {
    TObject * oa = ((TKey*)k)->ReadObj();
    TObject * ob = b->Get(k->GetName());
    if (oa->InheritsFrom(TDirectory::Class())) { compareDir((TDirectory*)oa, (TDirectory*)ob); continue; }
    if (!oa->InheritsFrom(TH1::Class())) continue;
    TH1 * ha = (TH1*)oa, * hb = (TH1*)ob;
    nH++;
    if (!hb || ha->GetNcells() != hb->GetNcells()) { nHdiff++; continue; }
    bool diff = false;
    for (int i = 0; i < ha->GetNcells(); i++) {
      double ca = ha->GetBinContent(i), cb = hb->GetBinContent(i);
      if (ca != cb || ha->GetBinError(i) != hb->GetBinError(i)) {
        diff = true;
        if (ca != 0) maxRel = max(maxRel, fabs(cb/ca - 1));
      }
    }
    if (diff) nHdiff++;
  }
}

void compareTrees(TDirectory * a, TDirectory * b, long & nTrees, long & nDiffEntries, long & nEntries) {
  for (TObject * k : *a->GetListOfKeys()) {
    TObject * oa = ((TKey*)k)->ReadObj();
    if (oa->InheritsFrom(TDirectory::Class())) { compareTrees((TDirectory*)oa, b->GetDirectory(k->GetName()), nTrees, nDiffEntries, nEntries); continue; }
    if (!oa->InheritsFrom(TTree::Class())) continue;
    TTree * ta = (TTree*)oa, * tb = (TTree*)b->Get(k->GetName());
    nTrees++;
    if (!tb || ta->GetEntries() != tb->GetEntries()) { printf("tree %s: entries %lld vs %lld\n", k->GetName(), ta->GetEntries(), tb ? tb->GetEntries() : -1); nDiffEntries++; continue; }
    vector<string> names;
    for (TObject * l : *ta->GetListOfLeaves()) if (string(((TLeaf*)l)->GetTypeName()) == "Float_t" && ((TLeaf*)l)->GetLen() == 1) names.push_back(l->GetName());
    vector<float> va(names.size()), vb(names.size());
    for (size_t i = 0; i < names.size(); i++) { ta->SetBranchAddress(names[i].c_str(), &va[i]); tb->SetBranchAddress(names[i].c_str(), &vb[i]); }
    for (Long64_t e = 0; e < ta->GetEntries(); e++) {
      ta->GetEntry(e); tb->GetEntry(e); nEntries++;
      if (va != vb) nDiffEntries++;
    }
  }
}

void compare_files(const char * ha, const char * hb, const char * ia, const char * ib) {
  TFile fa(ha), fb(hb);
  compareDir(&fa, &fb);
  printf("histograms: %d compared, %d differ (max relative bin difference %.2e)\n", nH, nHdiff, maxRel);
  TFile ga(ia), gb(ib);
  long nTrees = 0, nDiff = 0, nEnt = 0;
  compareTrees(&ga, &gb, nTrees, nDiff, nEnt);
  printf("in-situ trees: %ld trees, %ld entries, %ld entries differ\n", nTrees, nEnt, nDiff);
}
