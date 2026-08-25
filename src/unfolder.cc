#include "unfolder.h"
#include "insitu_utility.h"
using namespace std;

unfolder::~unfolder() {}

bool unfolder::check_pair(jet_object jet, int ir, pho_object pho, bool isreco, float testPt, float floorScale) {
  float dphi = jet.deltaPhi(pho);
  int iabcd = ana::findabcdBin(pho.iso4, pho.bdt, 0);

  int ptbin = ana::findPtBin(pho.pt);
  // testPt lets a caller test the xJ floor against a pt other than jet.pt (e.g. the
  // insitu tree's raw, uncorrected jet pt - see the in-situ test tree fill below);
  // floorScale loosens/tightens the floor itself by a multiplicative factor around
  // that same jet.pt-based lowbin edge.
  float val = (testPt >= 0 ? testPt : jet.pt)/pho.pt;
  float lowval = ana::jet_calib_pt_cut[ir]/ana::ptBins[ptbin];
  float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(lowval)+1];

  if (ptbin == -1) return false;
  if (val < lowbin*floorScale) return false;
  //if (iabcd != 0) return false;
  if (fabs(pho.eta) > ana::etacut) return false;
  if (fabs(jet.eta) > ana::etacut - ana::JetRs[ir]) return false;
  if (dphi < ana::oppcut) return false;

  return true;
}

bool unfolder::check_match(pho_object p1, pho_object p2, jet_object j1, jet_object j2) {
  if (p1.deltaR(p2) < 0.1 && j1.deltaR(j2) < 0.3) return true;
  else return false;
}
bool unfolder::check_match(pho_object p1, pho_object p2) {
  if (p1.deltaR(p2) < 0.1) return true;
  else return false;
}
bool unfolder::check_match(jet_object j1, jet_object j2) {
  if (j1.deltaR(j2) < 0.3) return true;
  else return false;
}
void unfolder::fill_matrix() {

  nentries = t->GetEntriesFast();
  int nsys = (int)systags.size();
  cout << "running (systags=";
  for (int isys = 0; isys < nsys; isys++) cout << (isys?",":"") << systags[isys];
  cout << ")..." << endl;

  // Per-systag constants - each depends only on that systag's name string, not on the
  // event, so computed once per systag here rather than re-derived every event. See the
  // constructor comment in unfolder.h for the full list/rationale of what each one means.
  vector<int> systagAbcdBinArr(nsys);
  vector<bool> systagThreejetVetoArr(nsys);
  vector<float> systagEmscaleShiftArr(nsys);
  vector<float> jesCorrectionArr(nsys);
  for (int isys = 0; isys < nsys; isys++) {
    const string & systag = systags[isys];
    systagAbcdBinArr[isys] = (systag == "narrowBDT") ? 1 : (systag == "narrowISO") ? 2 : 0;
    systagThreejetVetoArr[isys] = (systag == "threejet");
    // emscale_high/emscale_low shift the EM-calorimeter energy scale by +-1.1%: applied in
    // full to the (entirely-EM) photon cluster, and to only the EM-fraction portion of the
    // jet (the non-EM/hadronic portion of jet_pt_smear_truth[ir] is left untouched). MC-only, same
    // convention as JERhigh/JERlow above - jet_pt_smear_truth[ir] has no Data equivalent, so Data
    // always falls back to its nominal reco pT regardless of systag.
    systagEmscaleShiftArr[isys] = (systag == "emscale_high") ? 0.011 :
                                   (systag == "emscale_low")  ? -0.011 : 0.0;
    // JES: Data's reconstructed jet pT has a residual calibration gap relative to MC (found
    // via the in-situ jet-photon pT-balance study - see insitu/), corrected by dividing by
    // 0.9446 (grid_insitu.C's purity-corrected best-fit p_a, nominal systag - see
    // insitu/output/grid_insitu_nominal.root's "R04/results" tree). jes_high/jes_low vary that
    // correction factor by +-0.03; this +-0.03 has NOT been re-derived from the
    // purity-corrected scan's own (smaller, asymmetric) error - +0.008/-0.006 per the
    // FINAL RESULT printout - so it's carried over unchanged from the old 0.9729 central
    // value pending an explicit decision to update it. Data-only - this corrects a gap
    // specific to real Data's calibration, so MC (already on-scale by construction) is
    // untouched regardless of systag, the reverse convention from
    // JERhigh/JERlow/emscale_high/emscale_low above. Applied to recoJetPt (and everything
    // downstream: unfold response matrices, purity histograms, pairing) below, but
    // deliberately NOT to insitu_jet_pt - see rawJetPt.
    jesCorrectionArr[isys] = 0.9446 + ((systag == "jes_high") ? -0.03 :
                                       (systag == "jes_low")  ?  0.03 : 0.0);
  }

  TCanvas * c = new TCanvas("c","",500,1000);
  gStyle->SetOptStat(0);
  if (dodraw) c->SaveAs(Form("/home/samson72/sphnx/gammajet_unfold/pdfs/event_displays_%s.pdf[",trigger.c_str()));
  int ndraw = 0;

  for (Long64_t e = 0; e < nentries; e++) {
    t->GetEntry(e);
    bool use_half = rand.Integer(2) % 2;
    if (e % 1000 == 0)
      std::cout << "entry " << e << "/" << nentries
        << " (" << (float)e/nentries*100. << "%)\t\r" << std::flush;

    if (fabs(vz) > ana::vzcut) continue;

    // -----------------------
    // Event selection - shared across every systag (check_keep_MC never depends on
    // systag), computed once per event.
    // -----------------------

    vector<bool> keepMC = check_keep_MC(truth_cluster_pt, cluster_pt, truth_jet_pt, jet_pt_smear_truth, trigger);
    if (isMC && !keepMC.at(keepMC.size()-1)) continue;

    // Data/MC vz and cluster-pT reweighting (src/reweight_utility.h) - one weight per
    // event, from the raw (un-systag-shifted) vz/cluster_pt, applied identically across
    // every systag reprocessing below: it corrects an orthogonal Data-vs-MC mismatch
    // (vertex profile, pT-threshold stitching), not something a JES/JER/emscale
    // systematic shift should itself perturb. Data always gets weight 1 (the correction
    // reweights MC to match Data, not the other way around).
    float mcWeight = isMC ? rw.GetWeight(vz, cluster_pt) : 1.0f;

    // -----------------------
    // Per-systag reprocessing of this same event - one full pass over the tree fills
    // every systag's histograms, instead of re-reading the tree once per systag. Every
    // quantity below that varies with systag (recoClusterPt, recoJetPt, ispaired,
    // iabcd_reco) is a closed-form shift or a selection among branches already read by
    // the single t->GetEntry(e) above - see unfolder.h's constructor comment.
    // -----------------------
    for (int isys = 0; isys < nsys; isys++) {
      const string & systag = systags[isys];
      int systagAbcdBin = systagAbcdBinArr[isys];
      bool systagThreejetVeto = systagThreejetVetoArr[isys];
      float systagEmscaleShift = systagEmscaleShiftArr[isys];
      float jesCorrection = jesCorrectionArr[isys];

      // -----------------------
      // Leading photon & isolation
      // -----------------------
      float recoClusterPt = isMC ? cluster_pt * (1.0 + systagEmscaleShift) : cluster_pt;
      pho_object maxpho = pho_object(
          recoClusterPt,
          cluster_e,
          cluster_eta,
          cluster_phi,
          cluster_showershape[10],
          cluster_showershape[11],
          cluster_time,
          cluster_bdt_scores[9],
          pho_object::get_showershape(cluster_showershape, recoClusterPt)
      );
      pho_object maxpho_truth = (isMC ? pho_object(
          truth_cluster_pt,
          truth_cluster_e,
          truth_cluster_eta,
          truth_cluster_phi,
          //cluster_showershape[8], // TEMPORARY!!!!!
          //cluster_showershape[9],
          truth_cluster_iso3,
          truth_cluster_iso4,
          0, // no time object for truth
          0.99, // truth photon is guaranteed a photon
          2 // truth photon is guaranteed a photon
      ) : maxpho);


      vector<jet_object> maxjet(ana::nJetR);
      vector<jet_object> maxjet_truth(ana::nJetR);
      vector<bool> ispaired(ana::nJetR, false);
      vector<bool> ispaired_truth(ana::nJetR, false);
      // Looser than ispaired: tests the xJ floor against rawJetPt (not jesCorrection-
      // boosted recoJetPt) at floorScale=0.95, so the in-situ tree keeps every event a
      // grid_insitu.C scan point down to pa=0.95 could still scale above the floor -
      // see the in-situ test tree fill below and check_pair's testPt/floorScale params.
      vector<bool> ispairedInsitu(ana::nJetR, false);

      for (int ir = 0; ir < ana::nJetR; ir++) {
        if (isMC && !keepMC[ir]) continue;

        // JERhigh/JERlow/emscale_high/emscale_low only mean anything for MC (smearing is
        // applied to MC to match Data's resolution - Data has no smeared-high/low variant of
        // itself); jes_high/jes_low go the other way (Data-only - MC is already on-scale by
        // construction, so jesCorrection never enters the isMC branch below).
        float recoJetPt;
        if (isMC) {
          recoJetPt = (systag == "JERhigh") ? jet_pt_smear_high_truth[ir] :
                      (systag == "JERlow")  ? jet_pt_smear_low_truth[ir]  :
                      // scale only the EM-fraction portion of the jet by (1 + shift); the
                      // non-EM/hadronic portion, jet_pt_smear_truth[ir]*(1-jet_emfrac[ir]), is
                      // untouched.
                      (systagEmscaleShift != 0.0) ?
                          jet_pt_smear_truth[ir]*jet_emfrac[ir]
                            + jet_pt_smear_truth[ir]*jet_emfrac[ir]*systagEmscaleShift
                            + jet_pt_smear_truth[ir]*(1-jet_emfrac[ir]) :
                      jet_pt_smear_truth[ir];
        } else {
          recoJetPt = jet_pt_calib[ir] / jesCorrection;
        }
        // Uncorrected Data jet pT, for the insitu tree only (see insitu_jet_pt fill
        // below) - the insitu study is what jesCorrection (0.9446, see above) is itself
        // derived from, so baking that correction into insitu_jet_pt would make the
        // in-situ grid scan measure only the residual gap around an already-applied
        // guess instead of the actual Data/MC JES gap. MC never has jesCorrection
        // applied in the first place (see the isMC branch above), so this is just
        // recoJetPt there.
        float rawJetPt = isMC ? recoJetPt : jet_pt_calib[ir];
        maxjet[ir] = jet_object(
            recoJetPt,
            jet_e[ir],
            jet_eta[ir],
            jet_phi[ir],
            jet_emfrac[ir],
            0,
            0,
            jet_time[ir]
        );
        maxjet_truth[ir] = (isMC ? jet_object(
            truth_jet_pt[ir],
            truth_jet_e[ir],
            truth_jet_eta[ir],
            truth_jet_phi[ir],
            0, 0, 0, 0) : maxjet[ir]);

        hphodr[isys][ir]->Fill(maxpho.deltaR(maxpho_truth), mcWeight);
        hjetdr[isys][ir]->Fill(maxjet[ir].deltaR(maxjet_truth[ir]), mcWeight);

        hphopurden[isys][ir]->Fill(maxpho.pt, mcWeight);
        hphoeffden[isys][ir]->Fill(maxpho_truth.pt, mcWeight);
        if (check_match(maxpho, maxpho_truth)) {
          hphopurnum[isys][ir]->Fill(maxpho.pt, mcWeight);
          hphoeffnum[isys][ir]->Fill(maxpho_truth.pt, mcWeight);
          // Photon-only quantity, independent of jet radius - fill once (ir==1, same
          // convention as hpurity_num/den below) rather than nJetR times.
          if (ir == 1) {
            hphoIDeff_bdt[isys]->Fill(maxpho_truth.pt, maxpho.bdt, mcWeight);
            if (maxpho.iso4 > -999) hphoIDeff_iso[isys]->Fill(maxpho_truth.pt, maxpho.iso4, mcWeight);
          }
        }
        hjetpurden[isys][ir]->Fill(maxjet[ir].pt, mcWeight);
        hjeteffden[isys][ir]->Fill(maxjet_truth[ir].pt, mcWeight);
        if (check_match(maxjet[ir], maxjet_truth[ir])) {
          hjetpurnum[isys][ir]->Fill(maxjet[ir].pt, mcWeight);
          hjeteffnum[isys][ir]->Fill(maxjet_truth[ir].pt, mcWeight);
        }

        float xj = maxjet[ir].pt/maxpho.pt;
        float xj_truth = maxjet_truth[ir].pt/maxpho_truth.pt;
        int bin = ana::findUnfoldBin(xj,maxpho.pt);
        int bin_truth = ana::findUnfoldBin(xj_truth,maxpho_truth.pt);


        // -----------------------
        // Pairing
        // -----------------------
        // Both bounds must match ana::findPtBin's half-open [ptBins[i], ptBins[i+1]) convention,
        // not just the lower one: findUnfoldBin (via findPtBin) accepts pt == ptBins[0] (>=) but
        // rejects pt >= ptBins[nPtBins] (<), returning bin -1 outside that range. A pt that
        // doesn't match this exact convention can still be "paired" but land in no valid unfold
        // bin, which silently drops the event into RooUnfoldResponse's histogram underflow
        // (counted in _mes/_tru totals but not in the response matrix used for inversion, and
        // not routed through Miss()/Fake() either) - a bookkeeping leak.
        if (maxpho.pt >= ana::ptBins[0] && maxpho.pt < ana::ptBins[ana::nPtBins] && maxjet[ir].pt > ana::jet_calib_pt_cut[ir]) {
          ispaired[ir] = check_pair(maxjet[ir], ir, maxpho,1);
        }
        // Same pairing logic as ispaired above, but against rawJetPt with a floor scale
        // tied to insitu_utility::scanLow - the grid scan's lowest trial pa - so every
        // event a scan point down to scanLow could scale above the floor actually makes
        // it into insitutree (see insitu_utility.h; only the photon-pt bounds guard the
        // ptbin==-1 index inside check_pair, not the jet-pt pre-filter, which is a
        // nominal-floor shortcut and would needlessly exclude the very low-side events
        // this is meant to keep).
        if (maxpho.pt >= ana::ptBins[0] && maxpho.pt < ana::ptBins[ana::nPtBins]) {
          ispairedInsitu[ir] = check_pair(maxjet[ir], ir, maxpho, 1, rawJetPt, insitu_utility::scanLow);
        }
        // threejet: reco-only veto (no truth-level third-jet branch exists) - applied on
        // top of the nominal pairing requirement, before it feeds the ABCD fills below.
        if (systagThreejetVeto) ispaired[ir] = ispaired[ir] && !hasthirdjet[ir];
        if (systagThreejetVeto) ispairedInsitu[ir] = ispairedInsitu[ir] && !hasthirdjet[ir];
        if (maxpho_truth.pt >= ana::ptBins[0] && maxpho_truth.pt < ana::ptBins[ana::nPtBins] && maxjet_truth[ir].pt > ana::jet_calib_pt_cut[ir]) {
          ispaired_truth[ir] = check_pair(maxjet_truth[ir], ir, maxpho_truth,1);
        }
        // narrowBDT/narrowISO reselect the reco ABCD grid with a tighter working point
        // (ana::findabcdBin bin 1/2 instead of nominal bin 0) across all four regions, not
        // just region A, since purity correction downstream needs A and C together. The
        // truth-side ABCD stays pinned to bin 0 for every systag, including narrowBDT/
        // narrowISO: the response matrix's truth axis is the fixed fiducial definition
        // being measured, and a reconstruction/selection systematic should vary how well
        // that fixed target is reconstructed, not the target itself. (Truth photons do
        // have a real, non-trivial isolation spread - this is a deliberate choice, not an
        // invariant simplification.)
        int iabcd_reco = ana::findabcdBin(maxpho.iso4, maxpho.bdt, systagAbcdBin);
        int iabcd_truth = ana::findabcdBin(maxpho_truth.iso4, maxpho_truth.bdt, 0);
        if (ispaired[ir] && iabcd_reco != -1) hrecoxj_abcd[isys][ir][iabcd_reco]->Fill(bin, mcWeight);
        if (ispaired_truth[ir] && iabcd_truth != -1) htruthxj_abcd[isys][ir][iabcd_truth]->Fill(bin, mcWeight);
        if (ispaired[ir] && iabcd_reco != -1) hclusterpt_abcd[isys][ir][iabcd_reco]->Fill(maxpho.pt, mcWeight);

        // In-situ test tree: one per systag (not per radius - see insitu_tree's
        // construction in unfolder.h), an "ir" branch distinguishes which jet radius
        // each row is for. All four ABCD regions, gated on ispairedInsitu (eta/dphi/
        // ptbin cuts identical to the response matrix's ispaired, but the xJ floor is
        // tested against rawJetPt at insitu_utility::scanLow instead of the
        // jesCorrection-boosted recoJetPt at floorScale=1.0) - not yet narrowed to the
        // signal region below.
        //
        // insitu_weight carries the same vz/cluster_pt mcWeight used everywhere else in
        // this event loop (1.0 for Data, since mcWeight is only ever non-trivial for
        // isMC) - insitu/grid_insitu.C and friends must multiply it into every MC fill
        // they do from this tree so the in-situ reference shape matches the rest of the
        // pipeline's Data/MC reweighting instead of silently using an unweighted MC shape.
        if (ispairedInsitu[ir] && iabcd_reco != -1) {
          insitu_pho_pt[isys] = maxpho.pt;
          insitu_jet_pt[isys] = rawJetPt;
          insitu_abcd[isys] = iabcd_reco;
          insitu_weight[isys] = mcWeight;
          insitu_ir[isys] = ir;
          insitu_tree[isys]->Fill();
        }

        // hpurity_num/den (below) are NOT the real purity pipeline - they're written to
        // the unfolding output file but nothing downstream reads them. The actual
        // purity determination is hclusterpt_abcd[isys][i][j] (filled per radius above,
        // gated on ispaired[ir]) -> macros/puritymaker.C -> ana::getPurity/getPurityC
        // (which now take an ir argument - purity is "of paired photons", and pairing
        // genuinely differs by jet radius, so it's not one number reused everywhere).
        // hpurity_num/den stay a single representative radius (R=0.4, ir==2) since nothing
        // consumes them per-radius; kept only for whatever manual/future inspection they
        // were originally added for.
        if (ir == 2 && ispairedInsitu[ir] && iabcd_reco != -1) {
          bool maxpho_is_photon = check_match(maxpho, maxpho_truth);
          int iabcd = ana::findabcdBin(maxpho.iso4, maxpho.bdt, 0);
          float xj = maxjet[ir].pt/maxpho.pt;
          if (iabcd == 0) {
            if (maxpho_is_photon) {
              hpurity_num[isys]->Fill(maxpho.pt, xj, mcWeight);
              hpurity_num_1D[isys]->Fill(maxpho.pt, mcWeight);
            }
            hpurity_den[isys]->Fill(maxpho.pt,xj, mcWeight);
            hpurity_den_1D[isys]->Fill(maxpho.pt, mcWeight);
          }
        }


        // Response matrix / pair purity-efficiency counters below are restricted to the
        // signal region (abcd bin 0) on both reco and truth sides - the abcd-inclusive
        // ispaired/ispaired_truth above was only needed for the hrecoxj_abcd/htruthxj_abcd/
        // ispairedInsitu fills.
        ispaired[ir] = ispaired[ir] && iabcd_reco == 0;
        ispaired_truth[ir] = ispaired_truth[ir] && iabcd_truth == 0;

        // -----------------------
        // Fill response matrix
        // -----------------------
        bool ismatch = ispaired_truth[ir] && ispaired[ir] && check_match(maxjet[ir], maxjet_truth[ir]);
        if (ispaired_truth[ir] && ispaired[ir] && ismatch) {
          jet_response[isys][ir]->Fill(maxjet[ir].pt, maxjet_truth[ir].pt, mcWeight);
          pho_response[isys][ir]->Fill(maxpho.pt, maxpho_truth.pt, mcWeight);
          jet_response2D[isys][ir]->Fill(bin, bin_truth, mcWeight);

          hphomissfake[isys][ir]->Fill(maxpho.pt, maxpho_truth.pt, mcWeight);
          hjetmissfake[isys][ir]->Fill(maxjet[ir].pt, maxjet_truth[ir].pt, mcWeight);
          hpairmissfake[isys][ir]->Fill(bin, bin_truth, mcWeight);

          if (use_half) {
            jet_response_half[isys][ir]->Fill(maxjet[ir].pt, maxjet_truth[ir].pt, mcWeight);
            pho_response_half[isys][ir]->Fill(maxpho.pt, maxpho_truth.pt, mcWeight);
            jet_response_half2D[isys][ir]->Fill(bin, bin_truth, mcWeight);
          }
        }
        else if (ispaired_truth[ir] && ispaired[ir] && !ismatch) {
          jet_response[isys][ir]->Miss(maxjet_truth[ir].pt, mcWeight);
          jet_response[isys][ir]->Fake(maxjet[ir].pt, mcWeight);
          pho_response[isys][ir]->Miss(maxpho_truth.pt, mcWeight);
          pho_response[isys][ir]->Fake(maxpho.pt, mcWeight);
          jet_response2D[isys][ir]->Miss(bin_truth, mcWeight);
          jet_response2D[isys][ir]->Fake(bin, mcWeight);

          hphomissfake[isys][ir]->Fill(maxpho.pt, 100, mcWeight);
          hjetmissfake[isys][ir]->Fill(maxjet[ir].pt, 100, mcWeight);
          hpairmissfake[isys][ir]->Fill(bin, ana::nPtBins*(ana::nUnfoldXjBins+2), mcWeight);
          hphomissfake[isys][ir]->Fill(100, maxpho_truth.pt, mcWeight);
          hjetmissfake[isys][ir]->Fill(100, maxjet_truth[ir].pt, mcWeight);
          hpairmissfake[isys][ir]->Fill(ana::nPtBins*(ana::nUnfoldXjBins+2), bin_truth, mcWeight); // +2 for underflow and overflow bins

          if (use_half) {
            jet_response_half[isys][ir]->Miss(maxjet_truth[ir].pt, mcWeight);
            jet_response_half[isys][ir]->Fake(maxjet[ir].pt, mcWeight);
            pho_response_half[isys][ir]->Miss(maxpho_truth.pt, mcWeight);
            pho_response_half[isys][ir]->Fake(maxpho.pt, mcWeight);
            jet_response_half2D[isys][ir]->Miss(bin_truth, mcWeight);
            jet_response_half2D[isys][ir]->Fake(bin, mcWeight);
          }
        }
        else if (ispaired_truth[ir] && !ispaired[ir]) {
          jet_response[isys][ir]->Miss(maxjet_truth[ir].pt, mcWeight);
          pho_response[isys][ir]->Miss(maxpho_truth.pt, mcWeight);
          jet_response2D[isys][ir]->Miss(bin_truth, mcWeight);
          hphomissfake[isys][ir]->Fill(100, maxpho_truth.pt, mcWeight);
          hjetmissfake[isys][ir]->Fill(100, maxjet_truth[ir].pt, mcWeight);
          hpairmissfake[isys][ir]->Fill(ana::nPtBins*(ana::nUnfoldXjBins+2), bin_truth, mcWeight);
          if (use_half) {
            jet_response_half[isys][ir]->Miss(maxjet_truth[ir].pt, mcWeight);
            pho_response_half[isys][ir]->Miss(maxpho_truth.pt, mcWeight);
            jet_response_half2D[isys][ir]->Miss(bin_truth, mcWeight);
          }
        }
        else if (!ispaired_truth[ir] && ispaired[ir]) {
          jet_response[isys][ir]->Fake(maxjet[ir].pt, mcWeight);
          pho_response[isys][ir]->Fake(maxpho.pt, mcWeight);
          jet_response2D[isys][ir]->Fake(bin, mcWeight);
          hphomissfake[isys][ir]->Fill(maxpho.pt, 100, mcWeight);
          hjetmissfake[isys][ir]->Fill(maxjet[ir].pt, 100, mcWeight);
          hpairmissfake[isys][ir]->Fill(bin, ana::nPtBins*(ana::nUnfoldXjBins+2), mcWeight);
          if (use_half) {
            jet_response_half[isys][ir]->Fake(maxjet[ir].pt, mcWeight);
            pho_response_half[isys][ir]->Fake(maxpho.pt, mcWeight);
            jet_response_half2D[isys][ir]->Fake(bin, mcWeight);
          }
        }

        int ipt = ana::findPtBin(maxpho.pt);

        if (ispaired[ir]) {
          hpairpurden[isys][ir]->Fill(bin, mcWeight);

          hrecojetpt[isys][ir]->Fill(maxjet[ir].pt, mcWeight);
          hrecophopt[isys][ir]->Fill(maxpho.pt, mcWeight);
          hrecoxj[isys][ir]->Fill(bin, mcWeight);
          if (!use_half) {
            hrecojetpt_half[isys][ir]->Fill(maxjet[ir].pt, mcWeight);
            hrecophopt_half[isys][ir]->Fill(maxpho.pt, mcWeight);
            hrecoxj_half[isys][ir]->Fill(bin, mcWeight);
          }
        }

        if (ispaired_truth[ir]) {
          hpaireffden[isys][ir]->Fill(bin_truth, mcWeight);

          htruthjetpt[isys][ir]->Fill(maxjet_truth[ir].pt, mcWeight);
          htruthphopt[isys][ir]->Fill(maxpho_truth.pt, mcWeight);
          htruthxj[isys][ir]->Fill(bin_truth, mcWeight);
          if (!use_half) {
            htruthjetpt_half[isys][ir]->Fill(maxjet_truth[ir].pt, mcWeight);
            htruthphopt_half[isys][ir]->Fill(maxpho_truth.pt, mcWeight);
            htruthxj_half[isys][ir]->Fill(bin_truth, mcWeight);
          }
        }

        if (ismatch) {
          hpairpurnum[isys][ir]->Fill(bin, mcWeight);
          hpaireffnum[isys][ir]->Fill(bin_truth, mcWeight);
        }


        // Drawing event displays - only for the first systag in the list, so enabling
        // dodraw doesn't multiply the debug PDF output nsys-fold (event displays are a
        // manual debugging aid, not part of the physics output).
        if (dodraw && isys == 0 && ndraw < 100 && ir == 1 && (!ispaired_truth[ir] && ispaired[ir])) {
          if (!ispaired[ir]) {
            float dphi = maxjet[ir].deltaPhi(maxpho);
            int iabcd = ana::findabcdBin(maxpho.iso4, maxpho.bdt, 0);
            cout << "Event " << ndraw + 1 << " reco failed because: ";
            if (iabcd != 0) cout << endl << "abcd cut: BDT: " << maxpho.bdt << " ISO: " << maxpho.iso4;
            if (fabs(maxpho.eta) >= ana::etacut) cout << endl << "photon eta: " << fabs(maxpho.eta);
            if (fabs(maxjet[ir].eta) >= ana::etacut - ana::JetRs[ir]) cout << endl << "jet eta: " << fabs(maxjet[ir].eta);
            if (dphi <= ana::oppcut) cout << endl << "dphi: " << dphi;
            if (maxpho.pt <= ana::ptBins[0]) cout << endl << "cluster pt: " << maxpho.pt;
            if (maxjet[ir].pt <= ana::jet_calib_pt_cut[ir]) cout << endl << "jet pt: " << maxjet[ir].pt;
            cout << endl;
          }

          if (!ispaired_truth[ir]) {
            float dphi = maxjet_truth[ir].deltaPhi(maxpho_truth);
            int iabcd = ana::findabcdBin(maxpho_truth.iso4, maxpho_truth.bdt, 0);
            cout << "Event " << ndraw + 1 << " truth failed because: ";
            if (iabcd != 0) cout << endl << "abcd cut: BDT: " << maxpho_truth.bdt << " ISO: " << maxpho_truth.iso4;
            if (fabs(maxpho_truth.eta) >= ana::etacut) cout << endl << "photon eta: " << fabs(maxpho_truth.eta);
            if (fabs(maxjet_truth[ir].eta) >= ana::etacut - ana::JetRs[ir]) cout << endl << "jet eta: " << fabs(maxjet_truth[ir].eta);
            if (dphi <= ana::oppcut) cout << endl << "dphi: " << dphi;
            if (maxpho_truth.pt <= ana::ptBins[0]) cout << endl << "cluster pt: " << maxpho_truth.pt;
            if (maxjet_truth[ir].pt <= ana::jet_calib_pt_cut[ir]) cout << endl << "jet pt: " << maxjet_truth[ir].pt;
            cout << endl;
          }


          TH2D * h = new TH2D("heventdisplay",";eta;phi",100,-1.5,1.5,100,-M_PI,M_PI);
          h->Draw();

          TMarker * star = new TMarker(maxpho.eta, maxpho.phi, 29);
          star->SetMarkerColor(kRed);
          star->SetMarkerSize(2);
          TMarker * star_truth = new TMarker(maxpho_truth.eta, maxpho_truth.phi, 29);
          star_truth->SetMarkerColor(kBlack);
          star_truth->SetMarkerSize(3);

          TMarker * circle = new TMarker(maxjet[ir].eta, maxjet[ir].phi, 20);
          circle->SetMarkerColorAlpha(kRed,0.5);
          circle->SetMarkerSize(15);
          TMarker * circle_truth = new TMarker(maxjet_truth[ir].eta, maxjet_truth[ir].phi, 20);
          circle_truth->SetMarkerColorAlpha(kBlack, 0.5);
          circle_truth->SetMarkerSize(15);

          // For the legend
          TMarker * dummy_truth = new TMarker(maxpho_truth.eta, maxpho_truth.phi, 29);
          dummy_truth->SetMarkerColor(kBlack);
          dummy_truth->SetMarkerSize(2);
          TMarker * dummy_circle = new TMarker(maxpho_truth.eta, maxpho_truth.phi, 20);
          dummy_circle->SetMarkerColorAlpha(kRed,0.5);
          dummy_circle->SetMarkerSize(2);
          TMarker * dummy_circle_truth = new TMarker(maxpho_truth.eta, maxpho_truth.phi, 20);
          dummy_circle_truth->SetMarkerColorAlpha(kBlack,0.5);
          dummy_circle_truth->SetMarkerSize(2);

          if (maxpho_truth.pt > 0) star_truth->Draw();
          if (maxpho.pt > 0) star->Draw();
          if (maxjet[ir].pt > 0) circle->Draw();
          if (maxjet_truth[ir].pt > 0) circle_truth->Draw();

          TLegend * l = new TLegend(0,0.8,.3,1);
          l->AddEntry(star, "Reco Cluster");
          l->AddEntry(dummy_truth, "Truth Cluster");
          l->AddEntry(dummy_circle, "Reco Jet");
          l->AddEntry(dummy_circle_truth, "Truth Jet");
          l->Draw();

          TLatex latex;
          latex.SetNDC();

          if (!ispaired[ir]) latex.DrawLatex(0.3,0.95,"Reco Failed");
          if (!ispaired_truth[ir]) latex.DrawLatex(0.3,0.9,"Truth Failed");

          ndraw++;
          c->SaveAs(Form("/home/samson72/sphnx/gammajet_unfold/pdfs/event_displays_%s.pdf",trigger.c_str()));
          delete h;
          delete l;
          c->Clear();

          if (ndraw == 100) c->SaveAs(Form("/home/samson72/sphnx/gammajet_unfold/pdfs/event_displays_%s.pdf]",trigger.c_str()));
        }
      }
    }
  }
}

void unfolder::unfold() {
  // -----------------------
  // Unfold
  // -----------------------
  int nsys = (int)systags.size();

  for (int isys = 0; isys < nsys; isys++) {
    for (int ir = 0; ir < ana::nJetR; ir++) {
      // Full Closure
      RooUnfoldBayes jet_unfold(jet_response[isys][ir], hrecojetpt[isys][ir], niterate, 0, 1);
      jet_unfold.SetVerbose(-1);
      hunfoldjetpt[isys][ir] = (TH1D*)jet_unfold.Hreco();
      hunfoldjetpt[isys][ir]->SetName(Form("hunfoldjetpt%i", ir));
      hjetresponse[isys][ir] = (TH2D*)jet_response[isys][ir]->Hresponse();
      hjetresponse[isys][ir]->SetName(Form("hjetresponse%i", ir));

      RooUnfoldBayes pho_unfold(pho_response[isys][ir], hrecophopt[isys][ir], niterate, 0, 1);
      pho_unfold.SetVerbose(-1);
      hunfoldphopt[isys][ir] = (TH1D*)pho_unfold.Hreco();
      hunfoldphopt[isys][ir]->SetName(Form("hunfoldphopt%i",ir));
      hphoresponse[isys][ir] = (TH2D*)pho_response[isys][ir]->Hresponse();
      hphoresponse[isys][ir]->SetName(Form("hphoresponse%i",ir));

      RooUnfoldBayes jet_unfold2D(jet_response2D[isys][ir], hrecoxj[isys][ir], niterate, 0, 1);
      jet_unfold2D.SetVerbose(-1);
      hunfoldxj[isys][ir] = (TH1D*)jet_unfold2D.Hreco();
      hunfoldxj[isys][ir]->SetName(Form("hunfoldxj%i",  ir));
      hxjresponse[isys][ir] = (TH2D*)jet_response2D[isys][ir]->Hresponse();
      hxjresponse[isys][ir]->SetName(Form("hxjresponse%i", ir));

      // Half closure
      RooUnfoldBayes jet_unfold_half(jet_response_half[isys][ir], hrecojetpt_half[isys][ir], niterate, 0, 1);
      jet_unfold_half.SetVerbose(-1);
      hunfoldjetpt_half[isys][ir] = (TH1D*)jet_unfold_half.Hreco();
      hunfoldjetpt_half[isys][ir]->SetName(Form("hunfoldjetpt_half%i", ir));
      hjetresponse_half[isys][ir] = (TH2D*)jet_response_half[isys][ir]->Hresponse();
      hjetresponse_half[isys][ir]->SetName(Form("hjetresponse_half%i", ir));

      RooUnfoldBayes pho_unfold_half(pho_response_half[isys][ir], hrecophopt_half[isys][ir], niterate, 0, 1);
      pho_unfold_half.SetVerbose(-1);
      hunfoldphopt_half[isys][ir] = (TH1D*)pho_unfold_half.Hreco();
      hunfoldphopt_half[isys][ir]->SetName(Form("hunfoldphopt_half%i",ir));
      hphoresponse_half[isys][ir] = (TH2D*)pho_response_half[isys][ir]->Hresponse();
      hphoresponse_half[isys][ir]->SetName(Form("hphoresponse_half%i",ir));


      RooUnfoldBayes jet_unfold_half2D(jet_response_half2D[isys][ir], hrecoxj_half[isys][ir], niterate, 0, 1);
      jet_unfold_half2D.SetVerbose(-1);
      hunfoldxj_half[isys][ir] = (TH1D*)jet_unfold_half2D.Hreco();
      hunfoldxj_half[isys][ir]->SetName(Form("hunfoldxj_half%i", ir));
      hxjresponse_half[isys][ir] = (TH2D*)jet_response_half2D[isys][ir]->Hresponse();
      hxjresponse_half[isys][ir]->SetName(Form("hxjresponse_half%i", ir));
    }
  }
}

void unfolder::savehists(TH1D * h[], int n) {
  for (int i = 0; i < n; i++) {
    h[i]->Write();
  }
}
void unfolder::savehists(TH2D * h[], int n) {
  for (int i = 0; i < n; i++) {
    h[i]->Write();
  }
}
void unfolder::savehists(RooUnfoldResponse * h[], int n) {
  for (int i = 0; i < n; i++) {
    h[i]->Write();
  }
}
void unfolder::savehists(TEfficiency * h[], int n) {
  for (int i = 0; i < n; i++) {
    h[i]->Write();
  }
}

void unfolder::end() {
  int nsys = (int)systags.size();

  for (int isys = 0; isys < nsys; isys++) {
    insitu_file[isys]->cd();
    insitu_tree[isys]->Write();
    insitu_file[isys]->Close();

    const string & systag = systags[isys];
    const char * wfilename = Form("/home/samson72/sphnx/gammajet_unfold/hists/%s_%s_%s_unfolding.root",trigger.c_str(),sim.c_str(),systag.c_str());
    if (!isMC) wfilename = Form("/home/samson72/sphnx/gammajet_unfold/hists/%s_%s_unfolding.root",trigger.c_str(),systag.c_str());
    cout << "Writing files to " << wfilename << endl;
    TFile * fout = TFile::Open(wfilename, "RECREATE");

    savehists(hphodr[isys].data(),ana::nJetR);
    savehists(hjetdr[isys].data(),ana::nJetR);

    savehists(hrecojetpt[isys].data(),ana::nJetR);
    savehists(htruthjetpt[isys].data(),ana::nJetR);
    savehists(hunfoldjetpt[isys].data(),ana::nJetR);
    savehists(hjetresponse[isys].data(),ana::nJetR);
    savehists(hrecophopt[isys].data(),ana::nJetR);
    savehists(htruthphopt[isys].data(),ana::nJetR);
    savehists(hunfoldphopt[isys].data(),ana::nJetR);
    savehists(hphoresponse[isys].data(),ana::nJetR);

    savehists(hrecojetpt_half[isys].data(),ana::nJetR);
    savehists(htruthjetpt_half[isys].data(),ana::nJetR);
    savehists(hunfoldjetpt_half[isys].data(),ana::nJetR);
    savehists(hjetresponse_half[isys].data(),ana::nJetR);
    savehists(hrecophopt_half[isys].data(),ana::nJetR);
    savehists(htruthphopt_half[isys].data(),ana::nJetR);
    savehists(hunfoldphopt_half[isys].data(),ana::nJetR);
    savehists(hphoresponse_half[isys].data(),ana::nJetR);

    savehists(hrecoxj[isys].data(),ana::nJetR);
    savehists(htruthxj[isys].data(),ana::nJetR);
    savehists(hunfoldxj[isys].data(),ana::nJetR);
    savehists(hxjresponse[isys].data(),ana::nJetR);

    savehists(hrecoxj_half[isys].data(),ana::nJetR);
    savehists(htruthxj_half[isys].data(),ana::nJetR);
    savehists(hunfoldxj_half[isys].data(),ana::nJetR);
    savehists(hxjresponse_half[isys].data(),ana::nJetR);

    savehists(jet_response2D[isys].data(),ana::nJetR);
    savehists(jet_response_half2D[isys].data(),ana::nJetR);

    savehists(hphopurden[isys].data(),ana::nJetR);
    savehists(hphopurnum[isys].data(),ana::nJetR);
    savehists(hphoeffden[isys].data(),ana::nJetR);
    savehists(hphoeffnum[isys].data(),ana::nJetR);
    savehists(hjetpurden[isys].data(),ana::nJetR);
    savehists(hjetpurnum[isys].data(),ana::nJetR);
    savehists(hjeteffden[isys].data(),ana::nJetR);
    savehists(hjeteffnum[isys].data(),ana::nJetR);
    savehists(hpairpurden[isys].data(),ana::nJetR);
    savehists(hpairpurnum[isys].data(),ana::nJetR);
    savehists(hpaireffden[isys].data(),ana::nJetR);
    savehists(hpaireffnum[isys].data(),ana::nJetR);

    for (int i = 0; i < ana::nJetR; i++) {
      for (int j = 0; j < 4; j++) {
        hrecoxj_abcd[isys][i][j]->Write();
        htruthxj_abcd[isys][i][j]->Write();
        hclusterpt_abcd[isys][i][j]->Write();
      }
    }

    for (int ir = 0; ir < ana::nJetR; ir++) {
      hphoeff[isys][ir] = new TEfficiency(*hphoeffnum[isys][ir], *hphoeffden[isys][ir]);
      hphopur[isys][ir] = new TEfficiency(*hphopurnum[isys][ir], *hphopurden[isys][ir]);
      hjeteff[isys][ir] = new TEfficiency(*hjeteffnum[isys][ir], *hjeteffden[isys][ir]);
      hjetpur[isys][ir] = new TEfficiency(*hjetpurnum[isys][ir], *hjetpurden[isys][ir]);
      hpaireff[isys][ir] = new TEfficiency(*hpaireffnum[isys][ir], *hpaireffden[isys][ir]);
      hpairpur[isys][ir] = new TEfficiency(*hpairpurnum[isys][ir], *hpairpurden[isys][ir]);

      // num/den are filled with mcWeight (see fill_matrix()) - the default kFCP
      // (Clopper-Pearson) interval assumes unweighted integer pass/total counts, so it's
      // not valid once entries carry a weight. kFNormal is the standard fallback for
      // weighted TEfficiency; SetUseWeightedEvents() tells it to use the num/den
      // histograms' Sumw2 errors instead of raw bin content for the variance.
      for (TEfficiency * eff : {hphoeff[isys][ir], hphopur[isys][ir], hjeteff[isys][ir],
                                 hjetpur[isys][ir], hpaireff[isys][ir], hpairpur[isys][ir]}) {
        eff->SetStatisticOption(TEfficiency::kFNormal);
        eff->SetUseWeightedEvents();
      }

      hphoeff [isys][ir]->SetName(Form("hphoeff%i",ir));
      hphopur [isys][ir]->SetName(Form("hphopur%i",ir));
      hjeteff [isys][ir]->SetName(Form("hjeteff%i",ir));
      hjetpur [isys][ir]->SetName(Form("hjetpur%i",ir));
      hpaireff[isys][ir]->SetName(Form("hpaireff%i",ir));
      hpairpur[isys][ir]->SetName(Form("hpairpur%i",ir));
    }
    savehists(hphoeff[isys].data(),ana::nJetR);
    savehists(hphopur[isys].data(),ana::nJetR);
    savehists(hjeteff[isys].data(),ana::nJetR);
    savehists(hjetpur[isys].data(),ana::nJetR);
    savehists(hpaireff[isys].data(),ana::nJetR);
    savehists(hpairpur[isys].data(),ana::nJetR);

    savehists(hphomissfake[isys].data(),ana::nJetR);
    savehists(hjetmissfake[isys].data(),ana::nJetR);
    savehists(hpairmissfake[isys].data(),ana::nJetR);

    hpurity_num[isys]->Write();
    hpurity_den[isys]->Write();
    hpurity_num_1D[isys]->Write();
    hpurity_den_1D[isys]->Write();

    hphoIDeff_bdt[isys]->Write();
    hphoIDeff_iso[isys]->Write();

    fout->Close();
  }
}
