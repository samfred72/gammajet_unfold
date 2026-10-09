#include "unfolder.h"
#include "insitu_utility.h"
using namespace std;

unfolder::~unfolder() {}

bool unfolder::check_pair(jet_object jet, int ir, pho_object pho, bool isreco, float testPt, float floorScale) {
  float dphi = jet.deltaPhi(pho);
  int iabcd = ana::findabcdBin(pho.iso4, pho.bdt, 0);

  int ptbin = ana::findPtBin(pho.pt);
  // testPt: test the x_J floor against another pT (the in-situ tree's raw jet pT); floorScale scales the floor.
  float val = (testPt >= 0 ? testPt : jet.pt)/pho.pt;
  float lowval = ana::jet_calib_pt_cut[ir]/ana::ptBins[ptbin];
  float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(lowval)+1];

  if (ptbin == -1) return false;
  if (val < lowbin*floorScale) return false;
  if (fabs(pho.eta) > ana::photonEtaCut) return false;
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

  // Per-systag constants, computed once (see unfolder.h).
  vector<int> systagAbcdBinArr(nsys);
  vector<bool> systagThreejetVetoArr(nsys);
  vector<float> systagEmscaleShiftArr(nsys);
  vector<int> systagEmrVariantArr(nsys);
  vector<float> systagTimingWidenArr(nsys);
  vector<vector<float>> jesCorrectionArr(nsys, vector<float>(ana::nJetR));
  for (int isys = 0; isys < nsys; isys++) {
    const string & systag = systags[isys];
    systagAbcdBinArr[isys] = (systag == "narrowBDT") ? 1 : (systag == "narrowISO") ? 2 :
                              (systag == "narrowBDTbkg") ? 3 : (systag == "narrowISObkg") ? 4 :
                              (systag == "wideISObkg") ? 5 : 0;
    systagThreejetVetoArr[isys] = (systag == "threejet");
    // timingwide (Data only): both timing cuts widened by ana::timingSystWiden.
    systagTimingWidenArr[isys] = (systag == "timingwide") ? ana::timingSystWiden : 0.0;
    // emscale_high/low (MC only): +-ana::emscaleShift on the photon and on the EM fraction of the jet.
    systagEmscaleShiftArr[isys] = (systag == "emscale_high") ?  ana::emscaleShift :
                                   (systag == "emscale_low")  ? -ana::emscaleShift : 0.0;
    // EMRhigh/EMRlow (MC only): the extra cluster smearing uses the wider or no Data resolution
    // (ana::emResolutionSigma). Independent of emscale.
    systagEmrVariantArr[isys] = (systag == "EMRhigh") ? ana::emrHigh :
                                 (systag == "EMRlow")  ? ana::emrLow  : ana::emrNominal;
    // JES: Data jet pT is divided by the in-situ p_a (ana.h). A systag with its own in-situ scan uses
    // its own p_a (ana::jesForSystag); jes_high/low shift the nominal p_a by its statistical error.
    // Not applied to the in-situ tree's jet pT (rawJetPt).
    for (int ir = 0; ir < ana::nJetR; ir++) {
      jesCorrectionArr[isys][ir] =
          (systag == "jes_high") ? ana::jesNominal[ir] - ana::jesStatErrLow[ir] :
          (systag == "jes_low")  ? ana::jesNominal[ir] + ana::jesStatErrHigh[ir] :
          ana::jesForSystag(systag, ir);
    }
  }

  TCanvas * c = new TCanvas("c","",500,1000);
  gStyle->SetOptStat(0);
  if (dodraw) c->SaveAs(Form("%s/pdfs/event_displays_%s.pdf[", ana::dir(),trigger.c_str()));
  int ndraw = 0;

  for (Long64_t e = 0; e < nentries; e++) {
    t->GetEntry(e);
    bool use_half = rand.Integer(2) % 2;
    // One standard-normal draw per MC event for the photon EM-resolution smear, shared by every systag
    // (each scales it by its own sigma), so the variations differ only by what they change, not by a
    // new random realization of the smear.
    float emSmearZ = isMC ? rand.Gaus(0, 1) : 0;
    if (showProgress && e % 1000 == 0)
      std::cout << "entry " << e << "/" << nentries
        << " (" << (float)e/nentries*100. << "%)\t\r" << std::flush;

    if (fabs(vz) > ana::vzcut) continue;

    // -----------------------
    // Event selection (systag-independent)
    // -----------------------

    vector<bool> keepMC = check_keep_MC(truth_cluster_pt, cluster_pt, truth_jet_pt, jet_pt_smear_truth, trigger);
    if (isMC && !keepMC.at(keepMC.size()-1)) continue;

    // Data/MC vz and cluster-pT weight from the raw quantities, the same for every systag; 1 for Data.
    float mcWeight = isMC ? rw.GetWeight(vz, cluster_pt) : 1.0f;

    // -----------------------
    // Per-systag reprocessing of this event (closed-form shifts or branch choices; no re-read)
    // -----------------------
    for (int isys = 0; isys < nsys; isys++) {
      const string & systag = systags[isys];
      int systagAbcdBin = systagAbcdBinArr[isys];
      bool systagThreejetVeto = systagThreejetVetoArr[isys];
      float systagEmscaleShift = systagEmscaleShiftArr[isys];
      int systagEmrVariant = systagEmrVariantArr[isys];

      // -----------------------
      // Leading photon & isolation
      // -----------------------
      float recoClusterPt = isMC ? cluster_pt * (1.0 + systagEmscaleShift) : cluster_pt;
      // Additive smear N(0, sigma_extra(E_truth)*E_truth), PPG12 convention (ana.h), from the event's draw.
      recoClusterPt = isMC ? recoClusterPt + emSmearZ*ana::emResolutionSigma(truth_cluster_pt, systagEmrVariant)*truth_cluster_pt : cluster_pt;
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
      // In-situ pairing: x_J floor tested against the uncorrected jet pT at floorScale = scanLow, so the
      // tree keeps every event a scan point down to scanLow could use.
      vector<bool> ispairedInsitu(ana::nJetR, false);

      for (int ir = 0; ir < ana::nJetR; ir++) {
        if (isMC && !keepMC[ir]) continue;
        // Data timing cut (ana::timing*); a failing event is dropped for this radius.
        if (!isMC) {
          if (fabs(cluster_time - ana::timingClusterCenter) >= ana::timingClusterHalfWidth + systagTimingWidenArr[isys]) continue;
          if (jet_pt_calib[ir] > 0 && fabs(cluster_time - jet_time[ir]) >= ana::timingDeltaMax + systagTimingWidenArr[isys]) continue;
        }

        // JER/emscale variations are MC only; jes_high/low are Data only.
        float recoJetPt;
        if (isMC) {
          recoJetPt = (systag == "JERhigh") ? jet_pt_smear_high_truth[ir] :
                      (systag == "JERlow")  ? jet_pt_smear_low_truth[ir]  :
                      // Scale only the EM fraction of the jet.
                      (systagEmscaleShift != 0.0) ?
                          jet_pt_smear_truth[ir]*jet_emfrac[ir]
                            + jet_pt_smear_truth[ir]*jet_emfrac[ir]*systagEmscaleShift
                            + jet_pt_smear_truth[ir]*(1-jet_emfrac[ir]) :
                      jet_pt_smear_truth[ir];
        } else {
          recoJetPt = jet_pt_calib[ir] / jesCorrectionArr[isys][ir];
        }
        // Uncorrected Data jet pT for the in-situ tree: the scan measures the full Data/MC gap, not a
        // residual around the current correction.
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
          // Photon-only: fill once (ir == 1).
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
        // Both photon-pT bounds follow findPtBin's half-open [low, high) convention; otherwise a paired
        // event could land in no unfold bin and leak into RooUnfold's underflow.
        if (maxpho.pt >= ana::ptBins[0] && maxpho.pt < ana::ptBins[ana::nPtBins] && maxjet[ir].pt > ana::jet_calib_pt_cut[ir]) {
          ispaired[ir] = check_pair(maxjet[ir], ir, maxpho,1);
        }
        // Same as ispaired but with the in-situ floor (see above). Only the photon-pT bounds guard
        // check_pair's ptbin index; the nominal jet-pT pre-filter would drop events the scan needs.
        if (maxpho.pt >= ana::ptBins[0] && maxpho.pt < ana::ptBins[ana::nPtBins]) {
          ispairedInsitu[ir] = check_pair(maxjet[ir], ir, maxpho, 1, rawJetPt, insitu_utility::scanLow);
        }
        // threejet: veto if the third jet's pT (recoil-jet definition: Data in-situ corrected, MC nominal
        // smear) exceeds ana::thirdJetPtCut. Reco only: a truth veto would change the observable.
        if (systagThreejetVeto) {
          float thirdJetPt = isMC ? thirdjet_pt[ir] : thirdjet_pt[ir] / jesCorrectionArr[isys][ir];
          bool hasThirdJet = thirdJetPt > ana::thirdJetPtCut;
          ispaired[ir] = ispaired[ir] && !hasThirdJet;
          // Data in-situ tree: the veto depends on p_a, so it is deferred to the reader (thirdjet_pt).
          if (isMC) ispairedInsitu[ir] = ispairedInsitu[ir] && !hasThirdJet;
        }
        if (maxpho_truth.pt >= ana::ptBins[0] && maxpho_truth.pt < ana::ptBins[ana::nPtBins] && maxjet_truth[ir].pt > ana::jet_calib_pt_cut[ir]) {
          ispaired_truth[ir] = check_pair(maxjet_truth[ir], ir, maxpho_truth,1);
        }
        // ABCD variations move one boundary for all four regions (purity needs A and C together). The
        // truth-side ABCD stays nominal for every systag: the truth axis is the fixed fiducial definition.
        int iabcd_reco = ana::findabcdBin(maxpho.iso4, maxpho.bdt, systagAbcdBin);
        int iabcd_truth = ana::findabcdBin(maxpho_truth.iso4, maxpho_truth.bdt, 0);
        if (ispaired[ir] && iabcd_reco != -1) hrecoxj_abcd[isys][ir][iabcd_reco]->Fill(bin, mcWeight);
        if (ispaired_truth[ir] && iabcd_truth != -1) htruthxj_abcd[isys][ir][iabcd_truth]->Fill(bin, mcWeight);
        if (ispaired[ir] && iabcd_reco != -1) hclusterpt_abcd[isys][ir][iabcd_reco]->Fill(maxpho.pt, mcWeight);
        // MC only: truth-matched (dR < 0.1) subset, the signal template for puritymaker.C's leakage fractions.
        if (isMC && ispaired[ir] && iabcd_reco != -1 && check_match(maxpho, maxpho_truth))
          hclusterpt_abcd_truthmatched[isys][ir][iabcd_reco]->Fill(maxpho.pt, mcWeight);

        // In-situ tree: all four ABCD regions, gated on ispairedInsitu. insitu_weight is the MC weight
        // (1 for Data); MC readers must apply it.
        if (ispairedInsitu[ir] && iabcd_reco != -1) {
          insitu_pho_pt[isys] = maxpho.pt;
          insitu_jet_pt[isys] = rawJetPt;
          insitu_abcd[isys] = iabcd_reco;
          insitu_weight[isys] = mcWeight;
          insitu_ir[isys] = ir;
          insitu_third_pt[isys] = (!isMC && systagThreejetVeto) ? thirdjet_pt[ir] : -1;
          insitu_tree[isys]->Fill();
        }

        // hpurity_num/den are not used downstream (the purity comes from hclusterpt_abcd via
        // puritymaker.C); kept at R = 0.4 for inspection.
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

        // Response matrix and pair counters: signal region (abcd 0) only, on both sides.
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

        // Event displays for the first systag only (debugging aid).
        if (dodraw && isys == 0 && ndraw < 100 && ir == 1 && (!ispaired_truth[ir] && ispaired[ir])) {
          if (!ispaired[ir]) {
            float dphi = maxjet[ir].deltaPhi(maxpho);
            int iabcd = ana::findabcdBin(maxpho.iso4, maxpho.bdt, 0);
            cout << "Event " << ndraw + 1 << " reco failed because: ";
            if (iabcd != 0) cout << endl << "abcd cut: BDT: " << maxpho.bdt << " ISO: " << maxpho.iso4;
            if (fabs(maxpho.eta) >= ana::photonEtaCut) cout << endl << "photon eta: " << fabs(maxpho.eta);
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
            if (fabs(maxpho_truth.eta) >= ana::photonEtaCut) cout << endl << "photon eta: " << fabs(maxpho_truth.eta);
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
          c->SaveAs(Form("%s/pdfs/event_displays_%s.pdf", ana::dir(),trigger.c_str()));
          delete h;
          delete l;
          c->Clear();

          if (ndraw == 100) c->SaveAs(Form("%s/pdfs/event_displays_%s.pdf]", ana::dir(),trigger.c_str()));
        }
      }
    }
  }
}

// Extracts each RooUnfoldResponse's matrix (Hresponse). No unfolding here: the plots unfold via
// unfold_utility::unfoldOnce at the nominal niterate.
void unfolder::unfold() {
  int nsys = (int)systags.size();

  for (int isys = 0; isys < nsys; isys++) {
    for (int ir = 0; ir < ana::nJetR; ir++) {
      // Full closure
      hjetresponse[isys][ir] = (TH2D*)jet_response[isys][ir]->Hresponse();
      hjetresponse[isys][ir]->SetName(Form("hjetresponse%i", ir));

      hphoresponse[isys][ir] = (TH2D*)pho_response[isys][ir]->Hresponse();
      hphoresponse[isys][ir]->SetName(Form("hphoresponse%i",ir));

      hxjresponse[isys][ir] = (TH2D*)jet_response2D[isys][ir]->Hresponse();
      hxjresponse[isys][ir]->SetName(Form("hxjresponse%i", ir));

      // Half closure
      hjetresponse_half[isys][ir] = (TH2D*)jet_response_half[isys][ir]->Hresponse();
      hjetresponse_half[isys][ir]->SetName(Form("hjetresponse_half%i", ir));

      hphoresponse_half[isys][ir] = (TH2D*)pho_response_half[isys][ir]->Hresponse();
      hphoresponse_half[isys][ir]->SetName(Form("hphoresponse_half%i",ir));

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
    const char * wfilename = Form("%s/hists/%s_%s_%s_unfolding.root", ana::dir(),trigger.c_str(),sim.c_str(),systag.c_str());
    if (!isMC) wfilename = Form("%s/hists/%s_%s_unfolding.root", ana::dir(),trigger.c_str(),systag.c_str());
    cout << "Writing files to " << wfilename << endl;
    TFile * fout = TFile::Open(wfilename, "RECREATE");

    savehists(hphodr[isys].data(),ana::nJetR);
    savehists(hjetdr[isys].data(),ana::nJetR);

    savehists(hrecojetpt[isys].data(),ana::nJetR);
    savehists(htruthjetpt[isys].data(),ana::nJetR);
    savehists(hjetresponse[isys].data(),ana::nJetR);
    savehists(hrecophopt[isys].data(),ana::nJetR);
    savehists(htruthphopt[isys].data(),ana::nJetR);
    savehists(hphoresponse[isys].data(),ana::nJetR);

    savehists(hrecojetpt_half[isys].data(),ana::nJetR);
    savehists(htruthjetpt_half[isys].data(),ana::nJetR);
    savehists(hjetresponse_half[isys].data(),ana::nJetR);
    savehists(hrecophopt_half[isys].data(),ana::nJetR);
    savehists(htruthphopt_half[isys].data(),ana::nJetR);
    savehists(hphoresponse_half[isys].data(),ana::nJetR);

    savehists(hrecoxj[isys].data(),ana::nJetR);
    savehists(htruthxj[isys].data(),ana::nJetR);
    savehists(hxjresponse[isys].data(),ana::nJetR);

    savehists(hrecoxj_half[isys].data(),ana::nJetR);
    savehists(htruthxj_half[isys].data(),ana::nJetR);
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
        hclusterpt_abcd_truthmatched[isys][i][j]->Write();
      }
    }

    for (int ir = 0; ir < ana::nJetR; ir++) {
      hphoeff[isys][ir] = new TEfficiency(*hphoeffnum[isys][ir], *hphoeffden[isys][ir]);
      hphopur[isys][ir] = new TEfficiency(*hphopurnum[isys][ir], *hphopurden[isys][ir]);
      hjeteff[isys][ir] = new TEfficiency(*hjeteffnum[isys][ir], *hjeteffden[isys][ir]);
      hjetpur[isys][ir] = new TEfficiency(*hjetpurnum[isys][ir], *hjetpurden[isys][ir]);
      hpaireff[isys][ir] = new TEfficiency(*hpaireffnum[isys][ir], *hpaireffden[isys][ir]);
      hpairpur[isys][ir] = new TEfficiency(*hpairpurnum[isys][ir], *hpairpurden[isys][ir]);

      // Weighted entries: Clopper-Pearson does not apply, so use kFNormal with weighted-event errors.
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
