#include "unfolder.h"
using namespace std;

unfolder::~unfolder() {}

bool unfolder::check_pair(jet_object jet, int ir, pho_object pho, bool isreco) {
  float dphi = jet.deltaPhi(pho);
  int iabcd = ana::findabcdBin(pho.iso4, pho.bdt, 0);
  
  int ptbin = ana::findPtBin(pho.pt);
  float val = jet.pt/pho.pt;      
  float lowval = ana::jet_calib_pt_cut[ir]/ana::ptBins[ptbin];
  float lowbin = ana::unfoldXjBins[ana::findUnfoldXjBin(lowval)+1];

  if (ptbin == -1) return false;
  if (val < lowbin) return false;
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
  cout << "running (systag=" << systag << ")..." << endl;

  // systag selects which reco-level variation of the nominal selection/kinematics this
  // pass applies - see the constructor comment in unfolder.h for the full list. Computed
  // once here since it doesn't depend on event or jet radius.
  int systagAbcdBin = (systag == "narrowBDT") ? 1 : (systag == "narrowISO") ? 2 : 0;
  bool systagThreejetVeto = (systag == "threejet");
  // emscale_high/emscale_low shift the EM-calorimeter energy scale by +-1.1%: applied in
  // full to the (entirely-EM) photon cluster, and to only the EM-fraction portion of the
  // jet (the non-EM/hadronic portion of jet_pt_smear[ir] is left untouched). MC-only, same
  // convention as JERhigh/JERlow above - jet_pt_smear[ir] has no Data equivalent, so Data
  // always falls back to its nominal reco pT regardless of systag.
  float systagEmscaleShift = (systag == "emscale_high") ? 0.011 :
                              (systag == "emscale_low")  ? -0.011 : 0.0;
  // JES: Data's reconstructed jet pT has a residual calibration gap relative to MC (found
  // via the in-situ jet-photon pT-balance study - see insitu/), corrected by dividing by
  // 0.9729. jes_high/jes_low vary that correction factor by the study's own +-0.03
  // systematic uncertainty. Data-only - this corrects a gap specific to real Data's
  // calibration, so MC (already on-scale by construction) is untouched regardless of
  // systag, the reverse convention from JERhigh/JERlow/emscale_high/emscale_low above.
  // Applied to recoJetPt (and everything downstream: unfold response matrices, purity
  // histograms, pairing) below, but deliberately NOT to insitu_jet_pt - see rawJetPt.
  float jesCorrection = 0.9729 + ((systag == "jes_high") ? -0.03 :
                                  (systag == "jes_low")  ?  0.03 : 0.0);

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
    // Event selection
    // -----------------------


    vector<bool> keepMC = check_keep_MC(truth_cluster_pt, cluster_pt, truth_jet_pt, jet_pt_smear, trigger);
    if (isMC && !keepMC.at(keepMC.size()-1)) continue;

    // -----------------------
    // Leading photon & isolation
    // -----------------------
    float recoClusterPt = isMC ? cluster_pt * (1.0 + systagEmscaleShift) : cluster_pt;
    pho_object maxpho = pho_object(
        recoClusterPt,
        cluster_e,
        cluster_eta,
        cluster_phi,
        cluster_showershape[8],//10],
        cluster_showershape[9],//11],
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

    for (int ir = 0; ir < ana::nJetR; ir++) {
      if (isMC && !keepMC[ir]) continue;
    
      // JERhigh/JERlow/emscale_high/emscale_low only mean anything for MC (smearing is
      // applied to MC to match Data's resolution - Data has no smeared-high/low variant of
      // itself); jes_high/jes_low go the other way (Data-only - MC is already on-scale by
      // construction, so jesCorrection never enters the isMC branch below).
      float recoJetPt;
      if (isMC) {
        recoJetPt = (systag == "JERhigh") ? jet_pt_smear_high[ir] :
                    (systag == "JERlow")  ? jet_pt_smear_low[ir]  :
                    // scale only the EM-fraction portion of the jet by (1 + shift); the
                    // non-EM/hadronic portion, jet_pt_smear[ir]*(1-jet_emfrac[ir]), is
                    // untouched.
                    (systagEmscaleShift != 0.0) ?
                        jet_pt_smear[ir]*jet_emfrac[ir]
                          + jet_pt_smear[ir]*jet_emfrac[ir]*systagEmscaleShift
                          + jet_pt_smear[ir]*(1-jet_emfrac[ir]) :
                    jet_pt_smear[ir];
      } else {
        recoJetPt = jet_pt_calib[ir] / jesCorrection;
      }
      // Uncorrected Data jet pT, for the insitu tree only (see insitu_jet_pt fill
      // below) - the insitu study is what jesCorrection (0.9729, see above) is itself
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
      
      hphodr[ir]->Fill(maxpho.deltaR(maxpho_truth));
      hjetdr[ir]->Fill(maxjet[ir].deltaR(maxjet_truth[ir]));
      
      hphopurden[ir]->Fill(maxpho.pt);
      hphoeffden[ir]->Fill(maxpho_truth.pt);
      if (check_match(maxpho, maxpho_truth)) {
        hphopurnum[ir]->Fill(maxpho.pt);
        hphoeffnum[ir]->Fill(maxpho_truth.pt);
        // Photon-only quantity, independent of jet radius - fill once (ir==1, same
        // convention as hpurity_num/den below) rather than nJetR times.
        if (ir == 1) {
          hphoIDeff_bdt->Fill(maxpho_truth.pt, maxpho.bdt);
          if (maxpho.iso4 > -999) hphoIDeff_iso->Fill(maxpho_truth.pt, maxpho.iso4);
        }
      }
      hjetpurden[ir]->Fill(maxjet[ir].pt);
      hjeteffden[ir]->Fill(maxjet_truth[ir].pt);
      if (check_match(maxjet[ir], maxjet_truth[ir])) {
        hjetpurnum[ir]->Fill(maxjet[ir].pt);
        hjeteffnum[ir]->Fill(maxjet_truth[ir].pt);
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
      // threejet: reco-only veto (no truth-level third-jet branch exists) - applied on
      // top of the nominal pairing requirement, before it feeds the ABCD fills below.
      if (systagThreejetVeto) ispaired[ir] = ispaired[ir] && !hasthirdjet[ir];
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
      if (ispaired[ir] && iabcd_reco != -1) hrecoxj_abcd[ir][iabcd_reco]->Fill(bin);
      if (ispaired_truth[ir] && iabcd_truth != -1) htruthxj_abcd[ir][iabcd_truth]->Fill(bin);
      if (ispaired[ir] && iabcd_reco != -1) hclusterpt_abcd[ir][iabcd_reco]->Fill(maxpho.pt);

      // In-situ test tree: R=0.4 only, all four ABCD regions, gated on the same strict
      // pairing requirement (check_pair, with the pt-ratio/eta/dphi cuts) used for the
      // response matrix - not yet narrowed to the signal region below.
      if (ir == 2 && ispaired[ir] && iabcd_reco != -1) {
        insitu_pho_pt = maxpho.pt;
        insitu_jet_pt = rawJetPt;
        insitu_abcd = iabcd_reco;
        insitu_tree->Fill();
      }

      ispaired[ir] = ispaired[ir] && iabcd_reco == 0;
      ispaired_truth[ir] = ispaired_truth[ir] && iabcd_truth == 0;

      // -----------------------
      // Fill response matrix
      // -----------------------
      //if (bin < 0 && ispaired[ir]) cout << "ISSUE RECO!!! xj: " << xj << " pt: " << maxpho.pt << endl;
      //if (bin_truth < 0 && ispaired_truth[ir]) cout << "ISSUE TRUTH!!! xj: " << xj_truth << " pt: " << maxpho_truth.pt << endl;
      bool ismatch = ispaired_truth[ir] && ispaired[ir] && check_match(maxjet[ir], maxjet_truth[ir]);
      if (ispaired_truth[ir] && ispaired[ir] && ismatch ) {
        jet_response[ir]->Fill(maxjet[ir].pt, maxjet_truth[ir].pt);
        pho_response[ir]->Fill(maxpho.pt, maxpho_truth.pt);
        jet_response2D[ir]->Fill(bin, bin_truth);

        hphomissfake[ir]->Fill(maxpho.pt,maxpho_truth.pt);
        hjetmissfake[ir]->Fill(maxjet[ir].pt,maxjet_truth[ir].pt);
        hpairmissfake[ir]->Fill(bin,bin_truth);

        if (use_half) {
          jet_response_half[ir]->Fill(maxjet[ir].pt, maxjet_truth[ir].pt);
          pho_response_half[ir]->Fill(maxpho.pt, maxpho_truth.pt);
          jet_response_half2D[ir]->Fill(bin, bin_truth); 
        }
      }
      else if (ispaired_truth[ir] && ispaired[ir] && !ismatch) {
          jet_response[ir]->Miss(maxjet_truth[ir].pt);
          jet_response[ir]->Fake(maxjet[ir].pt);
          pho_response[ir]->Miss(maxpho_truth.pt);
          pho_response[ir]->Fake(maxpho.pt);
          jet_response2D[ir]->Miss(bin_truth);
          jet_response2D[ir]->Fake(bin);
        
          hphomissfake[ir]->Fill(maxpho.pt,100);
          hjetmissfake[ir]->Fill(maxjet[ir].pt,100);
          hpairmissfake[ir]->Fill(bin,ana::nPtBins*(ana::nUnfoldXjBins+2));
          hphomissfake[ir]->Fill(100,maxpho_truth.pt);
          hjetmissfake[ir]->Fill(100,maxjet_truth[ir].pt);
          hpairmissfake[ir]->Fill(ana::nPtBins*(ana::nUnfoldXjBins+2),bin_truth); // +2 for underflow and overflow bins
          
          
          if (use_half) {
            jet_response_half[ir]->Miss(maxjet_truth[ir].pt);
            jet_response_half[ir]->Fake(maxjet[ir].pt);
            pho_response_half[ir]->Miss(maxpho_truth.pt);
            pho_response_half[ir]->Fake(maxpho.pt);
            jet_response_half2D[ir]->Miss(bin_truth);
            jet_response_half2D[ir]->Fake(bin);
          }

      }
      else if (ispaired_truth[ir] && !ispaired[ir]) {
        jet_response[ir]->Miss(maxjet_truth[ir].pt);
        pho_response[ir]->Miss(maxpho_truth.pt);
        jet_response2D[ir]->Miss(bin_truth);
        hphomissfake[ir]->Fill(100,maxpho_truth.pt);
        hjetmissfake[ir]->Fill(100,maxjet_truth[ir].pt);
        hpairmissfake[ir]->Fill(ana::nPtBins*(ana::nUnfoldXjBins+2),bin_truth);
        if (use_half) {
          jet_response_half[ir]->Miss(maxjet_truth[ir].pt);
          pho_response_half[ir]->Miss(maxpho_truth.pt);
          jet_response_half2D[ir]->Miss(bin_truth);
        }
      }
      else if (!ispaired_truth[ir] && ispaired[ir]) {
        jet_response[ir]->Fake(maxjet[ir].pt);
        pho_response[ir]->Fake(maxpho.pt);
        jet_response2D[ir]->Fake(bin);
        hphomissfake[ir]->Fill(maxpho.pt,100);
        hjetmissfake[ir]->Fill(maxjet[ir].pt,100);
        hpairmissfake[ir]->Fill(bin,ana::nPtBins*(ana::nUnfoldXjBins+2));
        if (use_half) {
          jet_response_half[ir]->Fake(maxjet[ir].pt);
          pho_response_half[ir]->Fake(maxpho.pt);
          jet_response_half2D[ir]->Fake(bin);
        }
      }
      
      // ------------------
      // Fill Histograms
      // ------------------
      
      // Purity checking
      //if (ir == 1) { 
      if (ir == 1 && check_pair(maxjet[ir], ir, maxpho, 1)) {
        bool maxpho_is_photon = check_match(maxpho, maxpho_truth);
        int iabcd = ana::findabcdBin(maxpho.iso4, maxpho.bdt, 0);
        float xj = maxjet[ir].pt/maxpho.pt;
        if (iabcd == 0) {
          if (maxpho_is_photon) {
            hpurity_num->Fill(maxpho.pt, xj);
            hpurity_num_1D->Fill(maxpho.pt);
          }
          hpurity_den->Fill(maxpho.pt,xj);
          hpurity_den_1D->Fill(maxpho.pt);
        }
      }
      
     
      int ipt = ana::findPtBin(maxpho.pt); 
      
      if (ispaired[ir]) {
        hpairpurden[ir]->Fill(bin);
        
        hrecojetpt[ir]->Fill(maxjet[ir].pt);
        hrecophopt[ir]->Fill(maxpho.pt);
        hrecoxj[ir]->Fill(bin);
        if (!use_half) {
          hrecojetpt_half[ir]->Fill(maxjet[ir].pt);
          hrecophopt_half[ir]->Fill(maxpho.pt);
          hrecoxj_half[ir]->Fill(bin);
        }
      }

      if (ispaired_truth[ir]) {
        hpaireffden[ir]->Fill(bin_truth);

        htruthjetpt[ir]->Fill(maxjet_truth[ir].pt);
        htruthphopt[ir]->Fill(maxpho_truth.pt);
        htruthxj[ir]->Fill(bin_truth);
        if (!use_half) {
          htruthjetpt_half[ir]->Fill(maxjet_truth[ir].pt);
          htruthphopt_half[ir]->Fill(maxpho_truth.pt);
          htruthxj_half[ir]->Fill(bin_truth);
        }
      }

      if (ismatch) {
        hpairpurnum[ir]->Fill(bin);
        hpaireffnum[ir]->Fill(bin_truth);
      }


      // Drawing event displays
      if (dodraw && ndraw < 100 && ir == 1 && (!ispaired_truth[ir] && ispaired[ir])) {
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

void unfolder::unfold() {
  // -----------------------
  // Unfold
  // -----------------------

  for (int ir = 0; ir < ana::nJetR; ir++) {
    // Full Closure
    RooUnfoldBayes jet_unfold(jet_response[ir], hrecojetpt[ir], niterate, 0, 1);
    jet_unfold.SetVerbose(-1);
    hunfoldjetpt[ir] = (TH1D*)jet_unfold.Hreco();
    hunfoldjetpt[ir]->SetName(Form("hunfoldjetpt%i", ir));
    hjetresponse[ir] = (TH2D*)jet_response[ir]->Hresponse();
    hjetresponse[ir]->SetName(Form("hjetresponse%i", ir));
  
    RooUnfoldBayes pho_unfold(pho_response[ir], hrecophopt[ir], niterate, 0, 1);
    pho_unfold.SetVerbose(-1);
    hunfoldphopt[ir] = (TH1D*)pho_unfold.Hreco();
    hunfoldphopt[ir]->SetName(Form("hunfoldphopt%i",ir));
    hphoresponse[ir] = (TH2D*)pho_response[ir]->Hresponse();
    hphoresponse[ir]->SetName(Form("hphoresponse%i",ir));
    
    RooUnfoldBayes jet_unfold2D(jet_response2D[ir], hrecoxj[ir], niterate, 0, 1);
    jet_unfold2D.SetVerbose(-1);
    hunfoldxj[ir] = (TH1D*)jet_unfold2D.Hreco();
    hunfoldxj[ir]->SetName(Form("hunfoldxj%i",  ir));
    hxjresponse[ir] = (TH2D*)jet_response2D[ir]->Hresponse();
    hxjresponse[ir]->SetName(Form("hxjresponse%i", ir));
    
    // Half closure
    RooUnfoldBayes jet_unfold_half(jet_response_half[ir], hrecojetpt_half[ir], niterate, 0, 1);
    jet_unfold_half.SetVerbose(-1);
    hunfoldjetpt_half[ir] = (TH1D*)jet_unfold_half.Hreco();
    hunfoldjetpt_half[ir]->SetName(Form("hunfoldjetpt_half%i", ir));
    hjetresponse_half[ir] = (TH2D*)jet_response_half[ir]->Hresponse();
    hjetresponse_half[ir]->SetName(Form("hjetresponse_half%i", ir));
  
    RooUnfoldBayes pho_unfold_half(pho_response_half[ir], hrecophopt_half[ir], niterate, 0, 1);
    pho_unfold_half.SetVerbose(-1);
    hunfoldphopt_half[ir] = (TH1D*)pho_unfold_half.Hreco();
    hunfoldphopt_half[ir]->SetName(Form("hunfoldphopt_half%i",ir));
    hphoresponse_half[ir] = (TH2D*)pho_response_half[ir]->Hresponse();
    hphoresponse_half[ir]->SetName(Form("hphoresponse_half%i",ir));
    

    RooUnfoldBayes jet_unfold_half2D(jet_response_half2D[ir], hrecoxj_half[ir], niterate, 0, 1);
    jet_unfold_half2D.SetVerbose(-1);
    hunfoldxj_half[ir] = (TH1D*)jet_unfold_half2D.Hreco();
    hunfoldxj_half[ir]->SetName(Form("hunfoldxj_half%i", ir));
    hxjresponse_half[ir] = (TH2D*)jet_response_half2D[ir]->Hresponse();
    hxjresponse_half[ir]->SetName(Form("hxjresponse_half%i", ir));
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
  insitu_file->cd();
  insitu_tree->Write();
  insitu_file->Close();

  const char * wfilename = Form("/home/samson72/sphnx/gammajet_unfold/hists/%s_%s_%s_unfolding.root",trigger.c_str(),sim.c_str(),systag.c_str());
  if (!isMC) wfilename = Form("/home/samson72/sphnx/gammajet_unfold/hists/%s_%s_unfolding.root",trigger.c_str(),systag.c_str());
  cout << "Writing files to " << wfilename << endl;
  TFile::Open(wfilename, "RECREATE");

  savehists(hphodr,ana::nJetR);
  savehists(hjetdr,ana::nJetR);

  savehists(hrecojetpt,ana::nJetR);
  savehists(htruthjetpt,ana::nJetR);
  savehists(hunfoldjetpt,ana::nJetR);
  savehists(hjetresponse,ana::nJetR);
  savehists(hrecophopt,ana::nJetR);
  savehists(htruthphopt,ana::nJetR);
  savehists(hunfoldphopt,ana::nJetR);
  savehists(hphoresponse,ana::nJetR);

  savehists(hrecojetpt_half,ana::nJetR);
  savehists(htruthjetpt_half,ana::nJetR);
  savehists(hunfoldjetpt_half,ana::nJetR);
  savehists(hjetresponse_half,ana::nJetR);
  savehists(hrecophopt_half,ana::nJetR);
  savehists(htruthphopt_half,ana::nJetR);
  savehists(hunfoldphopt_half,ana::nJetR);
  savehists(hphoresponse_half,ana::nJetR);

  savehists(hrecoxj,ana::nJetR);
  savehists(htruthxj,ana::nJetR);
  savehists(hunfoldxj,ana::nJetR);
  savehists(hxjresponse,ana::nJetR);

  savehists(hrecoxj_half,ana::nJetR);
  savehists(htruthxj_half,ana::nJetR);
  savehists(hunfoldxj_half,ana::nJetR);
  savehists(hxjresponse_half,ana::nJetR);

  savehists(jet_response2D,ana::nJetR);
  savehists(jet_response_half2D,ana::nJetR);

  savehists(hphopurden,ana::nJetR);
  savehists(hphopurnum,ana::nJetR);
  savehists(hphoeffden,ana::nJetR);
  savehists(hphoeffnum,ana::nJetR);
  savehists(hjetpurden,ana::nJetR);
  savehists(hjetpurnum,ana::nJetR);
  savehists(hjeteffden,ana::nJetR);
  savehists(hjeteffnum,ana::nJetR);
  savehists(hpairpurden,ana::nJetR);
  savehists(hpairpurnum,ana::nJetR);
  savehists(hpaireffden,ana::nJetR);
  savehists(hpaireffnum,ana::nJetR);

  for (int i = 0; i < ana::nJetR; i++) {
    for (int j = 0; j < 4; j++) {
      hrecoxj_abcd[i][j]->Write();
      htruthxj_abcd[i][j]->Write();
      hclusterpt_abcd[i][j]->Write();
    }
  }

  for (int ir = 0; ir < ana::nJetR; ir++) {
    hphoeff[ir] = new TEfficiency(*hphoeffnum[ir], *hphoeffden[ir]);
    hphopur[ir] = new TEfficiency(*hphopurnum[ir], *hphopurden[ir]);
    hjeteff[ir] = new TEfficiency(*hjeteffnum[ir], *hjeteffden[ir]);
    hjetpur[ir] = new TEfficiency(*hjetpurnum[ir], *hjetpurden[ir]);
    hpaireff[ir] = new TEfficiency(*hpaireffnum[ir], *hpaireffden[ir]);
    hpairpur[ir] = new TEfficiency(*hpairpurnum[ir], *hpairpurden[ir]);

    hphoeff [ir]->SetName(Form("hphoeff%i",ir));
    hphopur [ir]->SetName(Form("hphopur%i",ir));
    hjeteff [ir]->SetName(Form("hjeteff%i",ir));
    hjetpur [ir]->SetName(Form("hjetpur%i",ir));
    hpaireff[ir]->SetName(Form("hpaireff%i",ir));
    hpairpur[ir]->SetName(Form("hpairpur%i",ir));
  }
  savehists(hphoeff,ana::nJetR);
  savehists(hphopur,ana::nJetR);
  savehists(hjeteff,ana::nJetR);
  savehists(hjetpur,ana::nJetR);
  savehists(hpaireff,ana::nJetR);
  savehists(hpairpur,ana::nJetR);

  savehists(hphomissfake,ana::nJetR);
  savehists(hjetmissfake,ana::nJetR);
  savehists(hpairmissfake,ana::nJetR);

  hpurity_num->Write();
  hpurity_den->Write();
  hpurity_num_1D->Write();
  hpurity_den_1D->Write();

  hphoIDeff_bdt->Write();
  hphoIDeff_iso->Write();
}
