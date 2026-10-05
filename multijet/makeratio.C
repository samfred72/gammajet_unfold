void makeratio() {
  const char * names[2] = {"pythia","herwig"};
  // reco-based JER smear variants only - matches analysis.cc's sysNames/fitFuncNames
  const char * sysnames[3] = {"RECO","HIGH","LOW"};
  const char * fitsuffix[3] = {"_JERreco","_JERhigh","_JERlow"};
  const int radii[7] = {2,3,4,5,6,7,8};

  for (int in = 0; in < 2; in++) {
    // one consolidated input file per sim, covering all radii
    TFile * f = TFile::Open(Form("multijet_analysis_%s.root", names[in]), "READ");
    if (!f || f->IsZombie()) {
      std::cout << "Could not open multijet_analysis_" << names[in] << ".root" << std::endl;
      continue;
    }

    for (int ir = 0; ir < 7; ir++) {
      int i = radii[ir];

      TF1 * func[3];
      TH1D * ratio[3];
      TH1D * hdata = (TH1D*)f->Get(Form("leadingJet_r%d", i));
      hdata->Scale(1.0/hdata->Integral("width"));

      for (int j = 0; j < 3; j++) {
        TH1D * hMC = (TH1D*)f->Get(Form("leadingJetPT_Pyth_r%d_%s", i, sysnames[j]));
        hMC->Scale(1.0/hMC->Integral("width"));

        ratio[j] = (TH1D*)hdata->Clone(Form("hratio%i_%s%s", i, names[in], fitsuffix[j]));
        ratio[j]->Reset("ICES");
        ratio[j]->Divide(hdata,hMC);
        ratio[j]->Fit("expo");
        func[j] = (TF1*)ratio[j]->GetListOfFunctions()->FindObject("expo");
        func[j]->SetName(Form("ratio_func%s", fitsuffix[j]));
      }

      TH1D * hdataz = (TH1D*)f->Get(Form("zvtx_data_r%d", i));
      TH1D * hMCz = (TH1D*)f->Get(Form("zvtx_MC_r%d", i));
      hdataz->Scale(1.0/hdataz->Integral(189,213,"width"));
      hMCz->Scale(1.0/hMCz->Integral(189,213,"width"));

      TH1D * ratioz = (TH1D*)hdataz->Clone(Form("hratio_zvtx%i_%s",i,names[in]));
      ratioz->Reset("ICES");
      ratioz->Divide(hdataz,hMCz);

      TFile * wf = TFile::Open(Form("aux/ratio%i_%s.root",i,names[in]),"RECREATE");
      ratioz->Write();
      ratio[0]->Write();
      ratio[1]->Write();
      ratio[2]->Write();
      func[0]->Write();
      func[1]->Write();
      func[2]->Write();
      wf->Close();
    }
    f->Close();
  }
}
