void AnaClustersForSiPM(const Int_t proc, const char *particle)
{
  const vector<Double_t> v_energy{65, 70, 80, 90, 100};
  const Double_t theta_min = 15;
  const Double_t theta_max = 35;

  const Int_t ntype = 4;
  const char *clus_name[ntype] = {"RecHits", "TruthClusters", "Clusters", "MergedClusters"};

  const char *dir_eic = "/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/zji/data/eic";
  auto f_out = new TFile(Form("%s/histos/clus_%s_theta_%g_%gdeg-%d.root", dir_eic, particle, theta_min, theta_max, proc), "RECREATE");

  for(auto energy : v_energy)
  {
    string energy_str(Form("%g%s", energy < 1 ? energy*1e3 : energy, energy < 1 ? "MeV" : "GeV"));

    f_out->cd();
    TH2 *h2_edep[ntype];
    for(Int_t it = 0; it < ntype; it++)
      h2_edep[it] = new TH2F(Form("h2_edep_%s_%s", clus_name[it], energy_str.c_str()), Form("%s; E [GeV];#theta", energy_str.c_str()), 100,energy*0.1,energy*1.1, 20,theta_min,theta_max);

    TString file_name;
    file_name.Form("%s/endcap/rec_%s_%s_theta_%g_%gdeg-%d.tree.edm4eic.root", dir_eic, particle, energy_str.c_str(), theta_min, theta_max, proc);

    const Int_t max_track = 1000;
    TFile *data_file = TFile::Open(file_name);
    if(!data_file)
    {
      cerr << "Cannot open " << file_name << endl;
      continue;
    }
    cout << "Opening " << file_name << endl;

    auto events = (TTree*)data_file->Get("events");
    Float_t pmc[3][max_track], edep[ntype][max_track];
    events->SetBranchAddress("MCParticles.momentum.x", (Float_t*)pmc[0]);
    events->SetBranchAddress("MCParticles.momentum.y", (Float_t*)pmc[1]);
    events->SetBranchAddress("MCParticles.momentum.z", (Float_t*)pmc[2]);
    for(Int_t it = 0; it < ntype; it++)
      events->SetBranchAddress(Form("EcalEndcapP%s.energy", clus_name[it]), (Float_t*)edep[it]);

    for(Long64_t i = 0; i < events->GetEntries(); i++)
    {
      events->GetEntry(i);
      TVector3 v3_pmc(pmc[0][2], pmc[1][2], pmc[2][2]);
      Double_t theta = v3_pmc.Theta() * 180. / TMath::Pi();
      if( theta < theta_min || theta > theta_max )
        continue;
      Float_t edep_max = 0.;
      for(Int_t it = 0; it < ntype; it++)
        for(Long64_t j = 0; j < events->GetLeaf(Form("EcalEndcapP%s.energy", clus_name[it]))->GetLen(); j++)
        {
          if(it == 0)
          {
            if(edep[it][j] > edep_max)
              edep_max = edep[it][j];
          }
          else
          {
            h2_edep[it]->Fill(edep[it][j], theta);
          }
        }
      h2_edep[0]->Fill(edep_max, theta);
    }

    data_file->Close();
  }

  f_out->Write();
  f_out->Close();
}
