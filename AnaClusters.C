void AnaClusters(const Int_t proc, const char *particle)
{
  const vector<pair<Double_t, Double_t>> v_energy{{1, 1.324}, {2, 1.306}, {4, 1.289}, {8, 1.325}, {16, 1.366}, {24, 1.539}};
  const Double_t theta_min = 23;
  const Double_t theta_max = 37;

  const Int_t ndet = 2;
  const Int_t ntype = 2;
  const char *det_name[ndet+1] = {"EcalEndcapP", "EcalBarrelScFi", "EcalSum"};
  //const char *clus_name[ntype] = {"RecHits", "TruthClusters", "Clusters", "MergedClusters"};
  const char *clus_name[ntype] = {"RecHits", "Clusters"};

  const char *dir_eic = "/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/zji/data/eic";
  auto f_out = new TFile(Form("%s/histos/clus_%s_theta_%g_%gdeg-%d.root", dir_eic, particle, theta_min, theta_max, proc), "RECREATE");

  for(const auto &[energy, scale] : v_energy)
  {
    string energy_str(Form("%g%s", energy < 1 ? energy*1e3 : energy, energy < 1 ? "MeV" : "GeV"));

    f_out->cd();
    TH2 *h2_edep[ndet+1][ntype];
    for(Int_t id = 0; id < ndet+1; id++)
      for(Int_t it = 0; it < ntype; it++)
        if(id < ndet || it == 0)
          h2_edep[id][it] = new TH2F(
              Form("h2_edep_%s%s_%s", det_name[id], clus_name[it], energy_str.c_str()),
              Form("%s;E [GeV];#theta", energy_str.c_str()),
              200, energy*0., energy*2.,
              static_cast<Int_t>(theta_max - theta_min), theta_min, theta_max);

    TString file_name;
    file_name.Form("%s/endcap/rec_%s_%gGeV_theta_%g_%gdeg-%d.tree.edm4eic.root", dir_eic, particle, energy, theta_min, theta_max, proc);
    //file_name.Form("s3https://dtn01.sdcc.bnl.gov:9000/eictest/EPIC/RECO/22.11.0/epic_brycecanyon/SINGLE/gamma/%s/3to50deg/%s_%s_3to50deg.%04d.eicrecon.tree.edm4eic.root", energy_str.c_str(), particle, energy_str.c_str(), proc+1);

    // Reserve enough space to avoid segmentation fault
    const Int_t max_track = 1000;
    TFile *data_file = TFile::Open(file_name);
    if(!data_file)
    {
      cerr << "Cannot open " << file_name << endl;
      continue;
    }
    cout << "Opening " << file_name << endl;

    auto events = (TTree*)data_file->Get("events");
    // Check the type of variables
    Float_t pmc[3][max_track], edep[ndet][ntype][max_track];
    events->SetBranchAddress("MCParticles.momentum.x", (Float_t*)pmc[0]);
    events->SetBranchAddress("MCParticles.momentum.y", (Float_t*)pmc[1]);
    events->SetBranchAddress("MCParticles.momentum.z", (Float_t*)pmc[2]);
    for(Int_t id = 0; id < ndet; id++)
      for(Int_t it = 0; it < ntype; it++)
        events->SetBranchAddress(Form("%s%s.energy", det_name[id], clus_name[it]), (Float_t*)edep[id][it]);

    for(Long64_t i = 0; i < events->GetEntries(); i++)
    {
      events->GetEntry(i);
      // Check MCParticles for the index of the generated particle
      ROOT::Math::XYZVector v3_pmc(pmc[0][2], pmc[1][2], pmc[2][2]);
      Double_t theta = v3_pmc.Theta() * 180. / TMath::Pi();
      if( theta < theta_min || theta > theta_max )
        continue;
      for(Int_t it = 0; it < ntype; it++)
      {
        Float_t edep_tot = 0.;
        for(Int_t id = 0; id < ndet; id++)
        {
          Float_t edep_sum = 0.;
          for(Long64_t j = 0; j < events->GetLeaf(Form("%s%s.energy", det_name[id], clus_name[it]))->GetLen(); j++)
          {
            //if(strcmp(det_name[id], "EcalBarrelScFi") == 0)
            //  edep[id][it][j] *= scale;
            if(it == 0)
            {
              edep_sum += edep[id][it][j];
              edep_tot += edep[id][it][j];
            }
            else
            {
              h2_edep[id][it]->Fill(edep[id][it][j], theta);
            }
          } // leaf
          if(it == 0)
            h2_edep[id][it]->Fill(edep_sum, theta);
        } // det
        if(it == 0)
          h2_edep[ndet][it]->Fill(edep_tot, theta);
      } // type
    } // event

    data_file->Close();
  } // energy

  f_out->Write();
  f_out->Close();
}
