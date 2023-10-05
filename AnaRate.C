void AnaRate(const Int_t proc)
{
  typedef tuple<Double_t, Double_t, Double_t, string> eBeam_t; // eElectron, eProton, minQ2, property
  const vector<eBeam_t> v_eBeam{{5,41,1,"_vtxfix"}, {5,41,10,""}, {5,41,100,""}, {10,100,1,""}, {10,100,10,""}, {10,100,100,""}, {10,100,1000,""}, {10,275,0,"local"}, {10,275,0.5,"local"}, {10,275,1,"local"}, {18,275,1,""}, {18,275,10,""}, {18,275,100,""}, {18,275,1000,""}};
  //const char *version = "23.07.1";
  const char *version = "23.08.0";
  const Float_t edep_th = 0.0015;
  const Float_t eta_min = 3.5;
  const Float_t eta_max = 3.7;

  const char *dir_eic = "/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/zji/data/eic";
  auto f_out = new TFile(Form("%s/histos/ecal-rate-%d.root", dir_eic, proc), "RECREATE");

  set<ULong_t> s_cellID;
  TH1 *h_twr = new TH1F("h_twr", "Number of towers;proc", 1000, -0.5, 999.5);

  for(auto eBeam : v_eBeam)
  {
    f_out->cd();
    TH1 *h_prob = new TH1F(Form("h_prob_%gx%g_%gQ2", get<0>(eBeam), get<1>(eBeam), get<2>(eBeam)), "Probability above threshold", 3, -0.5, 2.5);
    TH1 *h_mul = new TH1F(Form("h_mul_%gx%g_%gQ2", get<0>(eBeam), get<1>(eBeam), get<2>(eBeam)), "Multipicity above threshold", 10, -0.5, 9.5);

    TString file_name;
    if(get<3>(eBeam) == "local")
      file_name.Form("/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/zji/data/eic/endcap/rec_pythia8NCDIS_%gx%g_minQ2_%g-%d.tree.edm4eic.root", get<0>(eBeam), get<1>(eBeam), get<2>(eBeam), proc);
    else
      file_name.Form("s3https://dtn01.sdcc.bnl.gov:9000/eictest/EPIC/RECO/%s/epic_brycecanyon/DIS/NC/%gx%g/minQ2=%g/pythia8NCDIS_%gx%g_minQ2=%g_beamEffects_xAngle=-0.025_hiDiv%s_1.%04d.eicrecon.tree.edm4eic.root", version, get<0>(eBeam), get<1>(eBeam), get<2>(eBeam), get<0>(eBeam), get<1>(eBeam), get<2>(eBeam), get<3>(eBeam).c_str(), proc);

    // Reserve enough space to avoid segmentation fault
    const Int_t max_track = 10000;
    TFile *data_file = TFile::Open(file_name);
    if(!data_file)
    {
      cerr << "Cannot open " << file_name << endl;
      continue;
    }
    cout << "Opening " << file_name << endl;

    auto events = (TTree*)data_file->Get("events");
    // Check the type of variables
    Float_t pos[3][max_track], edep[max_track];
    ULong_t cellID[max_track];
    events->SetBranchAddress("EcalEndcapPRecHits.position.x", (Float_t*)pos[0]);
    events->SetBranchAddress("EcalEndcapPRecHits.position.y", (Float_t*)pos[1]);
    events->SetBranchAddress("EcalEndcapPRecHits.position.z", (Float_t*)pos[2]);
    events->SetBranchAddress("EcalEndcapPRecHits.energy", (Float_t*)edep);
    events->SetBranchAddress("EcalEndcapPRecHits.cellID", (ULong_t*)cellID);

    for(Long64_t i = 0; i < events->GetEntries(); i++)
    {
      events->GetEntry(i);

      h_prob->Fill(0.);
      bool filled = false;
      Int_t mul = 0;

      for(Long64_t j = 0; j < events->GetLeaf("EcalEndcapPRecHits.energy")->GetLen(); j++)
      {
        ROOT::Math::XYZVector v3_pos(pos[0][j], pos[1][j], pos[2][j]);
        Double_t eta = v3_pos.Eta();

        if( eta > eta_min && eta < eta_max )
        {
          if( s_cellID.count(cellID[j]) == 0 )
            s_cellID.insert(cellID[j]);

          if( edep[j] > edep_th )
          {
            mul++;
            h_prob->Fill(1.);
            if( !filled )
            {
              h_prob->Fill(2.);
              filled = true;
            } // filled
          } // edep
        } // eta
      } // leaf

      h_mul->Fill(mul);
    } // event

    data_file->Close();
  } // energy

  h_twr->SetBinContent(proc+1, (Double_t)s_cellID.size());

  f_out->Write();
  f_out->Close();
}
