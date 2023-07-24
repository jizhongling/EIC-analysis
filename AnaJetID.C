void AnaJetID(const Int_t proc)
{
  typedef pair<Double_t, Double_t> eBeam_t;
  //const vector<eBeam_t> v_eBeam{{5,41}, {5,100}, {10,100}, {10,275}, {18,275}};
  const vector<eBeam_t> v_eBeam{{18,275}};
  const Float_t Q2min = 10;
  const Float_t eta_min = 1.4;
  const Float_t eta_max = 4;
  const Float_t energy_cut = 1;

  const Int_t max_cst = 15;
  Int_t mc_pid;
  Float_t jet_mom[4], jet_mass, cst_mom[4][max_cst], cst_energy[3][max_cst];

  const char *dir_eic = "/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/zji/data/eic";
  auto f_out = new TFile(Form("%s/histos/training-%d.root", dir_eic, proc), "RECREATE");
  auto t_out = new TTree("T", "Jet ID");
  t_out->Branch("mc_pid", &mc_pid, "mc_pid/I");
  t_out->Branch("jet_px", &jet_mom[0], "jet_px/F");
  t_out->Branch("jet_py", &jet_mom[1], "jet_py/F");
  t_out->Branch("jet_pz", &jet_mom[2], "jet_pz/F");
  t_out->Branch("jet_energy", &jet_mom[3], "jet_energy/F");
  t_out->Branch("jet_mass", &jet_mass, "jet_mass/F");
  t_out->Branch("cst_px", (Float_t*)cst_mom[0], Form("cst_px[%d]/F", max_cst));
  t_out->Branch("cst_py", (Float_t*)cst_mom[1], Form("cst_py[%d]/F", max_cst));
  t_out->Branch("cst_pz", (Float_t*)cst_mom[2], Form("cst_pz[%d]/F", max_cst));
  t_out->Branch("cst_track", (Float_t*)cst_mom[0], Form("cst_track[%d]/F", max_cst));
  t_out->Branch("cst_ecal", (Float_t*)cst_mom[1], Form("cst_ecal[%d]/F", max_cst));
  t_out->Branch("cst_hcal", (Float_t*)cst_mom[2], Form("cst_hcal[%d]/F", max_cst));

  for(auto eBeam : v_eBeam)
  {
    TString file_name;
    file_name.Form("%s/endcap/rec_pythia8NCDIS_%gx%g_minQ2_%g-%d.tree.edm4eic.root", dir_eic, eBeam.first, eBeam.second, Q2min, proc);
    TFile *data_file = TFile::Open(file_name);
    if(!data_file)
    {
      cerr << "Cannot open " << file_name << endl;
      continue;
    }
    cout << "Opening " << file_name << endl;

    const Int_t ntrack = 10000;
    const Int_t ntype = 2;
    const char *type_name[ntype] = {"MCParticles", "ReconstructedJets"};
    Int_t pid[ntype][ntrack], status[ntype][ntrack];
    Float_t mom[ntype][4][ntrack], energy_ref[ntype][3][ntrack], mass[ntype][ntrack];

    auto events = (TTree*)data_file->Get("events");
    for(Int_t it = 0; it < ntype; it++)
    {
      events->SetBranchAddress(Form("%s.PDG", type_name[it]), (Int_t*)pid[it]);
      events->SetBranchAddress(Form("%s.%s", type_name[it], it==0?"generatorStatus":"type"), (Int_t*)status[it]);
      events->SetBranchAddress(Form("%s.momentum.x", type_name[it]), (Float_t*)mom[it][0]);
      events->SetBranchAddress(Form("%s.momentum.y", type_name[it]), (Float_t*)mom[it][1]);
      events->SetBranchAddress(Form("%s.momentum.z", type_name[it]), (Float_t*)mom[it][2]);

      if(it == 1)
      {
        events->SetBranchAddress(Form("%s.energy", type_name[it]), (Float_t*)mom[it][3]);
        events->SetBranchAddress(Form("%s.mass", type_name[it]), (Float_t*)mass[it]);
        events->SetBranchAddress(Form("%s.referencePoint.x", type_name[it]), (Float_t*)energy_ref[it][0]);
        events->SetBranchAddress(Form("%s.referencePoint.y", type_name[it]), (Float_t*)energy_ref[it][1]);
        events->SetBranchAddress(Form("%s.referencePoint.z", type_name[it]), (Float_t*)energy_ref[it][2]);
      }
    }

    for(Long64_t iEvent = 0; iEvent < events->GetEntries(); iEvent++)
    {
      events->GetEntry(iEvent);

      TVector3 v3_mc;
      for(Int_t j = 0; j < events->GetLeaf(Form("%s.PDG", type_name[0]))->GetLen(); j++)
      {
        if(status[0][j] == 0) break;
        if(abs(pid[0][j]) > 6 || status[0][j] != 23) continue;
        mc_pid = pid[0][j];
        v3_mc = TVector3(mom[0][0][j], mom[0][1][j], mom[0][2][j]);
        break;
      } // leaf

      for(Int_t it = 1; it < ntype; it++)
      {
        Int_t ic0 = -1;
        Int_t ic1 = -1;
        Int_t ncst = 0;
        Float_t min_angle = 1e6;

        for(Int_t j = 0; j < events->GetLeaf(Form("%s.PDG", type_name[it]))->GetLen(); j++)
        {
          // status = 1 means jet constituents
          if(status[it][j] == 1)
          {
            ncst++;
          } // status = 1

          // status = 0 means jet
          else if(status[it][j] == 0)
          {
            TVector3 v3_jet(mom[it][0][j], mom[it][1][j], mom[it][2][j]);
            Float_t this_angle = v3_jet.Angle(v3_mc);
            if(mom[it][3][j] > energy_cut && this_angle < min_angle)
            {
              jet_mom[0] = mom[it][0][j];
              jet_mom[1] = mom[it][1][j];
              jet_mom[2] = mom[it][2][j];
              jet_mom[3] = mom[it][3][j];
              jet_mass = mass[it][j];
              min_angle = this_angle;
              ic0 = j - ncst;
              ic1 = j;
            } // matched jet
            ncst = 0;
          } // status = 0
        } // leaf

        if(ic0 < 0 || ic1 < 0) continue;

        for(Int_t j = ic0; j < ic1; j++)
        {
          Int_t ic = j - ic0;
          // status = 1 mean jet constituents
          if(status[it][j] == 1 && ic  < max_cst)
          {
            cst_mom[0][ic] = mom[it][0][j];
            cst_mom[1][ic] = mom[it][1][j];
            cst_mom[2][ic] = mom[it][2][j];
            cst_mom[3][ic] = mom[it][3][j];
            cst_energy[0][ic] = energy_ref[it][0][j];
            cst_energy[1][ic] = energy_ref[it][1][j];
            cst_energy[2][ic] = energy_ref[it][2][j];
          } // status = 1
        } // leaf

        t_out->Fill();
      } // type
    } // event

    data_file->Close();
  } // beam

  f_out->Write();
  f_out->Close();
}
