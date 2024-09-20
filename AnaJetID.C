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
  t_out->Branch("cst_energy", (Float_t*)cst_mom[3], Form("cst_energy[%d]/F", max_cst));
  t_out->Branch("cst_track", (Float_t*)cst_energy[0], Form("cst_track[%d]/F", max_cst));
  t_out->Branch("cst_ecal", (Float_t*)cst_energy[1], Form("cst_ecal[%d]/F", max_cst));
  t_out->Branch("cst_hcal", (Float_t*)cst_energy[2], Form("cst_hcal[%d]/F", max_cst));

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

    auto v_mc = new vector<edm4hep::MCParticleData>;
    auto v_jet = new vector<edm4eic::ReconstructedParticleData>;
    auto v_cst = new vector<podio::ObjectID>;
    auto v_rec = new vector<edm4eic::ReconstructedParticleData>;
    auto v_clus = new vector<podio::ObjectID>;
    auto v_ecal = new vector<edm4eic::ClusterData>;

    const UInt_t vsize = 1024;
    v_mc->reserve(vsize);
    v_jet->reserve(vsize);
    v_cst->reserve(vsize);
    v_rec->reserve(vsize);
    v_clus->reserve(vsize);
    v_ecal->reserve(vsize);

    auto events = (TTree*)data_file->Get("events");
    events->SetBranchAddress("MCParticles", &v_mc);
    events->SetBranchAddress("ReconstructedJets", &v_jet);
    events->SetBranchAddress("_ReconstructedJets_particles", &v_cst);
    events->SetBranchAddress("ReconstructedParticles", &v_rec);
    events->SetBranchAddress("_ReconstructedParticles_clusters", &v_clus);
    events->SetBranchAddress("EcalEndcapPClusters", &v_ecal);

    for(Long64_t iEvent = 0; iEvent < events->GetEntries(); iEvent++)
    {
      v_mc->clear();
      v_jet->clear();
      v_cst->clear();
      v_rec->clear();
      v_clus->clear();
      v_ecal->clear();
      events->GetEntry(iEvent);

      TVector3 v3_mc;
      for(UInt_t j = 0; j < v_mc->size(); j++)
      {
        Int_t pid = v_mc->at(j).PDG;
        Int_t status = v_mc->at(j).generatorStatus;
        Float_t px = v_mc->at(j).momentum.x;
        Float_t py = v_mc->at(j).momentum.y;
        Float_t pz = v_mc->at(j).momentum.z;
        if(status == 0) break;
        if(abs(pid) > 6 || status != 23) continue;
        mc_pid = pid;
        v3_mc = TVector3(px, py, pz);
        break;
      }

      UInt_t ijet;
      Float_t min_angle = 1e6;
      bool has_jet = false;
      for(UInt_t j = 0; j < v_jet->size(); j++)
      {
        Float_t px = v_jet->at(j).momentum.x;
        Float_t py = v_jet->at(j).momentum.y;
        Float_t pz = v_jet->at(j).momentum.z;
        Float_t en = v_jet->at(j).energy;
        Float_t mass = v_jet->at(j).mass;
        TVector3 v3_jet(px, py, pz);
        Float_t this_angle = v3_jet.Angle(v3_mc);
        if(en > energy_cut && this_angle < min_angle)
        {
          jet_mom[0] = px;
          jet_mom[1] = py;
          jet_mom[2] = pz;
          jet_mom[3] = en;
          jet_mass = mass;
          min_angle = this_angle;
          ijet = j;
          has_jet = true;
        } // matched jet
      }
      if(!has_jet) continue;

      Int_t ic = 0;
      for(UInt_t j = v_jet->at(ijet).particles_begin; j < v_jet->at(ijet).particles_end; j++)
        if(ic < max_cst)
        {
          UInt_t irec = v_cst->at(j).index;
          cst_mom[0][ic] = v_rec->at(irec).momentum.x;
          cst_mom[1][ic] = v_rec->at(irec).momentum.y;
          cst_mom[2][ic] = v_rec->at(irec).momentum.z;
          cst_mom[3][ic] = v_rec->at(irec).energy;
          cst_energy[0][ic] = 0;
          cst_energy[1][ic] = 0;
          cst_energy[2][ic] = 0;
          for(UInt_t k = v_rec->at(irec).clusters_begin; k < v_rec->at(irec).clusters_end; k++)
          {
            UInt_t iclus = v_clus->at(k).index;
            cst_energy[1][ic] += v_ecal->at(iclus).energy;
          }
          ic++;
        }
      while(ic < max_cst)
      {
        cst_mom[0][ic] = 0;
        cst_mom[1][ic] = 0;
        cst_mom[2][ic] = 0;
        cst_mom[3][ic] = 0;
        cst_energy[0][ic] = 0;
        cst_energy[1][ic] = 0;
        cst_energy[2][ic] = 0;
        ic++;
      }

      t_out->Fill();
    } // event

    data_file->Close();
  } // beam

  f_out->Write();
  f_out->Close();
}
