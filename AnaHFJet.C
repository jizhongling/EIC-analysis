void AnaHFJet(const Int_t proc)
{
  typedef pair<Double_t, Double_t> eBeam_t;
  //const vector<eBeam_t> v_eBeam{{5,41}, {5,100}, {10,100}, {10,275}, {18,275}};
  const vector<eBeam_t> v_eBeam{{18,275}};
  const Float_t Q2min = 10;
  const Float_t eta_min = 1.4;
  const Float_t eta_max = 4;
  const Float_t energy_cut = 1;

  const char *dir_eic = "/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/zji/data/eic";
  auto f_out = new TFile(Form("%s/histos/jet-%d.root", dir_eic, proc), "RECREATE");

  for(auto eBeam : v_eBeam)
  {
    f_out->cd();
    TH3 *h3_mass[2];
    TH2 *h2_energy[2];
    TH2 *h2_res[2];
    for(Int_t id = 0; id < 2; id++)
    {
      h3_mass[id] = new TH3F(Form("h3_mass_%s_%gx%g", id==0?"b":"bbar", eBeam.first, eBeam.second),
          "Mass of jets;E_{jet}^{truth} (GeV);m_{jet}^{truth} (GeV);m_{jet}^{reco} (GeV)",
          10,0.,50., 40,0.,10., 40,0.,10.);

      h2_energy[id] = new TH2F(Form("h2_energy_%s_%gx%g", id==0?"b":"bbar", eBeam.first, eBeam.second),
          "Energy of jets;E_{jet}^{truth} (GeV);E_{jet}^{reco} (GeV)",
          50,0.,50., 50,0.,50.);

      h2_res[id] = new TH2F(Form("h2_res_%s_%gx%g", id==0?"b":"bbar", eBeam.first, eBeam.second),
          "Resolution of jets;E_{jet}^{truth} (GeV);resolution",
          20,0.,50., 100,-5.,5.);
    }

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
    const Int_t ntype = 3;
    const char *type_name[ntype] = {"MCParticles", "GeneratedJets", "ReconstructedJets"};
    Int_t pid[ntype][ntrack], status[ntype][ntrack];
    Float_t mom[ntype][3][ntrack], energy[ntype][ntrack], mass[ntype][ntrack];

    auto events = (TTree*)data_file->Get("events");
    for(Int_t it = 0; it < ntype; it++)
    {
      events->SetBranchAddress(Form("%s.PDG", type_name[it]), (Int_t*)pid[it]);
      events->SetBranchAddress(Form("%s.%s", type_name[it], it==0?"generatorStatus":"type"), (Int_t*)status[it]);
      events->SetBranchAddress(Form("%s.momentum.x", type_name[it]), (Float_t*)mom[it][0]);
      events->SetBranchAddress(Form("%s.momentum.y", type_name[it]), (Float_t*)mom[it][1]);
      events->SetBranchAddress(Form("%s.momentum.z", type_name[it]), (Float_t*)mom[it][2]);
      if(it > 0)
      {
        events->SetBranchAddress(Form("%s.energy", type_name[it]), (Float_t*)energy[it]);
        events->SetBranchAddress(Form("%s.mass", type_name[it]), (Float_t*)mass[it]);
      }
    }

    for(Long64_t iEvent = 0; iEvent < events->GetEntries(); iEvent++)
    {
      events->GetEntry(iEvent);

      TVector3 v3_mc[2];
      for(Int_t j = 0; j < events->GetLeaf(Form("%s.PDG", type_name[0]))->GetLen(); j++)
      {
        if(status[0][j] == 0) break;
        if(abs(pid[0][j]) != 5 || status[0][j] != 23) continue;
        Int_t id = (pid[0][j] > 0 ? 0 : 1);
        v3_mc[id] = TVector3(mom[0][0][j], mom[0][1][j], mom[0][2][j]);
      } // leaf

      Float_t energy_bbar[ntype][2] = {};
      Float_t mass_bbar[ntype][2] = {};
      for(Int_t it = 1; it < ntype; it++)
      {
        Float_t min_angle[2] = {1e6, 1e6};
        for(Int_t j = 0; j < events->GetLeaf(Form("%s.PDG", type_name[it]))->GetLen(); j++)
        {
          if(status[it][j] != 0) continue;
          TVector3 v3_jet(mom[it][0][j], mom[1][1][j], mom[1][2][j]);
          for(Int_t id = 0; id < 2; id++)
          {
            Float_t this_angle = v3_jet.Angle(v3_mc[id]);
            if(energy[it][j] > energy_cut && this_angle < min_angle[id])
            {
              energy_bbar[it][id] = energy[it][j];
              mass_bbar[it][id] = mass[it][j];
              min_angle[id] = this_angle;
            }
          } // id
        } // leaf
      } // type

      for(Int_t id = 0; id < 2; id++)
      {
        Float_t res = (energy_bbar[2][id] - energy_bbar[1][id]) / energy_bbar[1][id];
        h3_mass[id]->Fill(energy_bbar[1][id], mass_bbar[1][id], mass_bbar[2][id]);
        h2_energy[id]->Fill(energy_bbar[1][id], energy_bbar[2][id]);
        h2_res[id]->Fill(energy_bbar[1][id], res);
      }
    } // event

    data_file->Close();
  } // beam

  f_out->Write();
  f_out->Close();
}
