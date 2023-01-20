void AnaSiPM(const Int_t proc)
{
  typedef pair<Double_t, Double_t> eBeam_t;
  const vector<eBeam_t> v_eBeam{{5,41}, {5,100}, {10,100}, {10,275}, {18,275}};
  const Double_t Q2min = 1;
  const Double_t eta_min = 1;
  const Double_t eta_max = 6;

  const char *dir_eic = "/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/zji/data/eic";
  auto f_out = new TFile(Form("%s/histos/energy-%d.root", dir_eic, proc), "RECREATE");

  for(auto eBeam : v_eBeam)
  {
    f_out->cd();
    auto h3_ekin = new TH3F(Form("h3_ekin_%gx%g", eBeam.first, eBeam.second), "Kinetic energy of MC particles;PID;#theta (degree);eKin (GeV)", 4,-0.5,3.5, 40,0.,40., 100,0.,100.);
    auto h3_ehit = new TH3F(Form("h3_ehit_%gx%g", eBeam.first, eBeam.second), "Energy deposit in pECal;z (mm);#theta (degree);log10(ehit (GeV))", 17,3425.,3595., 40,0.,40., 100,-3.,2.);
    auto h2_ehit = new TH2F(Form("h2_ehit_%gx%g", eBeam.first, eBeam.second), "Energy weight in pECal;z (mm);#theta (degree)", 17,3425.,3595., 40,0.,40.);

    TString file_name;
    file_name.Form("%s/endcap/sim_pythia8NCDIS_%gx%g_minQ2_%g-%d.edm4hep.root", dir_eic, eBeam.first, eBeam.second, Q2min, proc);

    TFile *data_file = TFile::Open(file_name);
    if(!data_file)
    {
      cerr << "Cannot open " << file_name << endl;
      continue;
    }
    cout << "Opening " << file_name << endl;

    const Int_t max_track = 10000;
    auto events = (TTree*)data_file->Get("events");
    Int_t pid[max_track] = {}, status[max_track] = {};
    Float_t pmc[3][max_track], ehit[max_track], pos[3][max_track];
    events->SetBranchAddress("MCParticles.PDG", (Int_t*)pid);
    events->SetBranchAddress("MCParticles.generatorStatus", (Int_t*)status);
    events->SetBranchAddress("MCParticles.momentum.x", (Float_t*)pmc[0]);
    events->SetBranchAddress("MCParticles.momentum.y", (Float_t*)pmc[1]);
    events->SetBranchAddress("MCParticles.momentum.z", (Float_t*)pmc[2]);
    events->SetBranchAddress("EcalEndcapPHits.energy", (Float_t*)ehit);
    events->SetBranchAddress("EcalEndcapPHits.position.x", (Float_t*)pos[0]);
    events->SetBranchAddress("EcalEndcapPHits.position.y", (Float_t*)pos[1]);
    events->SetBranchAddress("EcalEndcapPHits.position.z", (Float_t*)pos[2]);

    for(Long64_t iEvent = 0; iEvent < events->GetEntries(); iEvent++)
    {
      events->GetEntry(iEvent);
      for(Int_t iMC = 0; iMC < max_track; iMC++)
      {
        if(status[iMC] == 0) break;
        if(status[iMC] != 1) continue;
        TVector3 v3_pmc(pmc[0][iMC], pmc[1][iMC], pmc[2][iMC]);
        Double_t eta = v3_pmc.Eta();
        if( eta < eta_min || eta > eta_max ) continue;
        Double_t eKin = v3_pmc.Mag();
        Double_t theta = v3_pmc.Theta() * 180. / TMath::Pi();
        Int_t type = -1;
        switch(abs(pid[iMC]))
        {
          case 22:
            type = 0;
            break;
          case 11:
            type = 1;
            break;
          case 211:
            type = 2;
            break;
          default:
            type = 3;
        }
        h3_ekin->Fill((Double_t)type, theta, eKin);
      } // iMC
      for(Long64_t iHit = 0; iHit < events->GetLeaf("EcalEndcapPHits.energy")->GetLen(); iHit++)
        if(ehit[iHit] > 0.)
        {
          TVector3 v3_pos(pos[0][iHit], pos[1][iHit], pos[2][iHit]);
          Double_t theta = v3_pos.Theta() * 180. / TMath::Pi();
          h3_ehit->Fill(pos[2][iHit], theta, TMath::Log10(ehit[iHit]));
          h2_ehit->Fill(pos[2][iHit], theta, ehit[iHit]);
        } // iHit
    } // iEvent

    data_file->Close();
  }

  f_out->Write();
  f_out->Close();
}
