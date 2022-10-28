void DrawClusters()
{
  auto c0 = new TCanvas("c0", "c0", 4*600, 2*600);
  c0->Divide(4, 2);
  Int_t ipad = 1;

  for(auto energy : vector<Double_t>{20, 30, 40, 50, 60 ,80, 100})
  {
    const Double_t theta_min = 15;
    const Double_t theta_max = 15;
    const string particle = "gamma";

    TString file_name;
    const char *data_dir = "/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/zji/data/eic/endcap";
    file_name.Form("%s/rec_%s_%gGeV_theta_%g_%gdeg.tree.edm4eic.root", data_dir, particle.c_str(), energy, theta_min, theta_max);

    const Int_t ntype = 3;
    const char *clus_name[ntype] = {"TruthClusters", "Clusters", "MergedClusters"};

    TH1 *h_edep[ntype];
    for(Int_t it = 0; it < ntype; it++)
      h_edep[it] = new TH1F(Form("h_edep_%s_%g", clus_name[it], energy), Form("%g GeV; p_{T} (GeV)", energy), 100, energy*0.8, energy*1.2);

    const Int_t max_track = 1000;
    cout << "Opening " << file_name << endl;
    auto data_file = new TFile(file_name);
    auto events = (TTree*)data_file->Get("events");
    Float_t edep[ntype][max_track];
    for(Int_t it = 0; it < ntype; it++)
      events->SetBranchAddress(Form("EcalEndcapP%s.energy", clus_name[it]), (Float_t*)edep[it]);

    for(Long64_t i = 0; i < events->GetEntries(); i++)
    {
      events->GetEntry(i);
      for(Int_t it = 0; it < ntype; it++)
        for(Long64_t j = 0; j < events->GetLeaf(Form("EcalEndcapP%s.energy", clus_name[it]))->GetLen(); j++)
          h_edep[it]->Fill(edep[it][j]);
    }

    c0->cd(ipad++);
    auto leg0 = new TLegend;
    for(Int_t it = 0; it < ntype; it++)
    {
      h_edep[it]->SetLineColor(it + 1);
      h_edep[it]->Draw(it==0 ? "" : "SAME");
      leg0->AddEntry(h_edep[it], clus_name[it], "L");
    }
    leg0->Draw();

    delete data_file;
  }

  c0->Print("results/clusters-edep.pdf");
}
