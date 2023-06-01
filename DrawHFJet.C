void DrawHFJet()
{
  typedef pair<Double_t, Double_t> eBeam_t;
  //const vector<eBeam_t> v_eBeam{{5,41}, {5,100}, {10,100}, {10,275}, {18,275}};
  const vector<eBeam_t> v_eBeam{{18,275}};

  const Int_t ntype = 2;
  const char *clus_name[ntype] = {"reco", "charged"};
  const char *leg_name[ntype] = {"Reconstructed particles", "Charged particles"};

  TFile *f[ntype];
  for(Int_t it=0; it<ntype; it++)
    f[it] = new TFile(Form("results/jet-hf-%s.root", clus_name[it]));

  auto c0 = new TCanvas("c0", "c0", 4*600, 3*600);
  c0->Divide(4, 3);

  c0->Print("results/jet-hf.pdf[");
  for(auto eBeam : v_eBeam)
  {
    for(Int_t it = 0; it < ntype; it++)
      for(Int_t id = 0; id < 2; id++)
      {
        //"Mass of jets;E_{jet}^{truth} (GeV);m_{jet}^{truth} (GeV);m_{jet}^{reco} (GeV)",
        //10,0.,50., 40,0.,10., 40,0.,10.);
        auto h3_mass = (TH3*)f[it]->Get(Form("h3_mass_%s_%gx%g", id==0?"b":"bbar", eBeam.first, eBeam.second));
        for(Int_t ipt = 0; ipt < 10; ipt++)
        {
          c0->cd(ipt+1);
          Double_t pTLow = h3_mass->GetXaxis()->GetBinLowEdge(ipt+1);
          Double_t pTUp = h3_mass->GetXaxis()->GetBinUpEdge(ipt+1);
          h3_mass->GetXaxis()->SetRange(ipt+1, ipt+1);
          auto h2_mass = h3_mass->Project3D("zy");
          h2_mass->SetTitle(Form("%gx%g GeV, %s %s jet mass, p_{T} = %.1f--%.1f GeV", eBeam.first, eBeam.second, clus_name[it], id==0?"b":"bbar", pTLow, pTUp));
          h2_mass->DrawCopy("COLZ");
        }
        c0->Print("results/jet-hf.pdf");
        c0->Clear("D");

        //"Resolution of jets;E_{jet}^{truth} (GeV);resolution",
        //20,0.,50., 100,-5.,5.);
        auto h2_res = (TH2*)f[it]->Get(Form("h2_res_%s_%gx%g", id==0?"b":"bbar", eBeam.first, eBeam.second));
        for(Int_t ipt = 0; ipt < 10; ipt++)
        {
          c0->cd(ipt+1);
          Double_t pTLow = h2_res->GetXaxis()->GetBinLowEdge(ipt*2+1);
          Double_t pTUp = h2_res->GetXaxis()->GetBinUpEdge((ipt+1)*2);
          TH1 *h_res = h2_res->ProjectionY(Form("h_res_ipt%d",ipt), ipt*2+1, (ipt+1)*2);
          h_res->SetTitle(Form("%gx%g GeV, %s %s jet res, p_{T} = %.1f--%.1f GeV", eBeam.first, eBeam.second, clus_name[it], id==0?"b":"bbar", pTLow, pTUp));
          h_res->DrawCopy();
        }

        //"Energy of jets;E_{jet}^{truth} (GeV);E_{jet}^{reco} (GeV)",
        //50,0.,50., 50,0.,50.);
        c0->cd(11);
        auto h2_energy = (TH2*)f[it]->Get(Form("h2_energy_%s_%gx%g", id==0?"b":"bbar", eBeam.first, eBeam.second));
        h2_energy->SetTitle(Form("%gx%g GeV, %s %s jet energy", eBeam.first, eBeam.second, clus_name[it], id==0?"b":"bbar"));
        h2_energy->DrawCopy("COLZ");
        c0->Print("results/jet-hf.pdf");
        c0->Clear("D");
      }
  } // beam
  c0->Print("results/jet-hf.pdf]");
}
