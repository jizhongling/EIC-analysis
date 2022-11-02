void DrawClusters()
{
  const Int_t ntype = 4;
  const char *clus_name[ntype] = {"RecHits", "TruthClusters", "Clusters", "MergedClusters"};

  TCanvas *c[ntype+1];
  TGraphErrors *g_res[ntype];
  for(Int_t it = 0; it <= ntype; it++)
  {
    c[it] = new TCanvas(Form("c%d", it), Form("c%d", it), 4*600, 2*600);
    c[it]->Divide(4, 2);
    if(it < ntype)
      g_res[it] = new TGraphErrors(7);
  }
  auto c_res = new TCanvas("c_res", "c_res", 600, 600);
  Int_t ipad = 1;

  for(auto energy : vector<Double_t>{20, 30, 40, 50, 60 ,80, 100})
  {
    const Double_t theta_min = 15;
    const Double_t theta_max = 15;
    const string particle = "gamma";

    TString file_name;
    const char *data_dir = "/gpfs/mnt/gpfs02/phenix/spin/spin1/phnxsp01/zji/data/eic/endcap";
    file_name.Form("%s/rec_%s_%gGeV_theta_%g_%gdeg.tree.edm4eic.root", data_dir, particle.c_str(), energy, theta_min, theta_max);

    TH1 *h_edep[ntype];
    for(Int_t it = 0; it < ntype; it++)
      h_edep[it] = new TH1F(Form("h_edep_%s_%g", clus_name[it], energy), Form("%g GeV; E [GeV]", energy), 100, energy*0.8, energy*1.2);

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
      Float_t edep_sum = 0.;
      for(Int_t it = 0; it < ntype; it++)
        for(Long64_t j = 0; j < events->GetLeaf(Form("EcalEndcapP%s.energy", clus_name[it]))->GetLen(); j++)
        {
          if(it == 0)
            edep_sum += edep[it][j];
          else
            h_edep[it]->Fill(edep[it][j]);
        }
      h_edep[0]->Fill(edep_sum);
    }

    c[0]->cd(ipad);
    gStyle->SetOptStat(0);
    auto leg0 = new TLegend;
    for(Int_t it = 0; it < ntype; it++)
    {
      h_edep[it]->SetLineColor(it + 1);
      h_edep[it]->DrawCopy(it==0 ? "" : "SAME");
      leg0->AddEntry(h_edep[it], clus_name[it], "L");
    }
    leg0->Draw();

    for(Int_t it = 0; it < ntype ; it++)
    {
      c[it+1]->cd(ipad);
      gStyle->SetOptFit(1111);
      auto f_gaus = new TF1("f_gaus", "gaus", energy*0.8, energy*1.2);
      f_gaus->SetParameter(0, h_edep[it]->GetMaximum());
      f_gaus->SetParameter(1, energy);
      f_gaus->SetParameter(2, 0.03*energy);
      f_gaus->SetLineColor(kRed);
      f_gaus->SetLineWidth(1);
      h_edep[it]->Fit(f_gaus, "RQ0");
      h_edep[it]->DrawCopy();
      f_gaus->DrawCopy("SAME");

      Double_t mean = f_gaus->GetParameter(1);
      Double_t sigma = f_gaus->GetParameter(2);
      Double_t emean = f_gaus->GetParError(1);
      Double_t esigma = f_gaus->GetParError(2);
      Double_t res = sigma/mean;
      Double_t eres = res * sqrt(emean*emean/mean/mean + esigma*esigma/sigma/sigma);

      TLatex res_caption;
      res_caption.SetTextFont(62);
      res_caption.SetTextSize(.04);
      res_caption.SetNDC(kTRUE);
      res_caption.DrawLatex(.15, .8, Form("Res = %0.5f", res));

      g_res[it]->SetPoint(ipad-1, energy, res);
      g_res[it]->SetPointError(ipad-1, 0., eres);
    }

    ipad++;
    delete data_file;
  }

  c_res->cd();
  auto leg_res = new TLegend;
  for(Int_t it = 0; it < ntype; it++)
  {
    c[it+1]->Print(Form("results/clusters-res-%s.pdf", clus_name[it]));
    g_res[it]->SetTitle("Resolution");
    g_res[it]->GetXaxis()->SetTitle("E [GeV]");
    g_res[it]->SetMarkerStyle(it+20);
    g_res[it]->SetMarkerColor(it+1);
    g_res[it]->SetMarkerSize(1.6);
    g_res[it]->Draw(it==0 ? "AP" : "P");
    leg_res->AddEntry(g_res[it], clus_name[it], "PE");
  }
  leg_res->Draw();

  c[0]->Print("results/clusters-edep.pdf");
  c_res->Print("results/clusters-res.pdf");
}
