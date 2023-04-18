void DrawClusters(const char *particle = "gamma")
{
  const vector<Double_t> v_energy{1, 2, 3, 4, 6, 8, 10, 15, 20, 30, 40, 60, 80, 100};
  const Double_t theta_min = 15;
  const Double_t theta_max = 25;

  //const Int_t ntype = 4;
  //const char *clus_name[ntype] = {"RecHits", "TruthClusters", "Clusters", "MergedClusters"};

  //auto f = new TFile(Form("results/clus_%s_theta_%g_%gdeg-res-a316-b15.root", particle, theta_min, theta_max));

  //const Int_t ntype = 2;
  //const char *clus_name[ntype] = {"res-a316-b15", "pecal"};
  //const char *leg_name[ntype] = {"a = 0.00316, b = 0.0015", "a = 0.00340, b = 0.0009"};

  const Int_t ntype = 4;
  const char *clus_name[ntype] = {"full", "full-nosupport", "pecal", "pecal-nosupport"};
  const char *leg_name[ntype] = {"Full detector", "Full detector without Plexiglass & Fr4", "pECal only", "pECal only without Plexiglass & Fr4"};

  TFile *f[ntype];
  for(Int_t it=0; it<ntype; it++)
    f[it] = new TFile(Form("results/clus_%s_theta_%g_%gdeg-%s.root", particle, theta_min, theta_max, clus_name[it]));

  TCanvas *c[ntype+1];
  TGraphErrors *g_res[ntype];
  for(Int_t it = 0; it <= ntype; it++)
  {
    c[it] = new TCanvas(Form("c%d", it), Form("c%d", it), 4*600, 3*600);
    c[it]->Divide(4, 3);
    if(it < ntype)
      g_res[it] = new TGraphErrors(v_energy.size());
  }
  auto c_res = new TCanvas("c_res", "c_res", 600, 600);
  Int_t ipad = 1;

  for(auto energy : v_energy)
  {
    string energy_str(Form("%g%s", energy < 1 ? energy*1e3 : energy, energy < 1 ? "MeV" : "GeV"));

    TH1 *h_edep[ntype];
    for(Int_t it = 0; it < ntype; it++)
      //h_edep[it] = (TH1*)f->Get(Form("h_edep_%s_%s", clus_name[it], energy_str.c_str()));
      h_edep[it] = (TH1*)f[it]->Get(Form("h_edep_%s_%s", "Clusters", energy_str.c_str()));

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
      auto f_gaus = new TF1("f_gaus", "gaus", energy*0.4, energy*1.6);
      f_gaus->SetParameter(0, h_edep[it]->GetMaximum());
      f_gaus->SetParameter(1, energy);
      f_gaus->SetParameter(2, energy*0.03);
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
  }

  c_res->cd();
  gStyle->SetOptFit(0);
  auto leg_res = new TLegend(0.2, 0.6, 0.85, 0.85);
  leg_res->SetBorderSize(0);
  for(Int_t it = 0; it < ntype; it++)
  {
    c[it+1]->Print(Form("results/clusters-res-%s.pdf", clus_name[it]));
    g_res[it]->SetTitle("#sigma/E");
    g_res[it]->GetXaxis()->SetTitle("E [GeV]");
    g_res[it]->GetYaxis()->SetRangeUser(0., 0.13);
    g_res[it]->SetMarkerStyle(it+20);
    g_res[it]->SetMarkerColor(it+1);
    g_res[it]->SetMarkerSize(1.6);
    g_res[it]->Draw(it==0 ? "AP" : "P");

    auto f_res = new TF1("f_res", "TMath::Sqrt([0]*[0]/x+[1]*[1])", 0.5, 101.);
    f_res->SetParameter(0, 0.1);
    f_res->SetParameter(1, 0.03);
    f_res->SetLineWidth(1);
    f_res->SetLineColor(it+1);
    g_res[it]->Fit(f_res, "RQ");
    Double_t a = f_res->GetParameter(0);
    Double_t b = f_res->GetParameter(1);
    leg_res->AddEntry(g_res[it], Form("%s: %.1f%/#sqrt{E} #oplus %.1f%", leg_name[it], a*100., b*100.), "PE");
  }
  leg_res->Draw();

  c[0]->Print("results/clusters-edep.pdf");
  c_res->Print("results/clusters-res.pdf");
}
