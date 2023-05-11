inline Float_t theta2eta(Float_t theta)
{
  return -log( tan( theta / 2. / 180. * M_PI ) );
}

void DrawClustersOverlap(const char *particle = "gamma")
{
  const vector<Double_t> v_energy{1, 2, 8, 16, 24};
  const Double_t theta_min = 23;
  const Double_t theta_max = 37;

  const Int_t ntype = 2;
  const char *clus_name[ntype] = {"image", "noimage"};
  const char *leg_name[ntype] = {"With BarrelImage", "Without BarrelImage"};

  const Int_t nth = 16;
  const Int_t dth = 1;

  TFile *f[ntype];
  for(Int_t it=0; it<ntype; it++)
    f[it] = new TFile(Form("results/clus_%s_theta_%g_%gdeg-%s.root", particle, theta_min, theta_max, clus_name[it]));

  auto c_fit = new TCanvas("c_fit", "c_fit", 4*600, 4*600);
  c_fit->Divide(4, 4);
  auto c_res = new TCanvas("c_res", "c_res", 600, 600);

  Int_t ig_res[ntype][nth] = {};
  TGraphErrors *g_res[ntype][nth];
  for(Int_t it = 0; it < ntype; it++)
    for(Int_t ith = 0; ith < nth; ith++)
      g_res[it][ith] = new TGraphErrors(v_energy.size());

  c_fit->Print("results/clusters-fit-overlap.pdf[");
  for(auto energy : v_energy)
  {
    string energy_str(Form("%g%s", energy < 1 ? energy*1e3 : energy, energy < 1 ? "MeV" : "GeV"));

    for(Int_t it = 0; it < ntype; it++)
    {
      Int_t ipad = 1;
      TH2 *h2_edep = (TH2*)f[it]->Get(Form("h2_edep_%s%s_%s", "EcalSum", "RecHits", energy_str.c_str()));

      for(Int_t ith = 0; ith < nth; ith++)
      {
        c_fit->cd(ipad++);
        gStyle->SetOptFit(1111);
        TH1 *h_edep = h2_edep->ProjectionX(Form("h_edep_%d_%d", it, ith), ith*dth+1, (ith+1)*dth);
        h_edep->SetTitle(Form("%s %s: %.2f--%.2f", leg_name[it], energy_str.c_str(), theta2eta(theta_min+ith*dth), theta2eta(theta_min+(ith+1)*dth)));
        auto f_gaus = new TF1("f_gaus", "gaus", energy*0.4, energy*1.6);
        f_gaus->SetParameter(0, h_edep->GetMaximum());
        f_gaus->SetParameter(1, energy);
        f_gaus->SetParameter(2, energy*0.03);
        f_gaus->SetLineColor(kRed);
        f_gaus->SetLineWidth(1);
        h_edep->Fit(f_gaus, "RQ0");
        h_edep->DrawCopy();
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

        g_res[it][ith]->SetPoint(ig_res[it][ith], energy, res);
        g_res[it][ith]->SetPointError(ig_res[it][ith], 0., eres);
        ig_res[it][ith]++;
      }
      c_fit->Print("results/clusters-fit-overlap.pdf");
      c_fit->Clear("D");
    }
  }
  c_fit->Print("results/clusters-fit-overlap.pdf]");

  c_res->cd();
  gStyle->SetOptFit(0);
  auto leg_res = new TLegend(0.4, 0.4, 0.85, 0.85);
  leg_res->SetBorderSize(0);
  bool first = true;
  for(Int_t it = 1; it < ntype; it++)
    for(Int_t ith = 0; ith < nth-4; ith++)
    {
      g_res[it][ith]->SetTitle("#sigma/E");
      g_res[it][ith]->GetXaxis()->SetTitle("E [GeV]");
      g_res[it][ith]->GetYaxis()->SetRangeUser(0., 0.2);
      g_res[it][ith]->SetMarkerStyle(24);
      g_res[it][ith]->SetMarkerColor(1+ith);
      g_res[it][ith]->SetMarkerSize(1.6);
      g_res[it][ith]->Draw(first ? "AP" : "P");
      first = false;

      auto f_res = new TF1("f_res", "TMath::Sqrt([0]*[0]/x+[1]*[1])", 0.5, 101.);
      f_res->SetParameter(0, 0.1);
      f_res->SetParameter(1, 0.03);
      f_res->SetLineWidth(1);
      f_res->SetLineColor(kRed);
      g_res[it][ith]->Fit(f_res, "RQ");
      Double_t a = f_res->GetParameter(0);
      Double_t b = f_res->GetParameter(1);
      leg_res->AddEntry(g_res[it][ith], Form("%.2f--%.2f: %.1f%/#sqrt{E} #oplus %.1f%", theta2eta(theta_min+ith*dth), theta2eta(theta_min+(ith+1)*dth), a*100., b*100.), "P");
    }
  leg_res->Draw();
  c_res->Print("results/clusters-res-overlap.pdf");
}
