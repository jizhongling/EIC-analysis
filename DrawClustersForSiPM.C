void DrawClustersForSiPM(const char *particle = "pi0")
{
  const vector<Double_t> v_energy{65, 70, 80, 90, 100};
  const Double_t theta_min = 4;
  const Double_t theta_max = 15;

  const Int_t dtheta = 2;
  const Int_t ntheta = 5;

  const Int_t ntype = 4;
  const char *clus_name[ntype] = {"RecHits", "TruthClusters", "Clusters", "MergedClusters"};

  auto f = new TFile(Form("results/clus_%s_theta_%g_%gdeg.root", particle, theta_min, theta_max));

  TCanvas *c[ntheta+2];
  TGraphErrors *g_rec[ntheta];
  for(Int_t ith = 0; ith < ntheta+2; ith++)
  {
    c[ith] = new TCanvas(Form("c%d", ith), Form("c%d", ith), 3*600, 2*600);
    c[ith]->Divide(3, 2);
    if(ith < ntheta)
      g_rec[ith] = new TGraphErrors(v_energy.size());
  }
  auto c_rec = new TCanvas("c_rec", "c_rec", 600, 600);
  Int_t ipad = 1;

  for(auto energy : v_energy)
  {
    string energy_str(Form("%g%s", energy < 1 ? energy*1e3 : energy, energy < 1 ? "MeV" : "GeV"));

    TH2 *h2_edep[ntype];
    for(Int_t it = 0; it < ntype; it++)
      h2_edep[it] = (TH2*)f->Get(Form("h2_edep_%s_%s", clus_name[it], energy_str.c_str()));

    c[0]->cd(ipad);
    gStyle->SetOptStat(0);
    auto leg0 = new TLegend;
    for(Int_t it = ntype-1; it >= 0; it--)
    {
      auto h_edep = (TH1*)h2_edep[it]->ProjectionX("h_edep")->Clone(Form("h_edep_%d", it));
      h_edep->SetTitle(Form("E_{truth} = %s, #theta = %g^{#circ} -- %g^{#circ}", energy_str.c_str(), theta_min, theta_max));
      h_edep->SetLineColor(it + 1);
      h_edep->DrawCopy(it==ntype-1 ? "HIST" : "HIST SAME");
      leg0->AddEntry(h_edep, clus_name[it], "L");
    }
    leg0->Draw();

    c[1]->cd(ipad);
    h2_edep[0]->SetTitle(Form("E_{truth} = %s", energy_str.c_str()));
    h2_edep[0]->DrawCopy("COLZ");

    for(Int_t ith = 0; ith < ntheta; ith++)
    {
      c[ith+2]->cd(ipad);
      auto h_edep = (TH1*)h2_edep[0]->ProjectionX("h_edep_theta", 1+ith*dtheta, dtheta+ith*dtheta)->Clone(Form("h_edep_theta_%d", ith));
      h_edep->SetTitle(Form("E_{truth} = %s, #theta = %g^{#circ} -- %g^{#circ}", energy_str.c_str(), theta_min+ith*dtheta, theta_min+dtheta+ith*dtheta));
      gStyle->SetOptFit(1111);
      auto f_gaus = new TF1("f_gaus", "gaus", energy*0.1, energy*1.1);
      f_gaus->SetParameter(0, h_edep->GetMaximum());
      f_gaus->SetParameter(1, energy*0.5);
      f_gaus->SetParameter(2, energy*0.2);
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

      g_rec[ith]->SetPoint(ipad-1, energy, mean/energy);
      g_rec[ith]->SetPointError(ipad-1, 0., sigma/energy);
    }

    ipad++;
  }

  c_rec->cd();
  auto leg_rec = new TLegend;
  for(Int_t ith = 0; ith < ntheta; ith++)
  {
    c[ith+2]->Print(Form("results/energy-SiPM-maxhit-res-theta_%g_%gdeg.pdf", theta_min+ith*dtheta, theta_min+dtheta+ith*dtheta));
    g_rec[ith]->SetTitle("E_{rec}/E_{truth}");
    g_rec[ith]->GetXaxis()->SetTitle("E [GeV]");
    g_rec[ith]->GetYaxis()->SetTitle("E_{rec}/E_{truth}");
    g_rec[ith]->GetYaxis()->SetRangeUser(0.1, 0.7);
    g_rec[ith]->SetMarkerStyle(ith+20);
    g_rec[ith]->SetMarkerColor(ith+1);
    g_rec[ith]->SetMarkerSize(1.6);
    g_rec[ith]->Draw(ith==0 ? "AP" : "P");
    leg_rec->AddEntry(g_rec[ith], Form("#theta: %g -- %g", theta_min+ith*dtheta, theta_min+dtheta+ith*dtheta), "PE");
  }
  leg_rec->Draw();

  c[0]->Print(Form("results/energy-SiPM-all-rec-theta_%g_%gdeg.pdf", theta_min, theta_max));
  c[1]->Print(Form("results/energy-SiPM-maxhit-rec2D-theta_%g_%gdeg.pdf", theta_min, theta_max));
  c_rec->Print(Form("results/energy-SiPM-maxhit-ratio-theta_%g_%gdeg.pdf", theta_min, theta_max));
}
