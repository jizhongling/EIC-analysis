void DrawRate()
{
  typedef tuple<Double_t, Double_t, Double_t, Int_t> eBeam_t; // eElectron, eProton, minQ2
  const vector<eBeam_t> v_eBeam{{5,41,1,0}, {5,41,10,0}, {5,41,100,0}, {10,100,1,1}, {10,100,10,1}, {10,100,100,1}, {10,100,1000,1}, {10,275,0,2}, {10,275,0.5,2}, {10,275,1,2}, {18,275,1,3}, {18,275,10,3}, {18,275,100,3}, {18,275,1000,3}};

  const Int_t nb = 4;
  Int_t ig[nb] = {};
  TGraph *g_rate[nb];
  for(Int_t ib=0; ib<nb; ib++)
    g_rate[ib] = new TGraph(4);
  string str_beam[nb];

  auto f = new TFile("results/ecal-rate.root");
  TH1 *h_twr = (TH1*)f->Get("h_twr");
  const Double_t ntwr = h_twr->GetMaximum();
  cout << "Number of involed towers = " << ntwr << endl;

  for(auto eBeam : v_eBeam)
  {
    TH1 *h_prob = (TH1*)f->Get(Form("h_prob_%gx%g_%gQ2", get<0>(eBeam), get<1>(eBeam), get<2>(eBeam)));
    TH1 *h_mul = (TH1*)f->Get(Form("h_mul_%gx%g_%gQ2", get<0>(eBeam), get<1>(eBeam), get<2>(eBeam)));
    h_mul->GetXaxis()->SetRange(2, 10);

    Int_t ib = get<3>(eBeam);
    Double_t Q2 = max(get<2>(eBeam), 0.25);
    Double_t rate1 = h_prob->GetBinContent(2) / h_prob->GetBinContent(1) / ntwr;
    Double_t rate2 = h_prob->GetBinContent(3) / h_prob->GetBinContent(1) / ntwr * h_mul->GetMean();
    if( fabs(rate2/rate1 - 1) > 1e-3 )
      cout << "Rate difference between " << rate1 << " and " << rate2 << endl;
    g_rate[ib]->SetPoint(ig[ib], Q2, rate1);
    ig[ib]++;

    str_beam[ib] = Form("%gx%g", get<0>(eBeam), get<1>(eBeam));
    cout << str_beam[ib] << "," << Q2 << ": " << h_prob->GetBinContent(2) << "/" << h_prob->GetBinContent(1) << " = " << rate1 << ", " << h_mul->GetMean() << endl;
  } // energy

  auto c0 = new TCanvas("c0", "c0", 600, 600);
  auto leg0 = new TLegend(0.1, 0.7, 0.3, 0.9);
  gPad->SetLogx();
  for(Int_t ib=0; ib<nb; ib++)
  {
    g_rate[ib]->SetTitle("Prob with E_{twr} > 15 MeV in #eta = 3.5--3.7");
    g_rate[ib]->GetXaxis()->SetTitle("Q^{2}_{min} (GeV^{2})");
    g_rate[ib]->GetXaxis()->SetLimits(0.1, 2e3);
    g_rate[ib]->GetYaxis()->SetRangeUser(0., 0.08);
    g_rate[ib]->SetMarkerStyle(20+ib);
    g_rate[ib]->SetMarkerColor(1+ib);
    g_rate[ib]->SetMarkerSize(1.6);
    g_rate[ib]->Draw(ib==0?"AP":"P");
    leg0->AddEntry(g_rate[ib], str_beam[ib].c_str(), "P");
  }
  leg0->Draw();
  c0->Print("results/ecal-rate.pdf");
}
