void DrawJetID()
{
  const Int_t ntype = 4;
  const Int_t nw = 4;
  const Int_t npt = 3;

  const char *type_str[ntype] = {"Accuracy (%)", "Efficiency (%)", "Purity (%)", "Bkg rejection"};
  const Float_t ymin[ntype] = { 40, 20,  0,  0};
  const Float_t ymax[ntype] = {100, 80, 20, 15};

  const Float_t weight[nw] = {10, 12, 15, 20};
  const Float_t pt[npt] = {2, 5, 10};

  const Float_t perf[ntype][npt][nw] = {
    {{87, 78, 71, 49}, {87, 77, 72, 54}, {95, 93, 78, 59}},  // Accuracy (%)
    {{22, 39, 50, 73}, {17, 39, 46, 64}, { 0,  6, 32, 57}},  // Efficiency (%)
    {{13, 11, 10,  8}, { 9,  9,  9,  7}, { 4,  9,  7,  6}},  // Purity (%)
    {{10.65, 5.09, 3.56, 1.91}, {11.35, 4.67, 3.79, 2.14}, {506.04, 35.66, 4.96, 2.45}}  // Bkg rejection
  };

  auto c0 = new TCanvas("c0", "c0", 2*600, 2*600);
  c0->Divide(2, 2);
  auto leg0 = new TLegend(0.6, 0.6, 0.9, 0.9);

  for(Int_t it=0; it<ntype; it++)
  {
    mcd(0, it+1);
    for(Int_t iw=0; iw<nw; iw++)
    {
      auto g_perf = new TGraph(npt);
      for(Int_t ipt=0; ipt<npt; ipt++)
        g_perf->SetPoint(ipt, pt[ipt], perf[it][ipt][iw]);
      g_perf->SetTitle(type_str[it]);
      aset(g_perf, "Jet min p_{T} (GeV)","", 0.,12., ymin[it],ymax[it]);
      style(g_perf, iw+20, iw+1);
      g_perf->Draw(iw==0?"AP":"P");
      if(it == 3)
      {
        leg0->AddEntry(g_perf, Form("Sig. weight x%g", weight[iw]), "P");
        leg0->Draw();
      }
    } // weight
  } // type

  c0->Print("results/jet-id.pdf");
}
