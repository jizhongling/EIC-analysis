void DrawSiPM()
{
  typedef pair<Double_t, Double_t> eBeam_t;
  const vector<eBeam_t> v_eBeam{{5,41}, {5,100}, {10,100}, {10,275}, {18,275}};
  const Double_t Q2min = 1;

  const char *particle[4] = {"#gamma", "e^{#pm}", "#pi^{#pm}", "others"};
  auto f_seg = TFile::Open("results/energy-SiPM-seg.root");
  auto f_noseg = TFile::Open("results/energy-SiPM-noseg.root");

  auto c0 = new TCanvas("c0", "c0", 4*600, 3*600);
  c0->Divide(4, 3);
  gStyle->SetOptStat(0);

  for(auto eBeam : v_eBeam)
  {
    Int_t ipad = 1;

    auto h3_ekin = (TH3*)f_seg->Get(Form("h3_ekin_%gx%g", eBeam.first, eBeam.second));
    auto h2_ehit_seg = (TH2*)f_seg->Get(Form("h2_ehit_%gx%g", eBeam.first, eBeam.second));
    auto h3_ehit_noseg = (TH3*)f_noseg->Get(Form("h3_ehit_%gx%g", eBeam.first, eBeam.second));

    for(Int_t type = 0; type < 4; type++)
    {
      c0->cd(ipad++);
      h3_ekin->GetXaxis()->SetRange(1+type,1+type);
      h3_ekin->GetZaxis()->SetRangeUser(0.,10.);
      auto h2_ekin_theta = h3_ekin->Project3D("zy");
      h2_ekin_theta->SetTitle(Form("eP: %gx%g GeV, E_{kin}^{MC} for %s", eBeam.first, eBeam.second, particle[type]));
      h2_ekin_theta->DrawCopy("CONT1");
    }

    c0->cd(ipad++);
    auto h_ehit_seg = h2_ehit_seg->ProjectionX("h_ehit_seg");
    h_ehit_seg->SetTitle("#sum E_{hit} vs z");
    h_ehit_seg->DrawCopy("HIST");

    for(Int_t iTheta = 0; iTheta < 6; iTheta++)
    {
      c0->cd(ipad++);
      gPad->SetLogy();
      auto h_ehit_noseg = h3_ehit_noseg->ProjectionZ("h_ehit_noseg", 0,-1, 6+5*iTheta,10+5*iTheta);
      h_ehit_noseg->SetTitle(Form("log_{10}(E_{hit}) for #theta from %d^{#circ} to %d^{#circ}", 5+5*iTheta,10+5*iTheta));
      h_ehit_noseg->DrawCopy("HIST");
    }

    c0->Print(Form("results/energy-SiPM-%gx%g.pdf", eBeam.first, eBeam.second));
    c0->Clear("D");
  }
}