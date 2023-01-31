void DrawSiPM()
{
  typedef pair<Double_t, Double_t> eBeam_t;
  //const vector<eBeam_t> v_eBeam{{5,41}, {5,100}, {10,100}, {10,275}, {18,275}};
  const vector<eBeam_t> v_eBeam{{18,275}};

  const double n_pythia = 100.;
  const char *particle[4] = {"#gamma", "e^{#pm}", "#pi^{#pm}", "others"};
  auto f_pythia = TFile::Open("results/energy-SiPM-pythia8.root");
  auto f_seg = TFile::Open("results/energy-SiPM-seg.root");
  auto f_noseg = TFile::Open("results/energy-SiPM-noseg.root");

  auto c0 = new TCanvas("c0", "c0", 4*600, 4*600);
  c0->Divide(4, 4);
  gStyle->SetOptStat(0);

  for(auto eBeam : v_eBeam)
  {
    Int_t ipad = 1;

    auto h3_xsec = (TH3*)f_pythia->Get(Form("h3_xsec_%gx%g", eBeam.first, eBeam.second));
    auto h2_ehit_seg = (TH2*)f_seg->Get(Form("h2_ehit_%gx%g", eBeam.first, eBeam.second));
    auto h3_ehit_noseg = (TH3*)f_noseg->Get(Form("h3_ehit_%gx%g", eBeam.first, eBeam.second));

    h3_xsec->Scale(1./n_pythia);

    for(Int_t type = 0; type < 4; type++)
    {
      c0->cd(ipad++);
      gPad->SetLogy();
      auto h_ekin = h3_xsec->ProjectionZ("h_ekin", 1+type,1+type, 16,35);
      h_ekin->Scale(1./0.4);
      h_ekin->SetTitle(Form("eP: %gx%g GeV, #eta: 1.3--1.7, %s", eBeam.first, eBeam.second, particle[type]));
      h_ekin->GetYaxis()->SetTitle("d#sigma/(dEd#eta) (fb/GB)");
      h_ekin->GetXaxis()->SetRangeUser(2.,100.);
      h_ekin->DrawCopy("HIST");
    }

    for(Int_t type = 0; type < 4; type++)
    {
      c0->cd(ipad++);
      gPad->SetLogz();
      h3_xsec->GetXaxis()->SetRange(1+type,1+type);
      auto h2_ekin_eta = h3_xsec->Project3D("yz");
      h2_ekin_eta->Scale(1./0.02);
      h2_ekin_eta->SetTitle(Form("COLZ: d#sigma/(dEd#eta) (fb/GB), %s", particle[type]));
      h2_ekin_eta->GetXaxis()->SetRangeUser(2.,100.);
      h2_ekin_eta->DrawCopy("COLZ");
    }

    c0->cd(ipad++);
    auto h_ehit_seg = h2_ehit_seg->ProjectionX("h_ehit_seg");
    h_ehit_seg->Scale(100./h_ehit_seg->Integral(0,-1));
    h_ehit_seg->SetTitle("Energy percent in each layer");
    h_ehit_seg->GetYaxis()->SetTitle("Energy percent (%)");
    h_ehit_seg->DrawCopy("HIST");

    c0->cd(ipad++);
    auto g_ehit_cum = new TGraphErrors(17);
    Double_t tot, etot;
    tot = h_ehit_seg->IntegralAndError(0,-1, etot);
    for(Int_t iLayer = 0; iLayer < 17; iLayer++)
    {
      Double_t cum, ecum;
      cum = h_ehit_seg->IntegralAndError(1,1+iLayer, ecum);
      Double_t frac = cum / tot;
      Double_t efrac = frac * sqrt(ecum*ecum/cum/cum + etot*etot/tot/tot);
      g_ehit_cum->SetPoint(iLayer, 0.5+iLayer, frac);
      g_ehit_cum->SetPointError(iLayer, 0.5, efrac);
    }
    g_ehit_cum->SetTitle("Cumalative energy vs length");
    g_ehit_cum->GetXaxis()->SetTitle("Length (cm)");
    g_ehit_cum->GetYaxis()->SetTitle("Energy fraction");
    g_ehit_cum->SetMarkerStyle(20);
    g_ehit_cum->SetMarkerColor(2);
    g_ehit_cum->SetMarkerSize(1.2);
    g_ehit_cum->Draw("APC");

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
