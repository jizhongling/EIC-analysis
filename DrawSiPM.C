void DrawSiPM()
{
  const Double_t n_pythia = 200.;  // number of combined Pythia files
  const Double_t lum = 10.;  // EIC integrated luminosity in fb^-1

  // Kinematic cuts from Pythia
  const Double_t mom_bin = 1.;
  const Double_t mom_min = 2.;
  const Double_t mom_max = 200.;

  const Double_t eta_bin = 0.02;
  const Double_t eta_min = 1.;
  const Double_t eta_max = 4.;

  const Int_t mom_nbins = static_cast<Int_t>( (mom_max - mom_min) / mom_bin );
  const Int_t eta_nbins = static_cast<Int_t>( (eta_max - eta_min) / eta_bin );

  const Int_t ntypes = 5;
  const char *particle[ntypes] = {"e^{#pm}", "#gamma", "#pi^{0}", "#pi^{#pm}", "others"};

  typedef pair<Double_t, Double_t> eBeam_t;
  //const vector<eBeam_t> v_eBeam{{5,41}, {5,100}, {10,100}, {10,275}, {18,275}};
  const vector<eBeam_t> v_eBeam{{18,275}};

  auto f_pythia = TFile::Open("results/energy-pythia.root");
  auto f_seg = TFile::Open("results/energy-dd4hep-seg.root");
  auto f_noseg = TFile::Open("results/energy-dd4hep-noseg.root");

  auto c0 = new TCanvas("c0", "c0", 2*600, 2*600);
  c0->Divide(2, 2);
  gStyle->SetOptStat(0);

  auto c1 = new TCanvas("c1", "c1", 4*600, 3*600);
  c1->Divide(4, 3);
  gStyle->SetOptStat(0);

  auto c2 = new TCanvas("c2", "c2", 4*600, 2*600);
  c2->Divide(4, 2);
  gStyle->SetOptStat(0);


  for(auto eBeam : v_eBeam)
  {
    Int_t ipad0 = 1;
    Int_t ipad1 = 1;
    Int_t ipad2 = 1;

    auto h3_xsec = (TH3*)f_pythia->Get(Form("h3_xsec_%gx%g", eBeam.first, eBeam.second));
    auto h2_ehit_seg = (TH2*)f_seg->Get(Form("h2_ehit_%gx%g", eBeam.first, eBeam.second));
    auto h3_ehit_noseg = (TH3*)f_noseg->Get(Form("h3_ehit_%gx%g", eBeam.first, eBeam.second));

    h3_xsec->Scale(lum/n_pythia);

    for(Int_t type = 0; type < ntypes - 1; type++)
    {
      c0->cd(ipad0++);
      gPad->SetLogz();
      h3_xsec->GetXaxis()->SetRange(1+type, 1+type);
      auto h2_ekin_eta = h3_xsec->Project3D("zy");
      h2_ekin_eta->Scale(1./mom_bin/eta_bin);
      h2_ekin_eta->SetTitle(Form("L = 10 fb^{-1}: dN/(dpd#eta) (Yield/GeV/rapidity), %s", particle[type]));
      h2_ekin_eta->DrawCopy("COLZ");
    }
    h3_xsec->GetXaxis()->SetRange(1, ntypes);

    c0->Print(Form("results/energy-pythia-2d-%gx%g.pdf", eBeam.first, eBeam.second));
    c0->Clear("D");

    const Double_t deta = 1.;
    for(Double_t eta_low = eta_min; eta_low + deta/2. < eta_max; eta_low += deta)
      for(Int_t type = 0; type < ntypes - 1; type++)
      {
        c1->cd(ipad1++);
        gPad->SetLogy();
        Int_t ieta_low = h3_xsec->GetYaxis()->FindBin(eta_low);
        Int_t ieta_up = h3_xsec->GetYaxis()->FindBin(eta_low + deta) - 1;
        auto h_ekin = h3_xsec->ProjectionZ("h_ekin", 1+type,1+type, ieta_low,ieta_up);
        h_ekin->Scale(1./mom_bin/deta);
        h_ekin->SetTitle(Form("eP: %gx%g GeV, L = 10 fb^{-1}, #eta: %g--%g, %s", eBeam.first,eBeam.second, eta_low,eta_low+deta, particle[type]));
        h_ekin->GetYaxis()->SetTitle("dN/(dpd#eta) (Yield/GeV/rapidity)");
        h_ekin->DrawCopy("HIST");
      }

    c1->Print(Form("results/energy-pythia-%gx%g.pdf", eBeam.first, eBeam.second));
    c1->Clear("D");
    continue;

    c2->cd(ipad2++);
    auto h_ehit_seg = h2_ehit_seg->ProjectionX("h_ehit_seg");
    h_ehit_seg->Scale(100./h_ehit_seg->Integral(0,-1));
    h_ehit_seg->SetTitle("Energy percent in each layer");
    h_ehit_seg->GetYaxis()->SetTitle("Energy percent (%)");
    h_ehit_seg->DrawCopy("HIST");

    c2->cd(ipad2++);
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
      c2->cd(ipad2++);
      gPad->SetLogy();
      auto h_ehit_noseg = h3_ehit_noseg->ProjectionZ("h_ehit_noseg", 0,-1, 6+5*iTheta,10+5*iTheta);
      h_ehit_noseg->SetTitle(Form("log_{10}(E_{hit}) for #theta from %d^{#circ} to %d^{#circ}", 5+5*iTheta,10+5*iTheta));
      h_ehit_noseg->DrawCopy("HIST");
    }

    c2->Print(Form("results/energy-dd4hep-%gx%g.pdf", eBeam.first, eBeam.second));
    c2->Clear("D");
  }
}
