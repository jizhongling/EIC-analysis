#include <iostream>
#include <TString.h>
#include <TFile.h>
#include <TH3.h>
#include <Pythia8/Pythia.h>

using namespace std;
using namespace Pythia8;

int get_seed()
{
  int fSeed;

  ifstream devrandom;
  devrandom.open("/dev/urandom", ios::binary);
  devrandom.read((char*)&fSeed, sizeof(fSeed));
  devrandom.close();
  if (fSeed != -1)
  {
    cout << "Got seed from /dev/urandom" << endl;
    fSeed = abs(fSeed)%900000000; // pythia take seed from -1 to 900000000
  }
  else fSeed = 0;
  cout << "seed is " << fSeed << endl;
  return fSeed;
}

int main(int argc, char *argv[])
{
  // Check that correct number of command-line arguments
  if (argc != 6)
  {
    cerr << "Unexpected number of command-line arguments.\n"
      << "Usage: " << argv[0] << " <outfile> <eProton> <eElectron> <Q2min> <nEvent>\n"
      << "Program stopped!" << endl;
    return 1;
  }

  // Beam energies, minimal Q2, number of events to generate
  string eProton(argv[2]);
  string eElectron(argv[3]);
  string Q2min(argv[4]);
  int nEvent = atoi(argv[5]);

  // Generator
  Pythia pythia;

  // Print out frequency
  pythia.readString("Next:numberCount = 50000");

  // Pick new random number seed for each run, based on /dev/urandom
  string seed = to_string(get_seed());
  cout << "string seed " << seed << endl;
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + seed);

  // Set up incoming beams, for frame with unequal beam energies
  pythia.readString("Beams:frameType = 2");
  // BeamA = proton
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:eA = " + eProton);
  // BeamB = electron
  pythia.readString("Beams:idB = 11");
  pythia.readString("Beams:eB = " + eElectron);

  // Set up DIS process within some phase space
  // Neutral current (with gamma/Z interference)
  pythia.readString("WeakBosonExchange:ff2ff(t:gmZ) = on");
  // Uncomment to allow charged current
  //pythia.readString("WeakBosonExchange:ff2ff(t:W) = on");
  pythia.readString("PromptPhoton:all = off");

  // Phase-space cut: minimal Q2 of process
  pythia.readString("PhaseSpace:Q2Min = " + Q2min);
  pythia.readString("PhaseSpace:mHatMin = 0.1");
  pythia.readString("PhaseSpace:pTHatMinDiverge = 0.5");

  // Set dipole recoil on. Necessary for DIS + shower
  pythia.readString("SpaceShower:dipoleRecoil = on");

  // Allow emissions up to the kinematical limit
  // since rate known to match well to matrix elements everywhere
  pythia.readString("SpaceShower:pTmaxMatch = 2");

  // Input PDF set
  pythia.readString("PDF:pSet = 8");

  // QED radiation off lepton not handled yet by the new procedure
  pythia.readString("PDF:lepton = off");
  pythia.readString("TimeShower:QEDshowerByL = off");

  // Radiation and hadronization settings
  pythia.readString("PartonLevel:all = on");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("PartonLevel:ISR = on");
  pythia.readString("PartonLevel:FSR = on");
  pythia.readString("HadronLevel:all = on");
  pythia.readString("HadronLevel:Decay = on");

  // Kinematic cuts
  const Double_t mom_bin = 1.;
  const Double_t mom_min = 2.;
  const Double_t mom_max = 200.;

  const Double_t eta_bin = 0.02;
  const Double_t eta_min = 1.;
  const Double_t eta_max = 4.;

  const Int_t mom_nbins = static_cast<Int_t>( (mom_max - mom_min) / mom_bin );
  const Int_t eta_nbins = static_cast<Int_t>( (eta_max - eta_min) / eta_bin );
  const Int_t ntypes = 5;

  // Output file
  auto f_out = new TFile(argv[1], "RECREATE");
  auto h3_xsec = new TH3F(("h3_xsec_" + eElectron + "x" + eProton).c_str(),
      "Cross section (fb);PID;#eta;p (GeV)",
      ntypes, -0.5, ntypes - 0.5,
      eta_nbins, eta_min, eta_max,
      mom_nbins, mom_min, mom_max);

  // Initialize
  pythia.init();

  // Event loop
  for (int iEvent = 0; iEvent < nEvent; ++iEvent)
  {
    if (!pythia.next()) continue;

    // Particle loop
    for (int iPar = 5; iPar < pythia.event.size(); iPar++)
    {
      const Particle &part = pythia.event[iPar];
      Double_t mom = part.pAbs();
      Double_t eta = part.eta();

      if (!(part.isFinal() || part.id() == 111) ||
          mom < mom_min || mom > mom_max ||
          eta < eta_min || eta > eta_max)
        continue;

      Double_t type = -1;
      switch(part.idAbs())
      {
        case 11:
          type = 0;
          break;
        case 22:
          type = 1;
          break;
        case 111:
          type = 2;
          break;
        case 211:
          type = 3;
          break;
        default:
          type = 4;
      }

      h3_xsec->Fill(type, eta, mom);
    } // End of particle loop
  } // End of event loop

  // Statistics
  pythia.stat();

  // Scale by total cross section (fb)
  Double_t weight = 1e12 * pythia.info.sigmaGen() / pythia.info.nAccepted();
  h3_xsec->Scale(weight);

  // Write to file
  f_out->cd();
  h3_xsec->Write();
  f_out->Close();

  return 0;
}
