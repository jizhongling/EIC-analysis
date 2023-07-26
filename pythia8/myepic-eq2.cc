#include <iostream>
#include <TString.h>
#include <Pythia8/Pythia.h>
#include <Pythia8Plugins/HepMC3.h>

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

  // Interface for conversion from Pythia8::Event to HepMC event
  HepMC3::Pythia8ToHepMC3 toHepMC;

  // Specify file where HepMC events will be stored
  HepMC3::WriterAscii ascii_io(argv[1]);

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
  // Uncomment and turn on PDF:lepton to enable photon-parton processes
  //pythia.readString("PhotonParton:ggm2bbbar = on");

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
  // Uncomment and turn on PDF:lepton to enable photon-parton processes
  //pythia.readString("PDF:beamB2gamma = on");
  pythia.readString("TimeShower:QEDshowerByL = off");

  // Radiation and hadronization settings
  pythia.readString("PartonLevel:all = on");
  pythia.readString("PartonLevel:MPI = off");
  pythia.readString("PartonLevel:ISR = on");
  pythia.readString("PartonLevel:FSR = on");
  pythia.readString("HadronLevel:all = on");
  pythia.readString("HadronLevel:Decay = on");

  // Kinematic cuts
  const Int_t iparton = 5;  // 5: inclusive; 7: bbar
  const Double_t mom_min = 2.;
  const Double_t mom_max = 100.;
  const Double_t eta_min = 1.4;
  const Double_t eta_max = 4.;

  // Initialize
  pythia.init();

  // Event loop
  for (int iEvent = 0; iEvent < nEvent; ++iEvent)
  {
    if (!pythia.next()) continue;

    if (pythia.event.size() < iparton) continue;
    const Particle &part = pythia.event[iparton];
    Double_t mom = part.pAbs();
    Double_t eta = part.eta();
    if (mom < mom_min || mom > mom_max ||
        eta < eta_min || eta > eta_max)
    {
      iEvent--;
      continue;
    }

    // Construct new empty HepMC event and fill it
    // Default units are (HepMC3::Units::GEV, HepMC3::Units::MM)
    // but can be changed in the GenEvent constructor
    HepMC3::GenEvent hepmcevt;
    toHepMC.fill_next_event(pythia, &hepmcevt);

    // Write the HepMC event to file
    ascii_io.write_event(hepmcevt);
  }

  // Statistics
  pythia.stat();

  return 0;
}
