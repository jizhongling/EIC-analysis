/// Framework include files
#include "DDG4/Geant4SensDetAction.h"
#include "DDG4/Geant4FastSimSpot.h"

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {

  /// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
  namespace sim {

    struct TestEnergyCut : public Geant4Filter  {
      /// Energy cut value
      double m_energyCut;
      public:
      /// Constructor.
      TestEnergyCut(Geant4Context* c, const std::string& n);
      /// Standard destructor
      virtual ~TestEnergyCut();
      /// Filter action. Return true if hits should be processed
      virtual bool operator()(const G4Step* step) const  override  final  {
        // DD4hep: DDG4/plugins/Geant4SDActions.cpp
        std::cout << std::setprecision(20) << std::scientific;
        std::cout << "Global time: "
          << "Pre (" << std::setw(24) << step->GetPreStepPoint()->GetGlobalTime() << ") "
          << "Post (" << std::setw(24) << step->GetPostStepPoint()->GetGlobalTime() << ") "
          << std::endl;
        std::cout << "Local time: "
          << " Pre (" <<std::setw(24) << step->GetPreStepPoint()->GetLocalTime()  << ") "
          << " Post (" <<std::setw(24) << step->GetPostStepPoint()->GetLocalTime() << ") "
          << std::endl;
        return step->GetTotalEnergyDeposit() > m_energyCut;
      }
      /// GFlash/FastSim interface: Filter action. Return true if hits should be processed
      virtual bool operator()(const Geant4FastSimSpot* spot) const  override  final  {
        return spot->energy() > m_energyCut;
      }
    };
  }
}

/// Framework include files
#include "DD4hep/InstanceCount.h"
#include "DDG4/Factories.h"

// Geant4 include files
#include "G4ParticleTable.hh"
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
#include "G4Track.hh"
#include "G4Step.hh"

using namespace dd4hep::sim;
using namespace dd4hep;
using namespace std;

DECLARE_GEANT4ACTION(TestEnergyCut)

/// Constructor.
TestEnergyCut::TestEnergyCut(Geant4Context* c, const std::string& n)
  : Geant4Filter(c,n) {
    InstanceCount::increment(this);
    declareProperty("Cut",m_energyCut=0.0);
  }

/// Standard destructor
TestEnergyCut::~TestEnergyCut() {
  InstanceCount::decrement(this);
}
