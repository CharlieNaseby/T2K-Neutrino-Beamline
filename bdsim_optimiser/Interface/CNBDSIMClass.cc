/* 
Beam Delivery Simulation (BDSIM) Copyright (C) Royal Holloway, 
University of London 2001 - 2024.

This file is part of BDSIM.

BDSIM is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published 
by the Free Software Foundation version 3 of the License.

BDSIM is distributed in the hope that it will be useful, but 
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BDSIM.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "CNBDSIMClass.h"
#include "fastlist.h"

#include "BDSExecOptions.hh"     // executable command line options 
#include "BDSGlobalConstants.hh" //  global parameters

#include <algorithm>
#include <csignal>
#include <cstdlib>
#include <cstdio>

#include "G4EventManager.hh" // Geant4 includes
#include "G4GenericBiasingPhysics.hh"
#include "G4GeometryManager.hh"
#include "G4GeometryTolerance.hh"
#include "G4PhysicsListHelper.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4SteppingManager.hh"
#include "G4TrackingManager.hh"
#include "G4Version.hh"
#include "G4VModularPhysicsList.hh"
#include "G4ParticleGun.hh"
#include "G4StateManager.hh"

#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"  
  
  
#include "CLHEP/Units/SystemOfUnits.h"

#include "BDSAcceleratorModel.hh"
#include "BDSAperturePointsLoader.hh"
#include "BDSBeamPipeFactory.hh"
#include "BDSBunch.hh"
#include "BDSBunchFactory.hh"
#include "BDSCavityFactory.hh"
#include "BDSColours.hh"
#include "BDSComponentFactoryUser.hh"
#include "BDSDebug.hh"
#include "BDSDetectorConstruction.hh"
#include "BDSEventAction.hh"
#include "BDSException.hh"
#include "BDSFieldFactory.hh"
#include "BDSFieldLoader.hh"
#include "BDSGeometryFactory.hh"
#include "BDSGeometryFactorySQL.hh"
#include "BDSGeometryWriter.hh"
#include "BDSIonDefinition.hh"
#include "BDSMaterials.hh"
#include "BDSOutput.hh"
#include "BDSOutputFactory.hh"
#include "BDSParallelWorldUtilities.hh"
#include "BDSParser.hh" // Parser
#include "BDSParticleDefinition.hh"
#include "BDSPhysicsUtilities.hh"
#include "BDSPrimaryGeneratorAction.hh"
#include "BDSRandom.hh" // for random number generator from CLHEP
#include "BDSRunAction.hh"
#include "BDSRunManager.hh"
#include "BDSSamplerRegistry.hh"
#include "BDSSDManager.hh"
#include "BDSSteppingAction.hh"
#include "BDSStackingAction.hh"
#include "BDSTemporaryFiles.hh"
#include "BDSTrackingAction.hh"
#include "BDSUtilities.hh"
#include "BDSVisManager.hh"
#include "BDSWarning.hh"
#include "BDSAcceleratorComponentRegistry.hh"
#include "BDSFieldBuilder.hh"
#include "BDSFieldMagQuadrupole.hh"
#include "BDSFieldMagDipole.hh"
#include "BDSFieldMagDipoleQuadrupole.hh"
#include "BDSFieldMagGlobal.hh"

#include "G4TransportationManager.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"


#include "BDSFieldObjects.hh"
#include "BDSFieldInfo.hh"
#include "BDSFieldMagGlobalPlacement.hh"
#include "BDSIntegratorQuadrupole.hh"
#include "BDSIntegratorDipoleQuadrupole.hh"
#include "BDSOutputROOTEventSampler.hh"
#include "BDSIntegratorDipoleRodrigues2.hh"
//wild optimism
class MyPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
    G4ParticleGun* fParticleGun;

public:
    MyPrimaryGeneratorAction() {
        fParticleGun = new G4ParticleGun(1);  // One particle at a time
        auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("proton");
        fParticleGun->SetParticleDefinition(particleDefinition);
        fParticleGun->SetParticleEnergy(100 * CLHEP::MeV);
        fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, -0.5 * CLHEP::m));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 1));
    }

    ~MyPrimaryGeneratorAction() override { delete fParticleGun; }

    void GeneratePrimaries(G4Event* anEvent) override {
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
};

CNBDSIM::CNBDSIM():
  ignoreSIGINT(false),
  usualPrintOut(true),
  initialised(false),
  initialisationResult(1),
  argcCache(0),
  argvCache(nullptr),
  parser(nullptr),
  bdsOutput(nullptr),
  bdsBunch(nullptr),
  runManager(nullptr),
  userComponentFactory(nullptr),
  userPhysicsList(nullptr),
  realWorld(nullptr)
{;}

CNBDSIM::CNBDSIM(int argc, char** argv, bool usualPrintOutIn):
  ignoreSIGINT(false),
  usualPrintOut(usualPrintOutIn),
  initialised(false),
  initialisationResult(1),
  argcCache(argc),
  argvCache(argv),
  parser(nullptr),
  bdsOutput(nullptr),
  bdsBunch(nullptr),
  runManager(nullptr),
  userComponentFactory(nullptr),
  userPhysicsList(nullptr),
  realWorld(nullptr)
{
  initialisationResult = Initialise();
}

int CNBDSIM::Initialise(int argc, char** argv, bool usualPrintOutIn)
{
  argcCache = argc;
  argvCache = argv;
  usualPrintOut = usualPrintOutIn;
  initialisationResult = Initialise();
  return initialisationResult;
}

int CNBDSIM::Initialise()
{
  /// Initialize executable command line options reader object
  const BDSExecOptions* execOptions = new BDSExecOptions(argcCache,argvCache);
  if (usualPrintOut)
    {execOptions->Print();}
  ignoreSIGINT = execOptions->IgnoreSIGINT(); // different sig catching for cmake

//  execOptions->PrintCopyright();
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "DEBUG mode is on." << G4endl;
#endif

  /// Parse lattice file
  parser = BDSParser::Instance(execOptions->InputFileName());
  /// Update options generated by parser with those from executable options.
  parser->AmalgamateOptions(execOptions->Options());
  parser->AmalgamateBeam(execOptions->Beam(), execOptions->Options().recreate);
  /// Check options for consistency
  parser->CheckOptions();



  /// Explicitly initialise materials to construct required materials before global constants.
  BDSMaterials::Instance()->PrepareRequiredMaterials(execOptions->Options().verbose);

  /// No longer needed. Everything can safely use BDSGlobalConstants from now on.
  delete execOptions;

  /// Force construction of global constants after parser has been initialised (requires
  /// materials first). This uses the options and beam from BDSParser.
  /// Non-const as we'll update the particle definition.
  //BDSGlobalConstants* 
  globals = BDSGlobalConstants::Instance();

  /// Initialize random number generator
  BDSRandom::CreateRandomNumberGenerator(globals->RandomEngine());
  BDSRandom::SetSeed(); // set the seed from options

  /// Construct output
  bdsOutput = BDSOutputFactory::CreateOutput(globals->OutputFormat(),
                                             globals->OutputFileName());

  /// Check geant4 exists in the current environment
  if (!BDS::Geant4EnvironmentIsSet())
    {
      G4cerr << "No Geant4 environmental variables found - please source geant4.sh environment" << G4endl;
      G4cout << "A common fault is the wrong Geant4 environment as compared to the one BDSIM was compiled with." << G4endl;
      return 1;
    }

  /// Construct mandatory run manager (the G4 kernel) and
  /// register mandatory initialization classes.
  runManager = new BDSRunManager();

  /// Register the geometry and parallel world construction methods with run manager.
  realWorld = new BDSDetectorConstruction(userComponentFactory);
  
  /// Here the geometry isn't actually constructed - this is called by the runManager->Initialize()
  auto parallelWorldsRequiringPhysics = BDS::ConstructAndRegisterParallelWorlds(realWorld,
                                                                                realWorld->BuildSamplerWorld(),
                                                                                realWorld->BuildPlacementFieldsWorld());
  runManager->SetUserInitialization(realWorld);

  /// For geometry sampling, phys list must be initialized before detector.
  /// BUT for samplers we use a parallel world and this HAS to be before the physics
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "> Constructing physics processes" << G4endl;
#endif
  G4String physicsListName = parser->GetOptions().physicsList;

#if G4VERSION_NUMBER > 1049
  // from 10.5 onwards they have a looping particle killer that warnings and kills particles
  // deemed to be looping that are <100 MeV. This is unrelated to the primary energy so troublesome.
  // set to the 'low' limits here ~10keV. This must be done before any physics is created as the
  // parameters are copied into the transportation physics process for each particle and it's very
  // hard to sift through and fix afterwards
  G4PhysicsListHelper::GetPhysicsListHelper()->UseLowLooperThresholds();
#endif
  // sampler physics process for parallel world tracking must be instantiated BEFORE
  // regular physics.
  // Note, we purposively don't create a parallel world process for the curvilinear
  // world as we don't need the track information from it - unreliable that way. We
  // query the geometry directly using our BDSAuxiliaryNavigator class.
  auto parallelWorldPhysics = BDS::ConstructParallelWorldPhysics(parallelWorldsRequiringPhysics);
  G4int physicsVerbosity = globals->PhysicsVerbosity();
  physList;
  if (userPhysicsList)
    {
      G4cout << "Using externally registered user defined physics list" << G4endl;
      physList = userPhysicsList;
    }
  else
    {physList = BDS::BuildPhysics(physicsListName, physicsVerbosity);}

  // create geometry sampler and register importance sampling biasing. Has to be here
  // before physicsList is "initialised" in run manager.
  if (BDSGlobalConstants::Instance()->UseImportanceSampling())
    {BDS::RegisterImportanceBiasing(parallelWorldsRequiringPhysics,physList);}

  // Construction of the physics lists defines the necessary particles and therefore
  // we can calculate the beam rigidity for the particle the beam is designed w.r.t. This
  // must happen before the geometry is constructed (which is called by
  // runManager->Initialize()).
  BDSParticleDefinition* designParticle = nullptr;
//  BDSParticleDefinition* beamParticle = nullptr;
  G4bool beamDifferentFromDesignParticle = false;
  BDS::ConstructDesignAndBeamParticle(BDSParser::Instance()->GetBeam(),
                                      globals->FFact(),
                                      designParticle,
                                      beamParticle,
                                      beamDifferentFromDesignParticle);
  G4double minEK = globals->MinimumKineticEnergy();
  if (beamParticle->KineticEnergy() < minEK && BDS::IsFinite(minEK))
    {throw BDSException("option, minimumKineticEnergy is higher than kinetic energy of the beam - all primary particles wil be killed!");}
  if (usualPrintOut)
    {
      G4cout << "Design particle properties: " << G4endl << *designParticle;
      if (beamDifferentFromDesignParticle)
        {G4cout << "Beam particle properties: " << G4endl << *beamParticle;}
    }
  // update rigidity where needed
  realWorld->SetDesignParticle(designParticle);
  BDSFieldFactory::SetDesignParticle(designParticle);
  BDSGeometryFactorySQL::SetDefaultRigidity(designParticle->BRho()); // used for sql field loading
  
  // Muon splitting - optional - should be done *after* biasing to work with it - TBC it's before...
  BDS::BuildMuonBiasing(physList);
  
  BDS::RegisterSamplerPhysics(parallelWorldPhysics, physList);
  auto biasPhysics = BDS::BuildAndAttachBiasWrapper(parser->GetBiasing());
  if (biasPhysics)//could be nullptr and can't be passed to geant4 like this
    {physList->RegisterPhysics(biasPhysics);}
  runManager->SetUserInitialization(physList);

  /// Instantiate the specific type of bunch distribution.
  bdsBunch = BDSBunchFactory::CreateBunch(beamParticle,
                                          parser->GetBeam(),
                                          globals->BeamlineTransform(),
                                          globals->BeamlineS(),
                                          globals->GeneratePrimariesOnly());
  G4cout << "Bunch distribution: \"" << bdsBunch->Name() << "\"" << G4endl;
  /// We no longer need beamParticle so delete it to avoid confusion. The definition is
  /// held inside bdsBunch (can be updated dynamically).
// CERN  delete beamParticle;
  /// Construct extra common particles for possible tracking if required without using a physics list.
  if (bdsBunch->ExpectChangingParticleType())
    {BDS::ConstructExtendedParticleSet();}
  
  /// Optionally generate primaries only and exit
  /// Unfortunately, this has to be here as we can't query the geant4 particle table
  /// until after the physics list has been constructed and attached a run manager.
  if (globals->GeneratePrimariesOnly())
    {
      GeneratePrimariesOnly(globals);
      return 0;
    }
  
  /// Print the geometry tolerance
  G4GeometryTolerance* theGeometryTolerance = G4GeometryTolerance::GetInstance();
  if (usualPrintOut)
    {
      G4cout << __METHOD_NAME__ << "Geometry Tolerances: "     << G4endl;
      G4cout << __METHOD_NAME__ << std::setw(12) << "Surface: " << std::setw(7) << theGeometryTolerance->GetSurfaceTolerance() << " mm"   << G4endl;
      G4cout << __METHOD_NAME__ << std::setw(12) << "Angular: " << std::setw(7) << theGeometryTolerance->GetAngularTolerance() << " rad"  << G4endl;
      G4cout << __METHOD_NAME__ << std::setw(12) << "Radial: "  << std::setw(7) << theGeometryTolerance->GetRadialTolerance()  << " mm"   << G4endl;
    }
  
  /// Set user action classes
  BDSEventAction* eventAction = new BDSEventAction(bdsOutput);
  runManager->SetUserAction(eventAction);
  
  BDSRunAction* runAction = new BDSRunAction(bdsOutput,
                                             bdsBunch,
                                             bdsBunch->ParticleDefinition()->IsAnIon(),
                                             eventAction,
                                             globals->StoreTrajectorySamplerID());
  runManager->SetUserAction(runAction);
  
  // Only add stepping action if it is actually used, so do check here (for performance reasons)
  G4int verboseSteppingEventStart = globals->VerboseSteppingEventStart();
  G4int verboseSteppingEventStop  = BDS::VerboseEventStop(verboseSteppingEventStart,
                                                          globals->VerboseSteppingEventContinueFor());
  if (globals->VerboseSteppingBDSIM())
    {
      runManager->SetUserAction(new BDSSteppingAction(true,
                                                      verboseSteppingEventStart,
                                                      verboseSteppingEventStop));
    }
  
  runManager->SetUserAction(new BDSTrackingAction(globals->Batch(),
                                                  globals->StoreTrajectory(),
                                                  globals->StoreTrajectoryOptions(),
                                                  eventAction,
                                                  verboseSteppingEventStart,
                                                  verboseSteppingEventStop,
                                                  globals->VerboseSteppingPrimaryOnly(),
                                                  globals->VerboseSteppingLevel()));

  runManager->SetUserAction(new BDSStackingAction(globals));
  
  primaryGeneratorAction = new BDSPrimaryGeneratorAction(bdsBunch, parser->GetBeam(), globals->Batch());



  // possibly updated after the primary generator as loaded a beam file
  eventAction->SetPrintModulo(BDSGlobalConstants::Instance()->PrintModuloEvents());
  runManager->SetUserAction(primaryGeneratorAction);
  BDSFieldFactory::SetPrimaryGeneratorAction(primaryGeneratorAction);


  /// Initialize G4 kernel
  runManager->Initialize();

  /// Create importance store for parallel importance world
  if (globals->UseImportanceSampling())
    {BDS::AddIStore(parallelWorldsRequiringPhysics);}

  /// Implement bias operations on all volumes only after G4RunManager::Initialize()
  realWorld->BuildPhysicsBias();

  if (usualPrintOut && globals->PhysicsVerbose())
    {
      BDS::PrintPrimaryParticleProcesses(bdsBunch->ParticleDefinition()->Name());
      BDS::PrintDefinedParticles();
    }

  /// Set verbosity levels at run and G4 event level. Per event and stepping are controlled
  /// in event, tracking and stepping action. These have to be done here due to the order
  /// of construction in Geant4.
  runManager->SetVerboseLevel(std::min(globals->VerboseRunLevel(), globals->PhysicsVerbosity()));
  G4EventManager::GetEventManager()->SetVerboseLevel(globals->VerboseEventLevel());
  G4EventManager::GetEventManager()->GetTrackingManager()->SetVerboseLevel(globals->VerboseTrackingLevel());
  
  /// Close the geometry in preparation for running - everything is now fixed.
  G4bool bCloseGeometry = G4GeometryManager::GetInstance()->CloseGeometry();
  if (!bCloseGeometry)
    {
      G4cerr << __METHOD_NAME__ << "error - geometry not closed." << G4endl;
      return 1;
    }

  if (globals->ExportGeometry())
    {
      BDSGeometryWriter geometrywriter;
      geometrywriter.ExportGeometry(globals->ExportType(),
                                    globals->ExportFileName());
    }

  initialised = true;

  randomState = BDSRandom::GetSeedState();

  return 0;
}


void CNBDSIM::BeamOn(int nGenerate)
{
  BDSRandom::SetSeed(1989); // set the seed at the start of each run
  for(auto tr : bdsOutput->samplerTrees) tr->FlushCache();  //clear the storage for the samplers
  for(auto tr : bdsOutput->samplerTrees) tr->FlushLocal();  //clear the storage for the samplers
//bdsfieldobjects includes both the field and the integrator both use the strength
//to calculate things but we only need to change integrator for our purposes

  for(int i=0; i<BDSAcceleratorModel::Instance()->fields.size(); i++){ //loops over all fields 
    auto *fieldInfo =  BDSAcceleratorModel::Instance()->fields[i]->GetInfo();
    if(fieldInfo->FieldType() == BDS::DetermineFieldType(G4String("quadrupole"))){
      auto *integrator = ((BDSIntegratorQuadrupole*)(BDSAcceleratorModel::Instance()->fields[i]->GetIntegrator()));
//      G4cout << "quadrupole field name " << BDSFieldBuilder::Instance()->lvs[i][0]->GetName() << G4endl;
//      if(integrator->bPrime > 0){
//        integrator->bPrime = fieldInfo->BRho() * 1.23/CLHEP::m2;
//        G4cout <<"setting quad strength to k1 = 1.23 bPrime = " << integrator->bPrime << G4endl;
        (*BDSFieldBuilder::Instance()->infos[i]->MagnetStrength())["k1"] *= 1.01;
//      }
    }
    else if(fieldInfo->FieldType() == BDS::DetermineFieldType(G4String("dipole"))){
//      G4cout << "dipole field " << *fieldInfo->MagnetStrength() << G4endl;
//      G4cout << "dipole field name " << BDSFieldBuilder::Instance()->lvs[i][0]->GetName() << G4endl;
      auto *integrator = ((BDSIntegratorDipoleQuadrupole*)(BDSAcceleratorModel::Instance()->fields[i]->GetIntegrator()));

//      G4cout << "dipole field ratio " << integrator->fieldRatio << G4endl;
      std::string name = BDSFieldBuilder::Instance()->lvs[i][0]->GetName();
//      G4cout << name <<" field strength according to the field obj" << ((BDSFieldMagDipoleQuadrupole*)(BDSAcceleratorModel::Instance()->fields[i]->GetField()))->localField << G4endl;
//      ((BDSFieldMag*)(((BDSFieldMagGlobal*)(BDSAcceleratorModel::Instance()->fields[i]->GetField()))->field))->isDipoleQuadrupole();
      if(name.find("BPD1") != std::string::npos){
        std::cout << "found BPD1" << std::endl;
//        integrator->nominalRho*=0.001;
//        integrator->dipole->fieldScale = 1.1;
//        G4cout << "BPD1 field strength according to the field obj" << ((BDSFieldMagDipoleQuadrupole*)(BDSAcceleratorModel::Instance()->fields[i]->GetField()))->dipole->localField << G4endl;
//        ((BDSFieldMagDipole*)(BDSAcceleratorModel::Instance()->fields[i]->GetField()))->localField *= 1.1;
//        ((BDSFieldMagDipole*)(((BDSFieldMagGlobal*)(BDSAcceleratorModel::Instance()->fields[i]->GetField()))->field))->localField *= 1.1;
        //lets try something wacky
//        auto *objs = (BDSFieldObjects*)(BDSFieldMagGlobal*)(BDSAcceleratorModel::Instance()->fields[i]);

        //objs->replaceEqofM((G4MagneticField*)(objs->GetField()));
        
        (*BDSFieldBuilder::Instance()->infos[i]->MagnetStrength())["field"] *= 1.01;
      }
    }
  }
  for(auto f : BDSAcceleratorModel::Instance()->fields) delete f;
  BDSAcceleratorModel::Instance()->fields.resize(0);
//  G4FieldManager* fieldManager = G4TransportationManager::GetTransportationManager()->GetFieldManager();
//  fieldManager->GetChordFinder()->ResetStepEstimate();
//  fieldManager->ClearChordFinders();
  realWorld->ConstructSDandField();
//  G4RunManager::GetRunManager()->GeometryHasBeenModified();

  if (initialisationResult > 1 || !initialised)
    {return;} // a mode where we don't do anything

  G4cout.precision(10);
  /// Catch aborts to close output stream/file. perhaps not all are needed.
  struct sigaction act;
  act.sa_handler = &BDS::HandleAborts;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;
  if (!ignoreSIGINT)
    {sigaction(SIGINT,  &act, nullptr);}
  sigaction(SIGABRT, &act, nullptr);
  sigaction(SIGTERM, &act, nullptr);
  sigaction(SIGSEGV, &act, nullptr);
  
  /// Run in either interactive or batch mode
  try
    {
      if (!BDSGlobalConstants::Instance()->Batch())   // Interactive mode
        {
          BDSVisManager visManager = BDSVisManager(BDSGlobalConstants::Instance()->VisMacroFileName(),
                                                   BDSGlobalConstants::Instance()->Geant4MacroFileName(),
                                                   realWorld,
                                                   BDSGlobalConstants::Instance()->VisVerbosity());
          visManager.StartSession(argcCache, argvCache);
        }
      else
        {// batch mode
          if (nGenerate < 0)
            {runManager->BeamOn(BDSGlobalConstants::Instance()->NGenerate());}
          else
            {runManager->BeamOn(nGenerate);}
        }
    }
  catch (const BDSException& exception)
    {
      // don't do this for now in case it's dangerous and we try tracking with open geometry
      //G4GeometryManager::GetInstance()->OpenGeometry();
      throw exception;
    }

}

std::vector<std::array<double, 4> > CNBDSIM::CalcBeamPars(){
    std::vector<std::array<double, 4> > ssemPred;
    std::cout<<"num events in ssem1 "<< bdsOutput->samplerTrees[0]->x.size() <<std::endl;
    ssemPred.reserve(bdsOutput->samplerTrees.size());
    for(int i=0; i<bdsOutput->samplerTrees.size(); i++){ //loop over ssems
      std::vector<float> x = (*bdsOutput->samplerTrees[i]).xcache;
      double meanx = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
      // Compute Standard Deviation using lambda
      double stddevx = std::sqrt(std::accumulate(x.begin(), x.end(), 0.0,
        [meanx](float acc, float val) { return acc + (val - meanx) * (val - meanx); }) / (x.size()-1));


      std::vector<float> y = bdsOutput->samplerTrees[i]->ycache;
      double meany = std::accumulate(y.begin(), y.end(), 0.0) / y.size();

      // Compute Standard Deviation using lambda
      double stddevy = std::sqrt(std::accumulate(y.begin(), y.end(), 0.0,
        [meany](float acc, float val) { return acc + (val - meany) * (val - meany); }) / (y.size()-1));
      
      ssemPred[i][0] = meanx;
      ssemPred[i][1] = meany;
      ssemPred[i][2] = stddevx;
      ssemPred[i][3] = stddevy;
    }
  return ssemPred;
}


CNBDSIM::~CNBDSIM()
{
  /// Termination & clean up.
  G4GeometryManager::GetInstance()->OpenGeometry();
  
#ifdef BDSDEBUG
  G4cout << __METHOD_NAME__ << "deleting..." << G4endl;
#endif
  delete bdsOutput;

  try
    {
      // order important here because of singletons relying on each other
      delete BDSSDManager::Instance();
      delete BDSBeamPipeFactory::Instance();
      delete BDSCavityFactory::Instance();
      delete BDSGeometryFactory::Instance();
      delete BDSAcceleratorModel::Instance();
      delete BDSTemporaryFiles::Instance();
      delete BDSFieldFactory::Instance(); // this uses BDSGlobalConstants which uses BDSMaterials
      delete BDSGlobalConstants::Instance();
      delete BDSMaterials::Instance();
      
      // instances not used in this file, but no other good location for deletion
      if (initialisationResult < 2)
        {
          delete BDSColours::Instance();
          delete BDSFieldLoader::Instance();
          delete BDSSamplerRegistry::Instance();
          BDSAperturePointsCache::Instance()->ClearCachedFiles();
        }
    }
  catch (...)
    {;} // ignore any exception as this is a destructor
  
  delete runManager;
  delete bdsBunch;
  delete parser;

  if (usualPrintOut)
    {G4cout << __METHOD_NAME__ << "End of Run. Thank you for using BDSIM!" << G4endl;}
}


void CNBDSIM::RegisterUserComponent(const G4String& componentTypeName,
                                  BDSComponentConstructor* componentConstructor)
{
  if (initialised)
    {BDS::Warning(__METHOD_NAME__, "BDSIM kernel already initialised - this component will not be available");}
  if (!userComponentFactory)
    {userComponentFactory = new BDSComponentFactoryUser();}

  userComponentFactory->RegisterComponent(componentTypeName,
                                          componentConstructor);
}

void CNBDSIM::GeneratePrimariesOnly(const BDSGlobalConstants* globals)
{
  // output creation is duplicated below but with this if loop, we exit so ok.
  bdsOutput->NewFile();
  const G4int nToGenerate = globals->NGenerate();
  const G4int printModulo = globals->PrintModuloEvents();
  bdsBunch->BeginOfRunAction(nToGenerate, globals->Batch());
  auto flagsCache(G4cout.flags());
  for (G4int i = 0; i < nToGenerate; i++)
    {
      if (i%printModulo == 0)
        {G4cout << "\r Primary> " << std::fixed << i << " of " << nToGenerate << G4endl;}
      BDSParticleCoordsFullGlobal coords = bdsBunch->GetNextParticleValid();
      // always pull particle definition in case it's updated
      const BDSParticleDefinition* pDef = bdsBunch->ParticleDefinition();
      bdsOutput->FillEventPrimaryOnly(coords, pDef);
    }
  G4cout.flags(flagsCache); // restore cout flags
  // Write options now the file is open
  const GMAD::OptionsBase* ob = BDSParser::Instance()->GetOptionsBase();
  bdsOutput->FillOptions(ob);
  
  // Write beam
  const GMAD::BeamBase* bb = BDSParser::Instance()->GetBeamBase();
  bdsOutput->FillBeam(bb);
  
  bdsOutput->CloseFile();
}
