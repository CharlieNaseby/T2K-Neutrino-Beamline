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
#ifndef CNBDSIMCLASS_H
#define CNBDSIMCLASS_H

class BDSBunch;
class BDSEventAction;
class BDSComponentConstructor;
class BDSComponentFactoryUser;
class BDSDetectorConstruction;
class BDSGlobalConstants;
class BDSOutput;
class BDSParser;
class BDSRunManager;
class G4VModularPhysicsList;
class BDSPrimaryGeneratorAction;
class BDSParticleDefinition;
class BDSGlobalConstants;
class BDSAcceleratorComponentRegistry;

#include "G4String.hh"
#include <vector>
#include <array>
#include <numeric>
#include <map>


/** 
 * @brief Interface class to use BDSIM.
 *
 * First way to use:
 * bds = new BDSIM(argc, argv);
 * bds->BeamOn();
 * 
 * Second way to use (delayed initialisation):
 * bds = new BDSIM();
 * // modify bds in some way
 * bds->Initialise(argc, argv);
 * bds->BeamOn();
 *
 * @author Laurie Nevay
 */

class CNBDSIM
{
public:
  /// Construct an instance but don't initialise. Requires initialisation with
  /// arguments argc and arv
  CNBDSIM();

  /// Initialise everything given these arguments.
  int Initialise(int argc, char** argv, bool usualPrintOut=true);

  double GetParameterValue(std::string key);

  /// Construct and initialise BDSIM.
  CNBDSIM(int argc, char** argv, bool usualPrintOut=true);

  /// The destructor opens the geometry in Geant4 and deletes everything.
  ~CNBDSIM();

  /// Accessor as to whether BDSIM kernel is initialised - ie all geometry and physics
  /// constructed.
  inline bool Initialised() const {return initialised;};

  inline int InitialisationResult() const {return initialisationResult;}

  /// Generate nGenerate events. If the default argument -1 is used, the number is taken
  /// from the standard input e.g. the executable option ngenerate and then the one specified
  /// in the input gmad files as an option.
  void BeamOn(int nGenerate=-1, std::map<std::string, double> pars = {{std::string("NULL"), -999.0}} );
  std::vector<std::array<double, 4> > CalcBeamProperties();
  /// Register a custom user beam line element by the type name you'd like it to have
  /// and the (user-provided) constructor that can construct it.
  void RegisterUserComponent(const G4String& componentTypeName,
			     BDSComponentConstructor* componentConstructor);

  /// Generate the primaries, fill in the output and close the output. In a function to
  /// simplify the main.
  void GeneratePrimariesOnly(const BDSGlobalConstants* globals);
  
  /// Provide a physics list that will be used inplace of the BDSIM generate one.
  void RegisterUserPhysicsList(G4VModularPhysicsList* userPhysicsListIn) {userPhysicsList = userPhysicsListIn;}
  G4VModularPhysicsList* UserPhysicsList() const {return userPhysicsList;} ///< Access user physics list.
  
//private:
  /// The main function where everything is constructed.
  int Initialise();
  void SetFileWriting(bool write);
  G4String randomState; 
  bool   ignoreSIGINT;         ///< For cmake testing.
  bool   usualPrintOut;        ///< Whether to allow the usual cout output.
  bool   initialised;          ///< Whether initialisation was completed safely
  int    initialisationResult; ///< Possible to not finish initialisation but have completed ok - flag for this.
  int    argcCache;            ///< Cache of argc.
  char** argvCache;            ///< Cache of argv.

  /// @{ Cache of main objects in BDSIM.
  BDSParser*     parser;
  BDSOutput*   bdsOutput;
  BDSBunch*      bdsBunch;
  BDSEventAction* eventAction;
  BDSRunManager* runManager;
  G4VModularPhysicsList* physList;
  BDSComponentFactoryUser* userComponentFactory; ///< Optional user registered component factory.
  G4VModularPhysicsList* userPhysicsList;        ///< Optional user registered physics list.
  BDSDetectorConstruction* realWorld;
  BDSPrimaryGeneratorAction* primaryGeneratorAction;
  BDSParticleDefinition* beamParticle = nullptr;
  BDSGlobalConstants* globals;
  /// @}
};

#endif
