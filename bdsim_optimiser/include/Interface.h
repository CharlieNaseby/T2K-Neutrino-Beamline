#include <iostream>
#include <fstream>
#include <iomanip>
#include <set>
#include <string>
#include <vector>
#include <array>

#include "Beam.hh"
#include "Config.hh"
#include "DataLoader.hh"
#include "EventAnalysis.hh"
#include "Options.hh"
#include "RBDSException.hh"

#include "BDSOutputROOTEventBeam.hh"
#include "BDSOutputROOTEventHeader.hh"
#include "BDSOutputROOTEventOptions.hh"
#include "BDSOutputROOTEventSampler.hh"
#include "AnalysisUtilities.hh"


#include "TFile.h"
#include "TVector.h"
#include "TChain.h"
#include "TTree.h"

#include "BDSIMClass.hh"
#include "BDSException.hh"
#include "CNBDSIMClass.h"


#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"


#define NSSEM 9

class Interface{
public:
  std::vector<std::array<double, 4> > dat;
  std::vector<double> s;
  std::vector<std::string> beamline;
  int nPars, nMagnetPars, nBeamPars;
  std::vector<double> internalPars;
  std::map<std::string, double> nominalPars;
  std::vector<double> magCurrent;
  std::vector<double> preFit;
  std::map<std::string, int> magMap;
  std::vector<std::string> magNames;
  std::vector<std::string> beamNames;
  std::vector<std::string> parNames;
  std::map<std::string, double> priorErrors;
  CNBDSIM *bds;

  int testval = 0;
  unsigned int fitMode=1+2+4+8;  //by default fit width and position
 
  Interface(std::string dataFile, std::string baseBeamlineFile, int npars, int nmagnetpars, int nbeampars);
  ~Interface();
  void SetInitialValues(bool usePrevBestFit, bool useFieldMaps, bool useFudgeFactor, bool useInputFile, double* pars);
  void SetInternalPars(const double *pars);
  std::map<std::string, double> GetNominalPars();
  void ParamScan(int param, TH1D *hist);
  void ParamScan(int param, int npoints, double xlo, double xup);
  bool CheckBounds(std::map<std::string, double> pars);
  double fcn(std::vector<double> pars);
  double fcn_wrapper(const double *pars);
  double fcn(const double *pars);
  void GenerateInputFile(const double *pars);
  void ParseInputFile(std::string baseBeamlineFile);
  std::vector<std::array<double, 4> > GetBeamPars();
  double CalcChisq(const double *pars);
  double CalcPrior(std::map<std::string, double> pars);
  void TestBdsim();
  void SetChisqMode(int mode){fitMode=mode;};
  double PrintTestVal(const double *pars){
    std::cout<<"testval = "<<testval<<std::endl;
    return 1.0;
  }

};
