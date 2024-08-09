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
#include "AnalysisUtilities.hh"


#include "TFile.h"
#include "TChain.h"

#include "BDSIMClass.hh"
#include "BDSException.hh"

#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"

class Interface{
public:
  std::vector<std::array<double, 4> > dat;
  std::vector<double> s;
  std::vector<std::string> beamline;
  int nPars;
  std::vector<double> internalPars;
  Interface(std::vector<std::array<double, 4> > data, std::string baseBeamlineFile);
  ~Interface();
  void SetNPars(int npars);
  void SetInternalPars(const double *pars);
  void ParamScan(int param, TH1D *hist);
  void ParamScan(int param, int npoints, double xlo, double xup);
  double fcn(std::vector<double> pars);
  double fcn_wrapper(const double *pars);
  double fcn(const double *pars);
  void GenerateInputFile(const double *pars);
  void ParseInputFile(std::string baseBeamlineFile);
  double CalcChisq();

};
