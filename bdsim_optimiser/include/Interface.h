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
#include "Parameter.h"


#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TRandom3.h"


#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#endif


#define NSSEM 9

class Interface{
private:
  std::shared_ptr<std::vector<std::shared_ptr<Parameter>>> parametersVec;
  std::shared_ptr<std::unordered_map<std::string, std::shared_ptr<Parameter>>> parametersMap;
public:
  std::vector<std::array<double, 4> > dat;
  std::vector<double> s;
  std::vector<double> fcnHistory, MAEHistory;
//  std::vector<std::string> beamline;
  int nPars, nMagnetPars, nBeamPars;
  std::vector<double> internalPars;
  std::map<std::string, double> nominalPars;

  std::vector<Parameter> parameters;

  std::vector<double> magCurrent;
  std::vector<double> preFit;
  std::map<std::string, int> magMap;
  std::vector<std::string> magNames;
  std::vector<std::string> parNames;
  CNBDSIM *bds;
  std::vector<int> ssemMask;

  int testval = 0;
  unsigned int fitMode=1+2+4+8;  //by default fit width and position
 
  Interface(std::string dataFile, std::string baseBeamlineFile,
    std::shared_ptr<std::unordered_map<std::string, std::shared_ptr<Parameter>>> parametersMapIn,
    std::shared_ptr<std::vector<std::shared_ptr<Parameter>>> parametersVecIn);
   
  ~Interface();

  void SetInitialValues(char *usePrevBestFit, bool useFieldMaps, bool useInputFile, double noise);
  void SetFileWriting(bool write);
  void SetSSEMData(std::string dataFile);
  void SetData(std::vector<std::array<double, 4> > target);
  void SetData(std::vector<double> x, std::vector<double> y, std::vector<double> wx, std::vector<double> wy);
  void SetInternalPars(const double *pars);
  std::map<std::string, double> GetNominalPars();
  void ParamScan(std::string param, TH1D *hist);
  void ParamScan(int param, int npoints, double xlo, double xup);
  bool CheckBounds(std::map<std::string, double> pars);
  double fcn();
  double fcn(std::vector<double> pars);
  double fcn_wrapper(const double *pars);
  double fcn(const double *pars);
//  double MAE(const double *pars);
//  double MAE(std::vector<double> pars);
//  void GenerateInputFile(const double *pars);
  void ParseInputFile(std::string baseBeamlineFile);
  std::vector<std::array<double, 4> > GetBeamProperties();
//  std::vector<double> FitToPhysical(std::vector<double> fitval);
  double FitToPhysical(std::string name, double fitval);
//  std::vector<double> PhysicalToFit(std::vector<double> physval);
  double PhysicalToFit(std::string name, double physval);
  double CalcChisq(const double *pars);
//  double CalcMAE(const double *pars);
  std::map<std::string, double> GetParMap(const double *pars);
  void BeamOn(int n, std::map<std::string, double> parmap);
  double CalcPrior(std::map<std::string, double> pars);
  void TestBdsim();
  void SetChisqMode(int mode){fitMode=mode;};
  void SetSSEMMask(std::vector<int> mask){ssemMask = mask;};
  std::vector<double> GetFcnHistory(){return fcnHistory;};
  std::vector<double> GetMAEHistory(){return MAEHistory;};

};

#ifdef PYBIND
PYBIND11_MODULE(Interface, m) {
    pybind11::class_<Interface>(m, "Interface")
        .def(pybind11::init<std::string &, std::string &, int, int, int>())
        .def(pybind11::init<std::vector<double> &, std::vector<std::array<double, 4> > &>())
        .def("fcn", static_cast<double (Interface::*)(std::vector<double>)> (&Interface::fcn))
        .def("MAE", static_cast<double (Interface::*)(std::vector<double>)> (&Interface::MAE))
        .def("SetInitialValues", static_cast<void (Interface::*)(std::vector<double> )> (&Interface::SetInitialValues))
        .def("GetBeamProperties", static_cast<std::vector<std::array<double, 4> > (Interface::*)()>(&Interface::GetBeamProperties))
        .def("GetFcnHistory", static_cast<std::vector<double> (Interface::*)()>(&Interface::GetFcnHistory))
        .def("SetData", static_cast<void (Interface::*)(std::vector<std::array<double, 4> > ) >(&Interface::SetData))
        .def("SetFileWriting", static_cast<void (Interface::*)(bool ) >(&Interface::SetFileWriting))
        .def("SetChisqMode", static_cast<void (Interface::*)(int)> (&Interface::SetChisqMode));
}
#endif