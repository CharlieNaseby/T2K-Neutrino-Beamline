#include <iostream>
#include <set>
#include <string>
#include <vector>

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
  std::vector<double> dat;
  Interface(double *data);
  ~Interface();
  double calc_chisq(const double *pars);

};
