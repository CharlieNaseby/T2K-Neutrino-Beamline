#include "Interface.h"
#include "TClassTable.h"
#include <chrono>
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"


int main(int argc, char **argv){
  auto starttime = std::chrono::high_resolution_clock::now();


  std::string baseBeamlineFile="../gmad/test.gmad";
  std::string ssemDataFile="./ssem_data/run0910580_gen.root";
  Interface inter(ssemDataFile, baseBeamlineFile);
  const double pars[3] = {1, 1, 1};
  int nPars = 3;
  inter.SetNPars(nPars);
  //scan a parameter

  inter.SetInternalPars(pars);

  ROOT::Math::Functor f(inter, &Interface::fcn_wrapper, nPars);
  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetStrategy(3);
  min->SetFunction(f);
  min->SetMaxFunctionCalls(10000);

  min->SetVariable(0, "BPD1", 0.839889, 0.1);
  min->SetVariable(1, "BPD2", 0.830461, 0.1);
  min->SetVariable(2, "QPQ4", 1.0, 0.1);
//  min->FixVariable(0);

  min->Minimize();
  min->PrintResults();
 
//  TFile *result = new TFile("result.root", "RECREATE");
//  TH1D *scan = new TH1D("param_0", "param_0", 20, 0.5, 1.0);
//
//  inter.ParamScan(0, scan);
//  result->cd();
//  scan->Write();
//  result->Close();

  auto endtime = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::microseconds>(endtime-starttime).count();
  std::cout<<"Took "<<time*1e-6<<"s to run"<<std::endl;

  return 0;

}
