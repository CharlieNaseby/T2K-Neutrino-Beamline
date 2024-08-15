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


  int nPars = 11;

  double magset[31];
  TFile *inf = new TFile(ssemDataFile.c_str(), "READ");
  TTree *anabeam = dynamic_cast<TTree*>(inf->Get("anabeam"));
  anabeam->SetBranchAddress("magset", magset);
  anabeam->GetEntry(0);

  std::string magnames[11];
  magnames[0] = "BPV1";
  magnames[1] = "BPH2";
  magnames[2] = "QPQ1";
  magnames[3] = "QPQ2";
  magnames[4] = "BPD1";
  magnames[5] = "BPD2";
  magnames[6] = "QPQ3";
  magnames[7] = "BPV2";
  magnames[8] = "QPQ4";
  magnames[9] = "BPH3";
  magnames[10] = "QPQ5";

  double fudgeFactor[11];
  for(int i=0; i<nPars; i++) fudgeFactor[i] = 1.0;

  fudgeFactor[1] = 1.402333;

  fudgeFactor[4] = 0.851299;
  fudgeFactor[5] = 0.830461;
  fudgeFactor[7] = 1.023666;


  double pars[11];
  pars[0] = magset[0]/100.;
  pars[1] = magset[1]/100.;
  pars[2] = magset[2]/100.;
  pars[3] = magset[4]/100.;
  pars[4] = magset[5]/100.;
  pars[5] = magset[6]/100.;
  pars[6] = magset[7]/100.;
  pars[7] = magset[8]/100.;
  pars[8] = magset[9]/100.;
  pars[9] = magset[10]/100.;
  pars[10] = magset[11]/100.;

  for(int i=0; i<nPars; i++) pars[i] *= fudgeFactor[i];

  inter.SetNPars(nPars);

  pars[0] = -0.05;  //CERN TODO forcefully set BPV1

  inter.SetInternalPars(pars);
  inter.SetNominalPars(pars);

  for(int i=0; i<nPars; i++) std::cout<<magnames[i]<<"\t"<<pars[i]<<std::endl;



  auto iterstarttime = std::chrono::high_resolution_clock::now();
  inter.fcn(pars);
  auto iterendtime = std::chrono::high_resolution_clock::now();
  auto itertime = std::chrono::duration_cast<std::chrono::microseconds>(iterendtime-iterstarttime).count();
  std::cout<<"Took "<<itertime*1e-6<<"s to run a single iteration"<<std::endl;


  exit(1);

  ROOT::Math::Functor f(inter, &Interface::fcn_wrapper, nPars);
  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetStrategy(3);
  min->SetFunction(f);
  min->SetMaxFunctionCalls(10000);
  min->SetTolerance(10);

  for(int i=0; i<nPars; i++) min->SetVariable(i, magnames[i].c_str(), pars[i], 0.1);

  std::vector<int> bMagVars = {0, 1, 2, 3, 6, 7, 8, 9, 10};
  std::vector<int> qMagVars = {2, 3, 6, 8, 10};

//  for(auto fixed : qMagVars) min->FixVariable(fixed);

  for(int i=0; i<nPars; i++) min->FixVariable(i);
  min->SetPrintLevel(2);

  min->ReleaseVariable(0); 
  min->ReleaseVariable(1);
  min->ReleaseVariable(4);
  min->ReleaseVariable(7);
  inter.SetChisqMode(1);
  min->Minimize();
  min->PrintResults();

//  for(int i=0; i<nPars; i++) min->ReleaseVariable(i);
//  for(auto fixed : bMagVars) min->FixVariable(fixed);
//
//  inter.SetChisqMode(2);
//  min->Minimize();
//  min->PrintResults();


 
//  TFile *result = new TFile("result.root", "RECREATE");
//  TH1D *scan = new TH1D("param_4", "param_4", 20, 0., 15);

//  inter.ParamScan(4, scan);
//  result->cd();
//  scan->Write();
//  result->Close();

  auto endtime = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::microseconds>(endtime-starttime).count();
  std::cout<<"Took "<<time*1e-6<<"s to run"<<std::endl;

  return 0;

}
