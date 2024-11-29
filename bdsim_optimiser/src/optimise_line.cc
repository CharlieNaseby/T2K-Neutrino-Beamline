#include "Interface.h"
#include "TClassTable.h"
#include <chrono>
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TMatrixD.h"
#include "TVectorD.h"

double performFit(ROOT::Math::Minimizer *min, Interface *inter, int nPars, double *pars, std::vector<int> fixedVars, int fitMode){
  inter->SetChisqMode(fitMode);
  min->SetStrategy(3);
  min->SetMaxFunctionCalls(10000);
  min->SetTolerance(10);

  for(int i=0; i<nPars; i++) min->SetVariable(i, inter->magNames[i], pars[i], 0.1);

  for(auto fixed : fixedVars) min->FixVariable(fixed);

  min->SetPrintLevel(2);

  min->FixVariable(9); //BPH3 and QPQ5 do not have and SSEMs afterwards in the sim so would have no data constraint
  min->FixVariable(10);
  min->Minimize();
  min->PrintResults();
  for(int i=0; i<nPars; i++) pars[i] = min->X()[i];
  return min->MinValue();
}

int main(int argc, char **argv){
  auto starttime = std::chrono::high_resolution_clock::now();

  std::string baseBeamlineFile="../gmad/test.gmad";
  std::string ssemDataFile="./ssem_data/run0910216_gen.root";

  const int nPars = 11;
  Interface inter(ssemDataFile, baseBeamlineFile, nPars);

  double pars[nPars];

  bool usePrevBestFit = false;
  bool useFieldMaps = false;
  bool useFudgeFactor = false;


  inter.SetInitialValues(usePrevBestFit, useFieldMaps, useFudgeFactor, pars);

  //pars is now set to the expected prefit values of the magnets based on the bools above
  inter.SetInternalPars(pars);
  //inter.SetNominalPars(pars);

  for(int i=0; i<nPars; i++) std::cout<<inter.magNames[i]<<"\t"<<pars[i]<<std::endl;

  inter.SetChisqMode(1+2+4+8);
 
//just see how long a single iteration takes
 
  auto iterstarttime = std::chrono::high_resolution_clock::now();
  inter.fcn(pars);
  auto iterendtime = std::chrono::high_resolution_clock::now();
  auto itertime = std::chrono::duration_cast<std::chrono::microseconds>(iterendtime-iterstarttime).count();
  std::cout<<"Took "<<itertime*1e-6<<"s to run a single iteration"<<std::endl;
 
 
 
//now for fitting
 
  std::vector<int> bMagVars = {0, 1, 4, 5, 7, 9};
  std::vector<int> qMagVars = {2, 3, 6, 8, 10};

  auto wrappedFcn = [&inter](const double* pars) {
      return inter.fcn_wrapper(pars);
  };
 
  ROOT::Math::Functor f(wrappedFcn, nPars);
  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

  min->SetFunction(f);
//  performFit(min, &inter, nPars, pars, qMagVars, 1+2);
//  performFit(min, &inter, nPars, pars, bMagVars, 4+8);
//
//  performFit(min, &inter, nPars, pars, qMagVars, 1+2);
//  performFit(min, &inter, nPars, pars, bMagVars, 4+8);

  performFit(min, &inter, nPars, pars, {}, 1+2+4+8);

//  for(auto fixed : qMagVars) min->FixVariable(fixed);
//  for(auto fixed : bMagVars) min->FixVariable(fixed);
 
//  for(int i=0; i<nPars; i++) min->FixVariable(i);
//  min->SetPrintLevel(2);
 
//  min->ReleaseVariable(0); 
  // min->FixVariable(9); //BPH3 and QPQ5 do not have and SSEMs afterwards in the sim so would have no data constraint
  // min->FixVariable(10);
  // min->Minimize();
  // min->PrintResults();

 
 
  TFile *outf = new TFile("fit_results.root", "RECREATE");
 
  double *covarray = new double[nPars*nPars];
  min->GetCovMatrix(covarray);
  TMatrixD *cov = new TMatrixD(nPars, nPars, covarray);
  cov->Write("postfit_covariance");
  const double *errors = min->Errors();
  const double *bestFit = min->X();
 
  TVectorD nom(nPars);
  TVectorD pre(nPars);
  TVectorD pos(nPars);
  TVectorD posError(nPars);
  
  for(int i=0; i<nPars; i++){
    nom[i] = inter.nominalPars[i];
    pre[i] = inter.preFit[i];
    pos[i] = bestFit[i];
    posError[i] = errors[i];
  }
 
  nom.Write("nominal");
  pre.Write("preFit");
  pos.Write("postFit");
  posError.Write("postFitError");
  outf->Close();
 
//  for(int i=0; i<nPars; i++) min->ReleaseVariable(i);
//  for(auto fixed : bMagVars) min->FixVariable(fixed);
 
  
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
