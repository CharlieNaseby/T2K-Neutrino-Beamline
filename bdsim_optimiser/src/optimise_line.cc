#include "Interface.h"
#include "TClassTable.h"
#include <chrono>
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TMatrixD.h"
#include "TVectorD.h"

double performFit(ROOT::Math::Minimizer *min, Interface *inter, int nPars, double *pars, std::initializer_list<std::vector<int> > fixedVars, int fitMode){
  inter->SetChisqMode(fitMode);
  min->SetStrategy(3);
  min->SetMaxFunctionCalls(10000);
  min->SetTolerance(1000);

  for(int i=0; i<nPars; i++) min->SetVariable(i, inter->parNames[i], pars[i], 0.1);
 

  for(auto fixedvec : fixedVars)
    for(auto fixed : fixedvec) 
      min->FixVariable(fixed);
  
  min->SetPrintLevel(2);

  min->FixVariable(0); //BPV1 is set to 0 current
  min->FixVariable(9); //BPH3 and QPQ5 do not have and SSEMs afterwards in the sim so would have no data constraint
  min->FixVariable(10);
  min->Minimize();
  min->PrintResults();
  for(int i=0; i<nPars; i++) pars[i] = min->X()[i];
  return min->MinValue();
}

int main(int argc, char **argv){
  auto starttime = std::chrono::high_resolution_clock::now();

  std::string baseBeamlineFile="../gmad/optimised_test.gmad";
  std::string ssemDataFile="./ssem_data/run0910216_gen.root";

  const int nMagnetPars = 11;
  const int nBeamPars = 10;
  const int nPars = nMagnetPars + nBeamPars;

  Interface inter(ssemDataFile, baseBeamlineFile, nPars, nMagnetPars, nBeamPars);

  double pars[nPars];

  bool usePrevBestFit = false;
  bool useFieldMaps = false; //currently unsupported
  bool useFudgeFactor = false;
  bool useInputFile = true;


  inter.SetInitialValues(usePrevBestFit, useFieldMaps, useFudgeFactor, useInputFile, pars);

  //set prior constraints
  //negative values are a fractional uncertainty on the nominal value
  //positive values are an absolute uncertainty in the same units as the parameter

  inter.priorErrors["BPV1"] = 0.00001; //nominal current is 0 maybe there should actually be no freedom here..
  inter.priorErrors["BPH2"] = -0.2;
  inter.priorErrors["QPQ1"] = -0.2;
  inter.priorErrors["QPQ2"] = -0.2;
  inter.priorErrors["BPD1"] = -0.2;
  inter.priorErrors["BPD2"] = -0.2;
  inter.priorErrors["QPQ3"] = -0.2;
  inter.priorErrors["BPV2"] = -0.2;
  inter.priorErrors["QPQ4"] = -0.2;
  inter.priorErrors["BPH3"] = -0.2;
  inter.priorErrors["QPQ5"] = -0.2;

  //really lose constraints on beam parameters, but does really help with them not exploding 
  inter.priorErrors["X0"] = 2;
  inter.priorErrors["emitx"] = -0.2;
  inter.priorErrors["betx"] = 10;
  inter.priorErrors["alfx"] = 5;

  inter.priorErrors["Y0"] = 2;
  inter.priorErrors["emity"] = -0.2;
  inter.priorErrors["bety"] = 10;
  inter.priorErrors["alfy"] = 5;

  //pars is now set to the expected prefit values of the magnets based on the bools above
  inter.SetInternalPars(pars);
  //inter.SetNominalPars(pars);

  for(int i=0; i<nPars; i++) std::cout<<inter.parNames[i]<<"\t"<<pars[i]<<std::endl;

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
  std::vector<int> beamParVars = {11, 12, 13, 14, 15, 16, 17, 18, 19, 20};

  auto wrappedFcn = [&inter](const double* pars) {
      return inter.fcn_wrapper(pars);
  };
 
  ROOT::Math::Functor f(wrappedFcn, nPars);
  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

  min->SetFunction(f);
  min->SetVariableLowerLimit(13, 0); //emitx
  min->SetVariableLowerLimit(14, 0); //betax
  min->SetVariableLowerLimit(18, 0); //emity
  min->SetVariableLowerLimit(19, 0); //betay

  for(int i=0; i<inter.parNames.size(); i++) std::cout<< i<<"  " << inter.parNames[i] <<std::endl;
//  TFile *result = new TFile("result.root", "RECREATE");
//  TH1D *hist;
//  for(int i=0; i<nPars; i++){
//    std::string name = "param_";
//    name += inter.parNames[i];
//    if(pars[i]>=0) hist = new TH1D(name.c_str(), name.c_str(), 11, pars[i]*0.5, pars[i]*1.5);
//    else hist = new TH1D(name.c_str(), name.c_str(), 11, pars[i]*1.5, pars[i]*0.5);
//
//    inter.ParamScan(i, hist);
//
//    result->cd();
//    hist->Write(name.c_str());
//  }
//  result->Close();



//  performFit(min, &inter, nPars, pars, {qMagVars}, 1+2);
//  performFit(min, &inter, nPars, pars, {bMagVars}, 4+8);
//
//  performFit(min, &inter, nPars, pars, {qMagVars}, 1+2);
//  performFit(min, &inter, nPars, pars, {bMagVars}, 4+8);

//  performFit(min, &inter, nPars, pars, {}, 1+2+4+8);

//  std::cout << "about to call perform fit with B magnets only"<<std::endl;
//  performFit(min, &inter, nPars, pars, {qMagVars, beamParVars}, 1+2);
//  std::cout << "about to call perform fit with Q magnets only"<<std::endl;
//  performFit(min, &inter, nPars, pars, {bMagVars, beamParVars}, 4+8);
//  std::cout << "about to call perform fit with beam only"<<std::endl;
//  performFit(min, &inter, nPars, pars, {bMagVars, qMagVars}, 1+2+4+8);
//
//  std::cout << "about to call perform fit with B magnets only take 2"<<std::endl;
//  performFit(min, &inter, nPars, pars, {qMagVars, beamParVars}, 1+2);
//  std::cout << "about to call perform fit with Q magnets only take 2"<<std::endl;
//  performFit(min, &inter, nPars, pars, {bMagVars, beamParVars}, 4+8);
//  std::cout << "about to call perform fit with beam only take 2"<<std::endl;
//  performFit(min, &inter, nPars, pars, {bMagVars, qMagVars}, 1+2+4+8);
  std::cout << "about to call perform fit with all parameters free (this may take a while)"<<std::endl;
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
    nom[i] = inter.nominalPars[inter.parNames[i]];
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
 
  
 
  auto endtime = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::microseconds>(endtime-starttime).count();
  std::cout<<"Took "<<time*1e-6<<"s to run"<<std::endl;
 
  return 0;

}
