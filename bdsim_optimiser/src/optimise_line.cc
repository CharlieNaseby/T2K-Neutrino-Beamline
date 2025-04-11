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

//  for(int i=0; i<nPars; i++) min->SetVariable(i, inter->parNames[i], pars[i], 0.1);
   for(int i=0; i<nPars; i++) min->SetVariable(i, inter->parNames[i], inter->PhysicalToFit(i, pars[i]), 0.1);


  for(auto fixedvec : fixedVars)
    for(auto fixed : fixedvec) 
      min->FixVariable(fixed);
  
  min->SetPrintLevel(2);

  min->FixVariable(0); //BPV1 is set to 0 current
  min->FixVariable(9); //BPH3 and QPQ5 do not have and SSEMs afterwards in the sim so would have no data constraint
  min->FixVariable(10);
  min->Minimize();
  min->PrintResults();
  for(int i=0; i<nPars; i++) pars[i] = inter->FitToPhysical(i, min->X()[i]);
  return min->MinValue();
}

void saveResult(ROOT::Math::Minimizer *min, Interface *inter, int nPars, const char *filename){
  double fcnmin = inter->fcn_wrapper(min->X());
  std::cout<<"saving minimum at fcn "<<fcnmin<<std::endl;
  TFile *outf = new TFile(filename, "RECREATE");
  TTree *paramtree = new TTree("parameters", "parameters");
  double parameter_fit, parameter_physical, error_fit, error_physical;
  int param_num;
  TString name;
  paramtree->Branch("fit_value", &parameter_fit, "fit_value/D");
  paramtree->Branch("physical_value", &parameter_physical, "physical_value/D");
  paramtree->Branch("name", &name);
  paramtree->Branch("param_num", &param_num, "param_num/I");
  paramtree->Branch("fit_error", &error_fit, "fit_error/D");
  paramtree->Branch("physical_error", &error_physical, "physical_error/D");

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
  TVectorD posFitBasis(nPars);
  TVectorD posErrorFitBasis(nPars);
  TVectorD posFcn(2);
  posFcn[0] = fcnmin;
  posFcn[1] = inter->CalcPrior(inter->GetParmap(min->X()));

  for(int i=0; i<nPars; i++){
    nom[i] = inter->nominalPars[inter->parNames[i]];
    pre[i] = inter->preFit[i];
    pos[i] = inter->FitToPhysical(i, bestFit[i]);
    posError[i] = inter->FitToPhysical(i, errors[i]);
    posFitBasis[i] = bestFit[i];
    posErrorFitBasis[i] = errors[i];

    param_num = i;
    parameter_fit = bestFit[i];
    parameter_physical = inter->FitToPhysical(i, bestFit[i]);
    error_fit = errors[i];
    error_physical = inter->FitToPhysical(i, errors[i]);
    name = TString(inter->parNames[i]);
    paramtree->Fill();
  }

  paramtree->Write();
  std::vector<std::array<double, 4> > allSSEMSimulation = inter->GetBeamPars();
  std::vector<std::array<double, 4> > dat = inter->dat;

  for(int i=0; i<4; i++){
    TVectorD SSEM_data(NSSEM);
    TVectorD SSEM_sim(NSSEM);
    const char *dataNames[4] = {"xmean_data", "ymean_data", "xwidth_data", "ywidth_data"};
    const char *simNames[4] = {"xmean_sim", "ymean_sim", "xwidth_sim", "ywidth_sim"};
    for(int j=0; j<NSSEM; j++){
      SSEM_data[j] = dat[j][i];
      SSEM_sim[j] = allSSEMSimulation[j][i];
    }
    SSEM_data.Write(dataNames[i]);
    SSEM_sim.Write(simNames[i]);
  }

  nom.Write("nominal");
  pre.Write("preFit");
  pos.Write("postFit");
  posError.Write("postFitError");
  posFitBasis.Write("postFitFitBasis");
  posErrorFitBasis.Write("postFitErrorFitBasis");
  posFcn.Write("chisq");
  outf->Close();
}


int main(int argc, char **argv){
  auto starttime = std::chrono::high_resolution_clock::now();

  std::string baseBeamlineFile="../survey/unoptimised.gmad";
  std::string ssemDataFile="./ssem_data/run0910216_gen.root";

  const int nMagnetPars = 11;
  const int nBeamPars = 10;
  const int nPars = nMagnetPars + nBeamPars;

  Interface inter(ssemDataFile, baseBeamlineFile, nPars, nMagnetPars, nBeamPars);
  inter.bds->SetFileWriting(false); //dont want to save a file for every simulation run, realllly slows things down

  double pars[nPars];

  //options for the initial magnet field strengths, all false means use ssem data file and estimates of magnet strengths
  char *usePrevBestFit = nullptr; //"./mar_2025_cm_fits/fit_910216_with_misalignments_no_noise_5pc_constraint.root";
  bool useFieldMaps = false; //currently unsupported
  bool useFudgeFactor = false;
  bool useInputFile = false;


  inter.SetInitialValues(usePrevBestFit, useFieldMaps, useFudgeFactor, useInputFile, pars, 0); //last arg is noise

  //set prior constraints
  //negative values are a fractional uncertainty on the nominal value
  //positive values are an absolute uncertainty in the same units as the parameter

//  inter.priorErrors["BPV1"] = 0.00001; //nominal current is 0 maybe there should actually be no freedom here..
  inter.priorErrors["BPH2"] = -0.05;
  inter.priorErrors["QPQ1"] = -0.05;
  inter.priorErrors["QPQ2"] = -0.05;
  inter.priorErrors["BPD1"] = -0.05;
  inter.priorErrors["BPD2"] = -0.05;
  inter.priorErrors["QPQ3"] = -0.05;
  inter.priorErrors["BPV2"] = -0.05;
  inter.priorErrors["QPQ4"] = -0.05;
//  inter.priorErrors["BPH3"] = -0.05;
//  inter.priorErrors["QPQ5"] = -0.05;

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
//  inter.fcn(pars);
  auto iterendtime = std::chrono::high_resolution_clock::now();
  auto itertime = std::chrono::duration_cast<std::chrono::microseconds>(iterendtime-iterstarttime).count();
  std::cout<<"Took "<<itertime*1e-6<<"s to run a single iteration"<<std::endl;

//now for fitting

  std::vector<int> bMagVars = {0, 1, 4, 5, 7, 9};
  std::vector<int> qMagVars = {2, 3, 6, 8, 10};
  std::vector<int> beamTwissVars = {13, 14, 15, 18, 19, 20};
  std::vector<int> beamPosVars = {11, 12, 16, 17};

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

  for(unsigned int i=0; i<inter.parNames.size(); i++) std::cout<< i<<"  " << inter.parNames[i] <<std::endl;
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


  std::cout << "about to call perform fit with B magnets and beam position only"<<std::endl;
  performFit(min, &inter, nPars, pars, {qMagVars, beamTwissVars}, 1+2);
  std::cout << "about to call perform fit with Q magnets and beam twiss parameters only"<<std::endl;
  performFit(min, &inter, nPars, pars, {bMagVars, beamPosVars}, 4+8);
  std::cout << "about to call perform fit with beam only"<<std::endl;
  performFit(min, &inter, nPars, pars, {bMagVars, qMagVars}, 1+2+4+8);
 
  saveResult(min, &inter, nPars, "fit_results_after_first_split_optimisation.root");

  std::cout << "about to call perform fit with B magnets only take 2"<<std::endl;
  performFit(min, &inter, nPars, pars, {qMagVars, beamTwissVars}, 1+2);
  std::cout << "about to call perform fit with Q magnets only take 2"<<std::endl;
  performFit(min, &inter, nPars, pars, {bMagVars, beamPosVars}, 4+8);
  std::cout << "about to call perform fit with beam only take 2"<<std::endl;
  performFit(min, &inter, nPars, pars, {bMagVars, qMagVars}, 1+2+4+8);
 
  saveResult(min, &inter, nPars, "fit_results_after_second_split_optimisation.root");

  std::cout << "about to call perform fit with all parameters free (this may take a while)"<<std::endl;
  performFit(min, &inter, nPars, pars, {}, 1+2+4+8);


  saveResult(min, &inter, nPars, "fit_results.root");
  if(argc > 1) saveResult(min, &inter, nPars, argv[1]);

  inter.GenerateInputFile(min->X());

  auto endtime = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::microseconds>(endtime-starttime).count();
  std::cout<<"Took "<<time*1e-6<<"s to run"<<std::endl;
 
  return 0;

}
