#include "Interface.h"
#include "Parameter.h"
#include "TClassTable.h"
#include <chrono>
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TMatrixD.h"
#include "TVectorD.h"

//it is sometimes helpful to access via index or name so keep a synchronised vector and map of parameters
auto parameterMap = std::make_shared<std::unordered_map<std::string, std::shared_ptr<Parameter>>>();
auto parameterVec = std::make_shared<std::vector<std::shared_ptr<Parameter>>>();


double performFit(ROOT::Math::Minimizer *min, Interface *inter, std::initializer_list<std::vector<std::string>> freeGroups, int fitMode){
  inter->SetChisqMode(fitMode);
  min->SetStrategy(3);
  min->SetMaxFunctionCalls(10000);
  min->SetTolerance(1000);

  for(auto &par : *parameterVec) par->setState("fixed");
  for(auto &group : freeGroups){
    for(auto &name : group){
      parameterMap->at(name)->setState("free");
    }
  }

  for(size_t i=0; i<parameterVec->size(); i++){
    auto &par = (*parameterVec)[i];
    min->SetVariable(i, par->name, inter->PhysicalToFit(par->name, par->getPhysicalValue()), 0.1);
    if(par->getState() == "fixed") min->FixVariable(i);
  }

  min->SetPrintLevel(2);

  min->Minimize();
  min->PrintResults();
  inter->SetInternalPars(min->X()); //the result is now saved in the parameterVec
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
  posFcn[1] = inter->CalcPrior(inter->GetParMap(min->X()));

  for(int i=0; i<parameterVec->size(); i++){
    nom[i] = parameterVec->at(i)->nominalValue;
    pre[i] = parameterVec->at(i)->prefitValue;
    pos[i] = parameterVec->at(i)->getPhysicalValue();

    posError[i] = inter->FitToPhysical(parameterVec->at(i)->name, errors[i]);
    posFitBasis[i] = bestFit[i];
    posErrorFitBasis[i] = errors[i];

    param_num = i;
    parameter_fit = bestFit[i];
    parameter_physical = inter->FitToPhysical(parameterVec->at(i)->name, bestFit[i]);
    error_fit = errors[i];
    error_physical = inter->FitToPhysical(parameterVec->at(i)->name, errors[i]);
    name = TString(parameterVec->at(i)->name);
    paramtree->Fill();
  }

  paramtree->Write();
  std::vector<std::array<double, 4> > allSSEMSimulation = inter->GetBeamProperties();
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
  std::string ssemDataFile="ssem_data/run0920332_gen.root";//"./ssem_data/run0910216_gen.root";


  auto allParameterVec = std::vector<std::shared_ptr<Parameter>>();
  //BPV1 has zero current so set it to fixed, but should add some freedom due to residual magnetisation in future
  //if we set them fixed here they will be permanently fixed no matter what we do in performFit
  allParameterVec.push_back(std::make_shared<Parameter>("BPV1", "fixed", "vBend", "preparation", 0));
  allParameterVec.push_back(std::make_shared<Parameter>("BPH2", "free", "hBend", "preparation", 1));
  allParameterVec.push_back(std::make_shared<Parameter>("QPQ1", "free", "quad", "preparation", 2));
  allParameterVec.push_back(std::make_shared<Parameter>("QPQ2", "free", "quad", "preparation", 4));
  allParameterVec.push_back(std::make_shared<Parameter>("BPD1", "free", "hBend", "preparation", 5));
  allParameterVec.push_back(std::make_shared<Parameter>("BPD2", "free", "hBend", "preparation", 6));
  allParameterVec.push_back(std::make_shared<Parameter>("QPQ3", "free", "quad", "preparation", 7));
  allParameterVec.push_back(std::make_shared<Parameter>("BPV2", "free", "vBend", "preparation", 8));
  allParameterVec.push_back(std::make_shared<Parameter>("QPQ4", "free", "quad", "preparation", 9));
  allParameterVec.push_back(std::make_shared<Parameter>("BPH3", "free", "hBend", "preparation", 10));
  allParameterVec.push_back(std::make_shared<Parameter>("QPQ5", "free", "quad", "preparation", 11));

  allParameterVec.push_back(std::make_shared<Parameter>("BAD1",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAF1",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAD2",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAF2",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAD3",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAF3",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAD4",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAF4",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAD5",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAF5",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAD6",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAF6",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAD7",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAF7",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAD8",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAF8",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAD9",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAF9",  "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAD10", "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAF10", "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAD11", "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAF11", "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAD12", "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAF12", "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAD13", "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAF13", "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAD14", "free", "arcBend", "arc"));
  allParameterVec.push_back(std::make_shared<Parameter>("BAF14", "free", "arcBend", "arc"));

  allParameterVec.push_back(std::make_shared<Parameter>("QFQ1",  "free", "quad", "finalFocus", 19));
  allParameterVec.push_back(std::make_shared<Parameter>("BFV1",  "free", "vBend", "finalFocus", 20));
  allParameterVec.push_back(std::make_shared<Parameter>("BFH1",  "free", "hBend", "finalFocus", 21));
  allParameterVec.push_back(std::make_shared<Parameter>("BFV2",  "free", "vBend", "finalFocus", 22));
  allParameterVec.push_back(std::make_shared<Parameter>("QFQ2",  "free", "quad", "finalFocus", 23));
  allParameterVec.push_back(std::make_shared<Parameter>("QFQ3",  "free", "quad", "finalFocus", 24));
  allParameterVec.push_back(std::make_shared<Parameter>("BFH2",  "free", "hBend", "finalFocus", 25));
  allParameterVec.push_back(std::make_shared<Parameter>("BFVD1", "free", "vBend", "finalFocus", 26));
  allParameterVec.push_back(std::make_shared<Parameter>("QFQ4",  "free", "quad", "finalFocus", 27));
  allParameterVec.push_back(std::make_shared<Parameter>("BFVD2", "free", "vBend", "finalFocus", 28));

  allParameterVec.push_back(std::make_shared<Parameter>("X0", "free", "beamPosition"));
  allParameterVec.push_back(std::make_shared<Parameter>("Xp0", "free", "beamPosition"));
  allParameterVec.push_back(std::make_shared<Parameter>("emitx", "free", "beamTwiss"));
  allParameterVec.push_back(std::make_shared<Parameter>("betx", "free", "beamTwiss"));
  allParameterVec.push_back(std::make_shared<Parameter>("alfx", "free", "beamTwiss"));
  allParameterVec.push_back(std::make_shared<Parameter>("Y0", "free", "beamPosition"));
  allParameterVec.push_back(std::make_shared<Parameter>("Yp0", "free", "beamPosition"));
  allParameterVec.push_back(std::make_shared<Parameter>("emity", "free", "beamTwiss"));
  allParameterVec.push_back(std::make_shared<Parameter>("bety", "free", "beamTwiss"));
  allParameterVec.push_back(std::make_shared<Parameter>("alfy", "free", "beamTwiss"));


  //now filter the parameters to only those we want to use in this optimisation (only necessary if you have a subset of the beamline)
  for(auto &par : allParameterVec){
    if(par->section == "preparation" || par->type.find("beam") != std::string::npos){ //just preparation section for testing
      parameterVec->push_back(par);
      (*parameterMap)[par->name] = par;
    }
  }

  int nPars = parameterVec->size();

  //can apply some filtering to parameters if we want but can do this at fit time, doesn't care about parameters that are unused

  Interface inter(ssemDataFile, baseBeamlineFile, parameterMap, parameterVec);
  inter.bds->SetFileWriting(false); //dont want to save a file for every simulation run, realllly slows things down


  //options for the initial magnet field strengths, all false means use ssem data file and estimates of magnet strengths
  char *usePrevBestFit = nullptr; //"./may_2025_hk_fits/misalignments_no_noise_tight_tolerance.root";
  bool useFieldMaps = false; //currently unsupported
  bool useInputFile = false;


  inter.SetInitialValues(usePrevBestFit, useFieldMaps, useInputFile, 0.0); //last arg is noise

  //set prior constraints
  //negative values are a fractional uncertainty on the nominal value
  //positive values are an absolute uncertainty in the same units as the parameter

//  inter.priorErrors["BPV1"] = 0.00001; //nominal current is 0 maybe there should actually be no freedom here..
  (*parameterMap)["BPH2"]->priorError = -0.07;
  (*parameterMap)["QPQ1"]->priorError = -0.07;
  (*parameterMap)["QPQ2"]->priorError = -0.07;
  (*parameterMap)["BPD1"]->priorError = -0.07;
  (*parameterMap)["BPD2"]->priorError = -0.07;
  (*parameterMap)["QPQ3"]->priorError = -0.07;
  (*parameterMap)["BPV2"]->priorError = -0.07;
  (*parameterMap)["QPQ4"]->priorError = -0.07;
  (*parameterMap)["BPH3"]->priorError = -0.07;
  (*parameterMap)["QPQ5"]->priorError = -0.07;

  //really lose constraints on beam parameters, but does really help with them not exploding 
  (*parameterMap)["X0"]->priorError = 2;
  (*parameterMap)["emitx"]->priorError = -0.2;
  (*parameterMap)["betx"]->priorError = 10;
  (*parameterMap)["alfx"]->priorError = 5;

  (*parameterMap)["Y0"]->priorError = 2;
  (*parameterMap)["emity"]->priorError = -0.2;
  (*parameterMap)["bety"]->priorError = 10;
  (*parameterMap)["alfy"]->priorError = 5;


  inter.SetChisqMode(1+2+4+8);
 
//just see how long a single iteration takes
//  auto iterstarttime = std::chrono::high_resolution_clock::now();
//  inter.fcn(pars);
//  auto iterendtime = std::chrono::high_resolution_clock::now();
//  auto itertime = std::chrono::duration_cast<std::chrono::microseconds>(iterendtime-iterstarttime).count();
//  std::cout<<"Took "<<itertime*1e-6<<"s to run a single iteration"<<std::endl;

//now for fitting

  std::map<std::string, std::vector<std::string>> freeGroups;
  freeGroups["all"] = {};
  freeGroups["B_PS"] = {};
  freeGroups["Q_PS"] = {};
  freeGroups["beamPosition"] = {};
  freeGroups["beamTwiss"] = {};

  for(auto &par : *parameterVec){
    freeGroups["all"].push_back(par->name);
    if((par->type == "hBend" || par->type == "vBend") && par->section == "preparation") freeGroups["B_PS"].push_back(par->name);
    if(par->type == "quad" && par->section == "preparation") freeGroups["Q_PS"].push_back(par->name);
    if(par->type == "beamPosition") freeGroups["beamPosition"].push_back(par->name);
    if(par->type == "beamTwiss") freeGroups["beamTwiss"].push_back(par->name);
  }


  auto wrappedFcn = [&inter](const double* pars) {
      return inter.fcn_wrapper(pars);
  };
 
  ROOT::Math::Functor f(wrappedFcn, nPars);
  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

  min->SetFunction(f);

//  TFile *result = new TFile("result.root", "RECREATE");
//  TH1D *hist;
//  for(int i=0; i<nPars; i++){
//    std::string name = parameterVec->at(i)->name;
//    double val = parameterVec->at(i)->getFitValue();
//    if(val>=0) hist = new TH1D(name.c_str(), name.c_str(), 11, val*0.5, val*1.5);
//    else hist = new TH1D(name.c_str(), name.c_str(), 11, val*1.5, val*0.5);
//
//    inter.ParamScan(parameterVec->at(i)->name, hist);
//
//    result->cd();
//    hist->Write(name.c_str());
//  }
//  result->Close();
//  exit(1);

  std::cout << "about to call perform fit with B magnets and beam position only"<<std::endl;
  performFit(min, &inter, {freeGroups["B_PS"], freeGroups["beamPosition"]}, 1+2);
  std::cout << "about to call perform fit with Q magnets and beam twiss parameters only"<<std::endl;
  performFit(min, &inter, {freeGroups["Q_PS"], freeGroups["beamTwiss"]}, 4+8);
  std::cout << "about to call perform fit with beam only"<<std::endl;
  performFit(min, &inter, {freeGroups["B_PS"], freeGroups["Q_PS"]}, 1+2+4+8);
 
  saveResult(min, &inter, nPars, "fit_results_after_first_split_optimisation.root");

  std::cout << "about to call perform fit with B magnets only take 2"<<std::endl;
  performFit(min, &inter, {freeGroups["B_PS"], freeGroups["beamPosition"]}, 1+2);
  std::cout << "about to call perform fit with Q magnets only take 2"<<std::endl;
  performFit(min, &inter, {freeGroups["Q_PS"], freeGroups["beamTwiss"]}, 4+8);
  std::cout << "about to call perform fit with beam only take 2"<<std::endl;
  performFit(min, &inter, {freeGroups["B_PS"], freeGroups["Q_PS"]}, 1+2+4+8);
 
  saveResult(min, &inter, nPars, "fit_results_after_second_split_optimisation.root");

  std::cout << "about to call perform fit with all parameters free (this may take a while)"<<std::endl;
  performFit(min, &inter, {freeGroups["all"]}, 1+2+4+8);


  saveResult(min, &inter, nPars, "fit_results.root");
  if(argc > 1) saveResult(min, &inter, nPars, argv[1]);

  inter.GenerateInputFile(min->X());

  auto endtime = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::microseconds>(endtime-starttime).count();
  std::cout<<"Took "<<time*1e-6<<"s to run"<<std::endl;
 
  return 0;

}
