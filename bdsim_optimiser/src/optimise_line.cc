#include "Interface.h"
#include "TClassTable.h"
#include <chrono>
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TMatrixD.h"

int main(int argc, char **argv){
  auto starttime = std::chrono::high_resolution_clock::now();

  std::string baseBeamlineFile="../gmad/test.gmad";
  std::string ssemDataFile="./ssem_data/run0910580_gen.root";

  const int nPars = 11;
  Interface inter(ssemDataFile, baseBeamlineFile, nPars);

  double pars[nPars];

  bool usePrevBestFit = false;
  bool useFieldMaps = false;


  inter.SetInitialValues(usePrevBestFit, useFieldMaps, pars);
  inter.SetInternalPars(pars);
  inter.SetNominalPars(pars);

  for(int i=0; i<nPars; i++) std::cout<<inter.magNames[i]<<"\t"<<pars[i]<<std::endl;

  inter.SetChisqMode(3);

  auto iterstarttime = std::chrono::high_resolution_clock::now();
  inter.fcn(pars);
  auto iterendtime = std::chrono::high_resolution_clock::now();
  auto itertime = std::chrono::duration_cast<std::chrono::microseconds>(iterendtime-iterstarttime).count();
  std::cout<<"Took "<<itertime*1e-6<<"s to run a single iteration"<<std::endl;
//  exit(1);

  ROOT::Math::Functor f(inter, &Interface::fcn_wrapper, nPars);
  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  min->SetStrategy(3);
  min->SetFunction(f);
  min->SetMaxFunctionCalls(10000);
  min->SetTolerance(10);

  for(int i=0; i<nPars; i++) min->SetVariable(i, inter.magNames[i], pars[i], 0.1);

  std::vector<int> bMagVars = {0, 1, 4, 5, 7, 9};
  std::vector<int> qMagVars = {2, 3, 6, 8, 10};

//  for(auto fixed : qMagVars) min->FixVariable(fixed);

//  for(int i=0; i<nPars; i++) min->FixVariable(i);
  min->SetPrintLevel(2);

//  min->ReleaseVariable(0); 
//  min->ReleaseVariable(1);
//  min->ReleaseVariable(4);
//  min->ReleaseVariable(7);
  min->FixVariable(9); //BPH3 and QPQ5 do not have and SSEMs afterwards in the sim so would have no data constraint
  min->FixVariable(10);
  min->Minimize();
  min->PrintResults();



  TFile *outf = new TFile("fit_results.root", "RECREATE");

  double *covarray = new double[nPars*nPars];
  min->GetCovMatrix(covarray);
  TMatrixD *cov = new TMatrixD(nPars, nPars, covarray);
  cov->Write("postfit_covariance");
  const double *errors = min->Errors();
  const double *bestFit = min->X();
  TH1D *fitResult = new TH1D("fit_result", "fit_result", nPars, 0, nPars);
  TH1D *prefit = new TH1D("prefit", "prefit", nPars, 0, nPars);
  TH1D *post_minus_pred_div_pred = new TH1D("post_minus_pred_div_pred", "(Post - Pre)/Pre", nPars, 0, nPars);
  const double *nominal = inter.GetNominalPars();
  for(int i=0; i<nPars; i++){
    prefit->SetBinContent(i+1, inter.preFit[i]);
    fitResult->SetBinContent(i+1, bestFit[i]);
    fitResult->SetBinError(i+1, errors[i]);
    if(TMath::Abs(inter.preFit[i])>1e-6){
      post_minus_pred_div_pred->SetBinContent(i+1, (bestFit[i]-inter.preFit[i])/inter.preFit[i]);
      post_minus_pred_div_pred->SetBinError(i+1, errors[i]/inter.preFit[i]);
    }
   
  }
  fitResult->Write("best_fit");
  prefit->Write("prefit");
  post_minus_pred_div_pred->Write("post_minus_pre_div_pre");
  outf->Close();

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
