#include "Interface.h"
#include "TClassTable.h"
#include <chrono>
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "BDSIMClass.hh"
#include <gperftools/profiler.h>

int main(){
  auto starttime = std::chrono::high_resolution_clock::now();

  std::string baseBeamlineFile="../survey/unoptimised.gmad";
  std::string ssemDataFile="./ssem_data/run0910216_gen.root";

  const int nMagnetPars = 11;
  const int nBeamPars = 10;
  const int nPars = nMagnetPars + nBeamPars;

  Interface inter(ssemDataFile, baseBeamlineFile, nPars, nMagnetPars, nBeamPars);

  inter.bds->SetFileWriting(false); //dont want to save a file for every simulation run, realllly slows things down

  double pars_array[nPars];
  inter.SetInitialValues(nullptr, false, false, false, pars_array);
  inter.SetInternalPars(pars_array);
  inter.SetChisqMode(1+2+4+8);

  double ones[nPars];
  for(int i=0; i<nPars; i++) ones[i] = 1.0;


  inter.fcn(ones);
//  inter2.fcn(ones);


//just see how long a single iteration takes
 
  auto iterstarttime = std::chrono::high_resolution_clock::now();
  std::cout<<"Initialised!!\n\n"<<std::endl;
  std::map<std::string, double> pars;
  pars["BPD1"] = -1.09;
//  pars["QPQ1"] = 0.055;

  inter.bds->BeamOn(100, pars);
  std::vector<std::array<double, 4> > simResult = inter.bds->CalcBeamPars();
  for(auto axis : simResult[8]) std::cout << "resulting beam parameters at ssem9 " << axis << std::endl;

  pars["BPD1"] = -1.11;
  inter.bds->BeamOn(100, pars);
  simResult = inter.bds->CalcBeamPars();
  for(auto axis : simResult[8]) std::cout << "resulting beam parameters at ssem9 " << axis << std::endl;

  pars["BPD1"] = -1.15;
  inter.bds->BeamOn(100, pars);
  simResult = inter.bds->CalcBeamPars();
  for(auto axis : simResult[8]) std::cout << "resulting beam parameters at ssem9 " << axis << std::endl;

//  inter2.bds->BeamOn();
  std::cout<<"Beam on done!!\n\n"<<std::endl;

  std::vector<int> times;


  for(int i=0; i<2000; i++){
    if(i==0) ProfilerStart("profile_first_200.prof");
    if(i==200) ProfilerStop();

    if(i==1799) ProfilerStart("profile_1800_to_2000.prof");
    if(i==1999) ProfilerStop();

//    auto loopstarttime = std::chrono::high_resolution_clock::now();
    inter.fcn(ones);
    //inter.bds->BeamOn(100, pars);
    //std::vector<std::array<double, 4> > simResult = inter.bds->CalcBeamPars();
    //for(auto axis : simResult[8]) std::cout << "resulting beam parameters at ssem9 " << axis << std::endl;
//    auto loopendtime = std::chrono::high_resolution_clock::now();
//    auto looptime = std::chrono::duration_cast<std::chrono::microseconds>(loopendtime-loopstarttime).count();
//    std::cout<<"\n"<<looptime<<"\n"<<std::endl;

//    times.push_back(looptime);
  }  
//  std::cout<<"Beam on done!!\n\n"<<std::endl;



  for(auto t : times) std::cout<<t*1e-6<<std::endl;

  auto iterendtime = std::chrono::high_resolution_clock::now();
  auto itertime = std::chrono::duration_cast<std::chrono::microseconds>(iterendtime-iterstarttime).count();
  std::cout<<"Took "<<itertime*1e-6<<"s to run a single iteration"<<std::endl;
 
 
  auto endtime = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::microseconds>(endtime-starttime).count();
  std::cout<<"Took "<<time*1e-6<<"s to run"<<std::endl;

  return 0;

}

