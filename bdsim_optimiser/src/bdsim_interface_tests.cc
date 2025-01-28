#include "Interface.h"
#include "TClassTable.h"
#include <chrono>
#include "Math/Minimizer.h"
#include "Math/Functor.h"
#include "Math/Factory.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "BDSIMClass.hh"


int main(int argc, char **argv){
  auto starttime = std::chrono::high_resolution_clock::now();

  std::string baseBeamlineFile="../gmad/test.gmad";
  std::string ssemDataFile="./ssem_data/run0910216_gen.root";

  const int nMagnetPars = 11;
  const int nBeamPars = 10;
  const int nPars = nMagnetPars + nBeamPars;

  Interface inter(ssemDataFile, baseBeamlineFile, nPars, nMagnetPars, nBeamPars);

 
//just see how long a single iteration takes
 
  auto iterstarttime = std::chrono::high_resolution_clock::now();
  char *dargv[6];
  char path[256] = "--file=/home/T2K-Neutrino-Beamline/gmad/test.gmad";
  char batch[256] = "--batch";
  char ngen[256] = "--ngenerate=200";
  char outfile[256] = "--outfile=/home/bdsim_output";
  char seed[256] = "--seed=1989";
  dargv[1] = path;
  dargv[2] = batch;
  dargv[3] = ngen;
  dargv[4] = outfile;
  dargv[5] = seed;


  inter.bds->Initialise(5, dargv);
  std::cout<<"Initialised!!\n\n"<<std::endl;
  inter.bds->BeamOn();
  std::cout<<"Beam on done!!\n\n"<<std::endl;
  inter.bds->BeamOn();
  std::cout<<"Beam on done!!\n\n"<<std::endl;



  auto iterendtime = std::chrono::high_resolution_clock::now();
  auto itertime = std::chrono::duration_cast<std::chrono::microseconds>(iterendtime-iterstarttime).count();
  std::cout<<"Took "<<itertime*1e-6<<"s to run a single iteration"<<std::endl;
 
 
  auto endtime = std::chrono::high_resolution_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::microseconds>(endtime-starttime).count();
  std::cout<<"Took "<<time*1e-6<<"s to run"<<std::endl;
 
  return 0;

}

