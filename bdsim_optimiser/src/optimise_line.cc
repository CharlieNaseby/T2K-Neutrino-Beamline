#include "Interface.h"
#include "TClassTable.h"
#define NSSEM 9

int main(){

  std::vector<std::array<double, 4> >  data;
  data.resize(NSSEM);
  for(int i=0; i<NSSEM; i++) data[i] = {0, 0, 0, 0};

  std::string baseBeamlineFile="../gmad/test.gmad";
  Interface inter(data, baseBeamlineFile);
  const double pars[2] = {1/7., 1};
  inter.SetNPars(1);
  //scan a parameter

  inter.SetInternalPars(pars);

  inter.ParamScan(0, 5, 0.5, 1.5);

  return 0;

}
