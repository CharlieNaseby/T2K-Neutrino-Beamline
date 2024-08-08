#include <iostream>
#include "Interface.h"

int main(){

  std::cout<<"Hello World!"<<std::endl;
  double data[2] = {3, 5};
  Interface inter(data);
  const double pars[2] = {1, 1};
  inter.calc_chisq(pars);

  return 0;

}

