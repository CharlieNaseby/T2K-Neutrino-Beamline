#include "Interface.h"
#include "Optics.h" 

//constructor to set data and the gmad file we'll use as a base to our fit

Interface::Interface(std::vector<std::array<double, 4> > data, std::string baseBeamlineFile){
  dat.resize(data.size());
  for(int i=0; i<dat.size(); i++){
    for(int j=0; j<4; j++) dat[i][j] = data[i][j];
  }
  ParseInputFile(baseBeamlineFile);
};

Interface::~Interface(){};

void Interface::SetNPars(int npars){
  if(npars != beamline.size()-1){
    std::cerr<<"Number of parameters set is not equal to the number in the gmad file supplied FILE "<<__FILE__<<":"<<__LINE__<<std::endl;
    throw  ;
  }
  internalPars.resize(npars);
  nPars = npars;
}

void Interface::ParamScan(int param, TH1D *hist){
  double store = internalPars[param];
  for(int i=1; i<=hist->GetNbinsX(); i++){
    double xpoint = hist->GetBinCenter(i);
    internalPars[param] = xpoint;
    hist->SetBinContent(i, fcn(internalPars));
  }
  internalPars[param] = store;
}

void Interface::ParamScan(int param, int npoints, double xlo, double xup){
  double step = (xup-xlo)/(double)(npoints-1);
  double store = internalPars[param];
  internalPars[param] = xlo;
  for(int i=0; i<npoints; i++){
    std::cout << "Param " << param << "=" << internalPars[param]<< " chisq = " << fcn(internalPars)<<std::endl;;
    internalPars[param]+=step;
  }
  internalPars[param] = store;
}

void Interface::SetInternalPars(const double *pars){
  for(int i=0; i<nPars; i++) internalPars[i] = pars[i];
}

//the cost function we'll minimise

double Interface::fcn(std::vector<double> pars){

  const double *tmp = pars.data();
  return fcn(tmp);
}
double Interface::fcn_wrapper(const double *pars){
  return fcn(pars);
}

double Interface::fcn(const double *pars){
  SetInternalPars(pars);
  GenerateInputFile(pars);
  return CalcChisq();
}

//I hate that this is the way to do this just as much as you do


//first need to write a gmad input file with the magnet strengths implied by pars
void Interface::GenerateInputFile(const double *pars){

  std::fstream out;
  out.open("../gmad/optimised.gmad", std::ios::out);
  for(int i=0; i<beamline.size()-1; i++){
    out<<beamline[i].c_str()<<std::setprecision(14)<<pars[i];
  }
  out<<beamline[beamline.size()-1];
  out.close();
}


//read an input file, splitting the file every time theres the string "bScaling=" 
//this indicates the start of a variable of out fit

void Interface::ParseInputFile(std::string baseBeamlineFile){
  std::ifstream infile;
  infile.open(baseBeamlineFile);
  beamline.resize(1);
  std::string line;

  while(std::getline(infile, line)){
    //does it have the string bScaling in it?
    size_t strpos = line.find("bScaling=");
    line.append("\n");
    if(strpos == std::string::npos) beamline[beamline.size()-1].append(line);
    else{
      std::string firstPart = line.substr(0, strpos);
      firstPart.append("bScaling=");
      beamline[beamline.size()-1].append(firstPart);
      std::string secondPart = line.substr(strpos);
      size_t comma = secondPart.find_first_of(",");
      std::string newSecondPart = secondPart.substr(comma);
      
      beamline.push_back(newSecondPart);
    }
  }
}

double Interface::CalcChisq(){

  std::system("bdsim --file=../gmad/optimised.gmad --batch --ngenerate=50 --outfile=/home/bdsim_output --seed=1989 > /dev/null");
  std::system("rebdsimOptics /home/bdsim_output.root /home/bdsim_output_optics.root > /dev/null");

  Optics beamOptics("/home/bdsim_output_optics.root");

  double chisq=0;
  for(int i=0; i<dat.size(); i++){
    beamOptics.fChain->GetEntry(i+1);
    std::array<double, 4> simulation = {beamOptics.Mean_x, beamOptics.Mean_y, beamOptics.Sigma_x, beamOptics.Sigma_y};
    for(int j=0; j<2; j++) chisq += simulation[j]*simulation[j]/(0.2/1000.);  //CERN TODO, just position first
//    std::cout<<"SSEM "<<i+1<<" position cumulative chisq contribution "<<chisq<<std::endl;
  }
  return chisq;
}

