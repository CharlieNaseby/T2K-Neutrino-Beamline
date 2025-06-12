#include "Interface.h"
#include "Optics.h" 

//constructor to set data and the gmad file we'll use as a base to our fit

//very simple setup for testing purposes and bayesian optimisation studies
Interface::Interface(std::vector<double> pars, std::vector<std::array<double, 4> > target){
  nPars = pars.size();
  nMagnetPars = pars.size();
  nBeamPars = 0;

  bds = new CNBDSIM();

  char *dargv[6];
  char path[256] = "--file=simple_beamline.gmad";
  char batch[256] = "--batch";
  char ngen[256] = "--ngenerate=100";
  char outfile[256] = "--outfile=/opt/bdsim_output";
  char seed[256] = "--seed=1989";
  //char verbose[256] = "--verbose=0";
  dargv[1] = path;
  dargv[2] = batch;
  dargv[3] = ngen;
  dargv[4] = outfile;
  dargv[5] = seed;
//  dargv[6] = verbose;
  bds->Initialise(5, dargv);
  internalPars.resize(nPars); //point about which LLH scans will be conducted
  magNames.resize(nMagnetPars);
  beamNames.resize(nBeamPars);
  preFit.resize(nPars); //the values of the parameters expected with scalings or fudge factors applied

  SetData(target);




  magMap["B1"] = 0;
  magMap["Q1"] = 1;
//  magMap["B1"] = 2;
//  magMap["Q3"] = 3;
//  magMap["B2"] = 4;


  for(auto [key, value] : magMap) magNames[value] = key;
  for(auto name : magNames) parNames.push_back(name);

  ParseInputFile("./simple_beamline.gmad");

}


//the full T2K beamline setup
Interface::Interface(std::string dataFile, std::string baseBeamlineFile, int npars, int nmagnetpars, int nbeampars){

  nPars = npars;
  nMagnetPars = nmagnetpars;
  nBeamPars = nbeampars;
  bds = new CNBDSIM();

  char *dargv[6];
  char path[256] = "--file=../survey/unoptimised.gmad";
  char batch[256] = "--batch";
  char ngen[256] = "--ngenerate=100";
  char outfile[256] = "--outfile=/opt/bdsim_output";
  char seed[256] = "--seed=1989";
  //char verbose[256] = "--verbose=0";
  dargv[1] = path;
  dargv[2] = batch;
  dargv[3] = ngen;
  dargv[4] = outfile;
  dargv[5] = seed;
//  dargv[6] = verbose;
  bds->Initialise(5, dargv);

  internalPars.resize(nPars); //point about which LLH scans will be conducted
//  nominalPars.resize(nPars); //the values of the parameters expected based on true magnet current
  magNames.resize(nMagnetPars);
  beamNames.resize(nBeamPars);
  preFit.resize(nPars); //the values of the parameters expected with scalings or fudge factors applied

  SetSSEMData(dataFile);


  magMap["BPV1"] = 0;
  magMap["BPH2"] = 1;
  magMap["QPQ1"] = 2;
  magMap["QPQ2"] = 3;
  magMap["BPD1"] = 4;
  magMap["BPD2"] = 5;
  magMap["QPQ3"] = 6;
  magMap["BPV2"] = 7;
  magMap["QPQ4"] = 8;
  magMap["BPH3"] = 9;
  magMap["QPQ5"] = 10;

  for(auto [key, value] : magMap) magNames[value] = key;

  beamNames ={"X0", "Xp0", "emitx", "betx", "alfx", "Y0", "Yp0", "emity", "bety", "alfy"};

  for(auto name : magNames) parNames.push_back(name);
  for(auto name : beamNames) parNames.push_back(name);

  ParseInputFile(baseBeamlineFile);
}

Interface::~Interface(){
  delete bds;
}

void Interface::SetSSEMData(std::string dataFile){
  TFile *ssemDataFile = new TFile(dataFile.c_str(), "READ");
  TTree *ssemData = dynamic_cast<TTree*> (ssemDataFile->Get("anabeam"));
  double ssemx[19], ssemax[19], ssemy[19], ssemwx[19], ssemwy[19], ct[5], magset[31];
  ssemData->SetBranchAddress("magset", magset);
  ssemData->SetBranchAddress("ssemx", ssemx);
  ssemData->SetBranchAddress("ssemax", ssemax);
  ssemData->SetBranchAddress("ssemy", ssemy);
  ssemData->SetBranchAddress("ssemwx", ssemwx);
  ssemData->SetBranchAddress("ssemwy", ssemwy);
  ssemData->SetBranchAddress("ct", ct);

  ssemData->GetEntry(0);
  for(int i=0; i<31; i++) magCurrent.push_back(magset[i]);

  std::vector<double> ssemxMean(NSSEM), ssemyMean(NSSEM), ssemwxMean(NSSEM), ssemwyMean(NSSEM);
  for(int i=0; i<NSSEM; i++){
    ssemxMean[i] = 0;
    ssemyMean[i] = 0;
    ssemwxMean[i] = 0;
    ssemwyMean[i] = 0;
  }
  int nShots=0;


  for(int i=0; i<ssemData->GetEntries(); i++){
    ssemData->GetEntry(i);
    if(ssemax[0] < 100) continue;  //if there was no beam then skip sadly for run910216 the CTs were disabled?
    for(int j=0; j<11; j++){ //CERN for some reason the later magnets do change from shot-to-shot
      if(magset[j] != magCurrent[j]) goto loopend; //select only cases with the same current each shot
                                                   //actually legit use of a goto!
//      std::cout<<j<<" magset "<<magset[j]<<" magcurrent "<<magCurrent[j]<<std::endl;
    }
    for(int j=0; j<NSSEM; j++){
      ssemxMean[j]+=ssemx[j];
      ssemyMean[j]+=ssemy[j];
      ssemwxMean[j]+=ssemwx[j];
      ssemwyMean[j]+=ssemwy[j];
    }
    nShots++;
loopend:;
  }
  for(int i=0; i<NSSEM; i++){
    ssemxMean[i] /= (double)nShots;
    ssemyMean[i] /= (double)nShots;
    ssemwxMean[i] /= (double)nShots;
    ssemwyMean[i] /= (double)nShots;
  }
  SetData(ssemxMean, ssemyMean, ssemwxMean, ssemwyMean);
}


void Interface::SetData(std::vector<std::array<double, 4> > target){
  std::vector<double> x, y, wx, wy;
  for(auto ssem : target){
    x.push_back(ssem[0]);
    y.push_back(ssem[1]);
    wx.push_back(ssem[2]);
    wy.push_back(ssem[3]);
  }
  SetData(x, y, wx, wy);
}


void Interface::SetData(std::vector<double> x, std::vector<double> y, std::vector<double> wx, std::vector<double> wy){
  dat.resize(x.size());
  for(unsigned int i=0; i<x.size(); i++){
    dat[i][0] = x[i];
    dat[i][1] = y[i];
    dat[i][2] = wx[i];
    dat[i][3] = wy[i];
  }
}

void Interface::SetInitialValues(char *usePrevBestFit, bool useFieldMaps, bool useFudgeFactor, bool useInputFile, double *pars, double noise){

  std::map<std::string, double> beamPars;
  beamPars["X0"]=-0.5e-3; //units of m
  beamPars["Xp0"]=3.5e-5;
  beamPars["emitx"]=0.0610768e-6; //m rad
  beamPars["betx"]=37.098; //m
  beamPars["alfx"]=-2.4187;
  beamPars["dispx"]=0.42373,
  beamPars["dispxp"]=0.07196,
  beamPars["Y0"]=-0.2e-3;
  beamPars["Yp0"]=7.8e-5;
  beamPars["emity"]=0.05976e-6;
  beamPars["bety"]=5.45;
  beamPars["alfy"]=0.178;
  beamPars["dispy"]=0.;
  beamPars["dispyp"]=0.0;


  std::map<std::string, double> kScaling;
  //derived from the "KIcurve sad file, scaling of K1 or K0 value in units of 100A
  kScaling["BPV1"] = -0.08760146181454302;
  kScaling["BPH2"] = -0.06372506685790527;
  kScaling["QPQ1"] = 0.010747477;
  kScaling["QPQ2"] = -0.012669113;
  kScaling["BPD1"] = -0.0926275917278186;
  kScaling["BPD2"] = -0.08940591194129406;
  kScaling["QPQ3"] = 0.01593171666666667;
  kScaling["BPV2"] = -0.12238521642828552;
  kScaling["QPQ4"] = -0.01539898;
  kScaling["BPH3"] = -0.11826754864042695;
  kScaling["QPQ5"] = 0.01768722;


  //values based on the fieldmap fields at 100A
  //note that for quads these strengths are 1/(B*rho) * dB/dx
  //for 30 GeV KE protons gamma*m_p*c^2 = 30.938 GeV
  //so B*rho = 30.938*beta/proton charge
  //=30.938 * 0.9995/0.30286   elementary charge in natural units from alpha=1/137
  //(1/103.101)*dB/dx
  kScaling["BPV1"] = -0.09456593004136546;

  kScaling["QPQ1"] = 0.010517526;
  kScaling["QPQ2"] = -0.01261134;
  kScaling["BPD1"] = -0.09289046446110355;
  kScaling["BPD2"] = -0.09289046446110355;
  kScaling["QPQ3"] = 0.01574514;
  kScaling["BPV2"] = -0.12654858254158613;
  kScaling["QPQ4"] = -0.01471238222852117;
  kScaling["BPH3"] = -0.12432151082458207;
  kScaling["QPQ5"] = 0.019445574477550798;


  for(auto &val : magCurrent) val += 1e-1; //annoyingly magnets with 0 strength create issues

  nominalPars["BPV1"] = magCurrent[0]/100.;
  nominalPars["BPH2"] = magCurrent[1]/100.;
  nominalPars["QPQ1"] = magCurrent[2]/100.;
  nominalPars["QPQ2"] = magCurrent[4]/100.;
  nominalPars["BPD1"] = magCurrent[5]/100.;
  nominalPars["BPD2"] = magCurrent[6]/100.;
  nominalPars["QPQ3"] = magCurrent[7]/100.;
  nominalPars["BPV2"] = magCurrent[8]/100.;
  nominalPars["QPQ4"] = magCurrent[9]/100.;
  nominalPars["BPH3"] = magCurrent[10]/100.;
  nominalPars["QPQ5"] = magCurrent[11]/100.;

  std::vector<double> fudgeFactor(nMagnetPars);

  fudgeFactor[0] = 1;
  fudgeFactor[1] = 1;
  fudgeFactor[2] = 1;
  fudgeFactor[3] = 1;
  fudgeFactor[4] = 1;
  fudgeFactor[5] = 1;
  fudgeFactor[6] = 1;
  fudgeFactor[7] = 1;
  fudgeFactor[8] = 1;
  fudgeFactor[9] = 1;
  fudgeFactor[10] = 1;

  for(auto mag : magNames){
    nominalPars[mag] *= kScaling[mag];//nominal is always the expected value for the parameters based on the currents
    preFit[magMap[mag]] = nominalPars[mag];
    if(useFudgeFactor) preFit[magMap[mag]] *= fudgeFactor[magMap[mag]];
  }

  if(useFieldMaps){
    std::cerr << "Using field maps for the fit CURRENTLY UNSUPPORTED" << std::endl;
    throw;
    for(int i=0; i<nMagnetPars; i++) fudgeFactor[i] = 1.0;
  
    pars[0] = magCurrent[0]/100.;
    pars[1] = -magCurrent[1]/100.;
    pars[2] = magCurrent[2]/100.;
    pars[3] = magCurrent[4]/100.;
    pars[4] = magCurrent[5]/100.;
    pars[5] = magCurrent[6]/100.;
    pars[6] = magCurrent[7]/100.;
    pars[7] = -magCurrent[8]/100.;
    pars[8] = -magCurrent[9]/100.;
    pars[9] = magCurrent[10]/100.;
    pars[10] = magCurrent[11]/100.;
    for(int i=0; i<nMagnetPars; i++){
      nominalPars[magNames[i]] = pars[i];
      if(useFudgeFactor) pars[i] *= fudgeFactor[i];
    }
    for(int i=0; i<nMagnetPars; i++) preFit[i] = pars[i];
  }
  else{
    for(int i=0; i<nMagnetPars; i++) pars[i] = preFit[i];
  }

  //now for the beam parameters
  for(int i=0; i<nBeamPars; i++){
    pars[i+nMagnetPars] = beamPars[beamNames[i]];
    preFit[i+nMagnetPars] = beamPars[beamNames[i]];
    nominalPars[beamNames[i]] = beamPars[beamNames[i]];
  }
  if(usePrevBestFit){
    std::cout << "Using previous fit result contained in file previous_fit.root" << std::endl;
    TFile inf(usePrevBestFit, "READ");
    TVectorT<double> filePostFit = *(TVectorT<double>*)inf.Get("postFit");
    TVectorT<double> filePreFit = *(TVectorT<double>*)inf.Get("preFit");
    TVectorT<double> fileNominal = *(TVectorT<double>*)inf.Get("nominal");
   
    if(filePostFit.GetNrows() != nPars){
      std::cerr << "invalid size of parameters supplied from previous fit file"<<std::endl;
      throw;
    }
    for(int i=0; i<filePostFit.GetNrows(); i++){
      pars[i] = filePostFit[i];
      preFit[i] = filePreFit[i];
      nominalPars[parNames[i]] = fileNominal[i];
    }
    inf.Close();
  }
  else if(useInputFile){ //use the values from the gmad file supplied
    for(int i=0; i<nPars; i++){
      pars[i] = bds->GetParameterValue(parNames[i]);
      preFit[i] = pars[i];
      nominalPars[parNames[i]] = pars[i];
    }
  }


  TRandom3 rand(1989);

  for(int i=1; i<nPars; i++){ //CERN TODO start at 1 to avoid changing BPV1, really needs a way of not throwing fixed values
    pars[i]*=rand.Gaus(1.0, noise);
    std::cout << parNames[i] << " is set to " << pars[i] << std::endl;
  }

}
//simple one for just setting them outside of the T2K construct
void Interface::SetInitialValues(std::vector<double> init){

  for(int i=0; i<init.size(); i++){
    preFit[i] = init[i];
    nominalPars[parNames[i]] = init[i];
  }
}

void Interface::SetFileWriting(bool write){
  bds->SetFileWriting(write);
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

std::map<std::string, double> Interface::GetNominalPars(){
  return nominalPars;
}
//the cost function we'll minimise

double Interface::fcn(std::vector<double> pars){

  const double *tmp = pars.data();
  return fcn(tmp);
}
//need a non-overloaded function that just takes const double*
double Interface::fcn_wrapper(const double *pars){
  return fcn(pars);
}

bool Interface::CheckBounds(std::map<std::string, double> pars){
  std::vector<std::string> limitedpars = {"emitx", "betx", "emity", "bety"};
  for(auto st : limitedpars){
    if(pars[st] < 0){
      std::cerr<<"found parameter outside of range "<<st<<" with value "<<pars[st]<<std::endl; 
      return false;
    }
  }
  return true;
}


double Interface::fcn(const double *pars){
  SetInternalPars(pars);
  return CalcChisq(pars);
}

double Interface::CalcPrior(std::map<std::string, double> pars){
  double chisq = 0;
  for(auto [key, value] : priorErrors){
    double contribution = 0;
    if(value >= 0) contribution = std::pow(pars[key] - nominalPars[key], 2)/(std::pow(value,2)); //before you efficiency freaks complain pow(x, 2) expands to x*x with -O1+
    else contribution = std::pow(pars[key] - nominalPars[key], 2)/(std::pow(nominalPars[key]*value, 2));
    chisq += contribution;
//    std::cout<<"adding to the prior "<<key << " diff "<<pars[key]-nominalPars[key]<<" adding "<<contribution<<" to chisq"<<std::endl;
  }
  return chisq;
}

//if we want to write a gmad input file with the best-fit parameters
void Interface::GenerateInputFile(const double *pars){

  std::fstream out;
  out.open("../survey/optimised.gmad", std::ios::out);
  for(unsigned int i=0; i<beamline.size()-1; i++){
    out<<beamline[i].c_str()<<std::setprecision(14)<<FitToPhysical(i, pars[i]);
  }
  out<<beamline[beamline.size()-1];
  out.close();
}


//read an input file, splitting the file every time theres a string indicative of a parameter 

void Interface::ParseInputFile(std::string baseBeamlineFile){
  std::ifstream infile;
  infile.open(baseBeamlineFile);
  beamline.resize(1);
  std::string line;
  std::vector<std::string> matchStrings = {"bScaling=", "B=", "k1="};
  for(auto str : beamNames) matchStrings.push_back(str+'=');


  while(std::getline(infile, line)){
    size_t labstart = std::string::npos;
    size_t labend = std::string::npos;
    for(auto name : matchStrings){
      labstart = line.find(name);
      if(labstart != std::string::npos){
        labend = labstart + name.length();
        break;
      }
    }

    line.append("\n");
    if(labstart == std::string::npos){
      beamline[beamline.size()-1].append(line);
      continue;
    }

    std::string firstPart = line.substr(0, labend);
    beamline[beamline.size()-1].append(firstPart);

    std::string secondPart = line.substr(labend);
    size_t endidx = secondPart.find_first_of("*");  //grab the units
    size_t comma = secondPart.find_first_of(","); //if its unitless
    if(endidx == std::string::npos) endidx = comma;
    std::string newSecondPart = secondPart.substr(endidx);
    beamline.push_back(newSecondPart);

  }

  if(nPars != beamline.size()-1){
    std::cerr<<"Number of parameters set is not equal to the number in the gmad file supplied "<<__FILE__<<":"<<__LINE__<<std::endl;
    std::cerr<<"Expected "<<nPars<<" but file has "<<beamline.size()-1<<std::endl;
    throw  ;
  }
}

std::vector<std::array<double, 4> > Interface::GetBeamProperties(){
  return bds->CalcBeamProperties();
}


void Interface::TestBdsim(){
  bds->BeamOn();
}

std::vector<double> Interface::FitToPhysical(double *fitval){
  std::vector<double> retval;
  for(unsigned int i=0; i<parNames.size(); i++){
    retval.push_back(FitToPhysical(i, fitval[i]));
  }
  return retval;
}

double Interface::FitToPhysical(int i, double fitval){
  return fitval*preFit[i];
}

std::vector<double> Interface::PhysicalToFit(double *physval){
  std::vector<double> retval;
  for(unsigned int i=0; i<parNames.size(); i++){
    retval.push_back(PhysicalToFit(i, physval[i]));
  }
  return retval;
}

double Interface::PhysicalToFit(int i, double physval){
  return physval/preFit[i];
}

std::map<std::string, double> Interface::GetParMap(const double *pars){
  std::map<std::string, double> parmap;
  int i=0;
  for(auto key : parNames){
    parmap[key] = FitToPhysical(i, pars[i]);
    i++;
  }
  return parmap;
}

void Interface::BeamOn(int n, std::map<std::string, double> parmap){
  bds->BeamOn(100, parmap);
}

double Interface::CalcChisq(const double *pars){
  std::map<std::string, double> parmap = GetParMap(pars);
  if(CheckBounds(parmap) == false) return 1234567891.0;
  BeamOn(100, parmap);

  std::vector<std::array<double, 4> > allSSEMSimulation = GetBeamProperties();

  double chisq = 0;
  double chisqx = 0;
  double chisqy = 0;
  double chisqwx = 0;
  double chisqwy = 0;

  for(unsigned int i=0; i<dat.size(); i++){
    std::array<double, 4> simulation = allSSEMSimulation[i];
//    std::cout<<"SSEM"<<i+1<<" beam sim postion = \t"<<simulation[0]<<", \t"<<simulation[1]<<" data \t"<<dat[i][0]<<", \t"<<dat[i][1]<<std::endl;
//    std::cout<<"SSEM"<<i+1<<" beam sim width   = \t"<<simulation[2]<<", \t"<<simulation[3]<<" data \t"<<dat[i][2]<<", \t"<<dat[i][3]<<std::endl;
    chisqx += (dat[i][0]-simulation[0])*(dat[i][0]-simulation[0])/(0.2*0.2);  //CERN 0.2mm uncert on ssem position x
    chisqy += (dat[i][1]-simulation[1])*(dat[i][1]-simulation[1])/(0.2*0.2);  //position y
    chisqwx += (dat[i][2]-simulation[2])*(dat[i][2]-simulation[2])/(0.2*0.2);  //CERN width with 0.2mm precision x
    chisqwy += (dat[i][3]-simulation[3])*(dat[i][3]-simulation[3])/(0.2*0.2);  //width y
//    std::cout<<"SSEM "<<i+1<<" position cumulative chisq contribution "<<chisq<<std::endl;
  }

  double prior = CalcPrior(parmap);
  if(fitMode & 0x01) chisq += chisqx;
  if(fitMode & 0x02) chisq += chisqy;
  if(fitMode & 0x04) chisq += chisqwx;
  if(fitMode & 0x08) chisq += chisqwy;
  chisq += prior;

//  std::cout<<"returning chisq = "<<chisq<<" of which prior contributes "<<prior<<std::endl;

//  std::cout<<"xpos chisq \t ypos chisq \t xwid chisq \t ywid chisq \t prior"<<std::endl;
//  std//::cout<<std::setprecision(4)<<chisqx<<"\t"<<chisqy<<"\t"<<chisqwx<<"\t"<<chisqwy<<"\t"<<prior<<std::endl;
  if(std::isnan(chisq)) return 123456789.0;
  return chisq;
}

