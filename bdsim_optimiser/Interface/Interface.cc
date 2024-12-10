#include "Interface.h"
#include "Optics.h" 

//constructor to set data and the gmad file we'll use as a base to our fit

Interface::Interface(std::string dataFile, std::string baseBeamlineFile, int npars, int nmagnetpars, int nbeampars){

  nPars = npars;
  nMagnetPars = nmagnetpars;
  nBeamPars = nbeampars;

  internalPars.resize(nPars); //point about which LLH scans will be conducted
  nominalPars.resize(nPars); //the values of the parameters expected based on true magnet current
  magNames.resize(nMagnetPars);
  beamNames.resize(nBeamPars);
  preFit.resize(nPars); //the values of the parameters expected with scalings or fudge factors applied


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

  double ssemxMean[NSSEM], ssemyMean[NSSEM], ssemwxMean[NSSEM], ssemwyMean[NSSEM];
  for(int i=0; i<NSSEM; i++){
    ssemxMean[i] = 0;
    ssemyMean[i] = 0;
    ssemwxMean[i] = 0;
    ssemwyMean[i] = 0;
  }
  int nShots=0;


  for(int i=0; i<ssemData->GetEntries(); i++){
    ssemData->GetEntry(i);
    if(ssemax[0] < 100) continue;  //if there was no beam skip
    for(int j=0; j<12; j++){ //CERN for some reason the later magnets do change from shot-to-shot
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
    std::cout<<magset[9]<<std::endl;
    nShots++;
loopend:;
  }
  for(int i=0; i<NSSEM; i++){
    ssemxMean[i] /= (double)nShots;
    ssemyMean[i] /= (double)nShots;
    ssemwxMean[i] /= (double)nShots;
    ssemwyMean[i] /= (double)nShots;
  }

  dat.resize(NSSEM);
  for(int i=0; i<dat.size(); i++){
    dat[i][0] = ssemxMean[i];
    dat[i][1] = ssemyMean[i];
    dat[i][2] = ssemwxMean[i];
    dat[i][3] = ssemwyMean[i];
  }

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
};

Interface::~Interface(){};

void Interface::SetInitialValues(bool usePrevBestFit, bool useFieldMaps, bool useFudgeFactor, double *pars){


  std::map<std::string, double> beamPars;
  beamPars["X0"]=-0.5;
  beamPars["Xp0"]=3.5e-5;
  beamPars["emitx"]=0.0610768;
  beamPars["betx"]=37.098;
  beamPars["alfx"]=-2.4187;
  beamPars["dispx"]=0.42373,
  beamPars["dispxp"]=0.07196,
  beamPars["Y0"]=-0.2;
  beamPars["Yp0"]=7.8e-5;
  beamPars["emity"]=0.05976;
  beamPars["bety"]=5.45;
  beamPars["alfy"]=0.178;
  beamPars["dispy"]=0.;
  beamPars["dispyp"]=0.0;





  std::map<std::string, double> kScaling;
  //derived from the KIcurve sad file, scaling of K1 or K0 value in units of 100A
  kScaling["BPV1"] = -0.08760146181454302;
  kScaling["BPH2"] = -0.06372506685790527;
  kScaling["QPQ1"] = 0.00953465;
  kScaling["QPQ2"] = -0.012669113;
  kScaling["BPD1"] = -0.0926275917278186;
  kScaling["BPD2"] = -0.08940591194129406;
  kScaling["QPQ3"] = 0.01593171666666667;
  kScaling["BPV2"] = -0.12238521642828552;
  kScaling["QPQ4"] = -0.01539898;
  kScaling["BPH3"] = -0.11826754864042695;
  kScaling["QPQ5"] = 0.01768722;


  //values based on the fieldmap fields at 100A
  kScaling["QPQ1"] = 0.011798078333463476;
  kScaling["QPQ2"] = -0.011492647647480913;
  kScaling["QPQ4"] = -0.01471238222852117;
  kScaling["QPQ5"] = 0.019445574477550798;

  kScaling["BPV1"] = -0.09456593004136546;
  kScaling["BPD1"] = -0.09289046446110355;
  kScaling["BPD2"] = -0.09289046446110355;
  kScaling["BPV2"] = -0.12654858254158613;
  kScaling["BPH3"] = -0.12432151082458207;

  nominalPars[0] = magCurrent[0]/100.;
  nominalPars[1] = magCurrent[1]/100.;
  nominalPars[2] = magCurrent[2]/100.;
  nominalPars[3] = magCurrent[4]/100.;
  nominalPars[4] = magCurrent[5]/100.;
  nominalPars[5] = magCurrent[6]/100.;
  nominalPars[6] = magCurrent[7]/100.;
  nominalPars[7] = magCurrent[8]/100.;
  nominalPars[8] = magCurrent[9]/100.;
  nominalPars[9] = magCurrent[10]/100.;
  nominalPars[10] = magCurrent[11]/100.;

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
    nominalPars[magMap[mag]] *= kScaling[mag];//nominal is always the expected value for the parameters based on the currents
    preFit[magMap[mag]] = nominalPars[magMap[mag]];
    if(useFudgeFactor) preFit[magMap[mag]] *= fudgeFactor[magMap[mag]];
  }

  if(usePrevBestFit){
    std::ifstream infile;
    infile.open("bestFitHardEdge.txt");
    std::string line;

    while(std::getline(infile, line)){ //TODO broken with addition of beam parameters
      size_t strpos = line.find(" =");
      std::string mag = line.substr(0, strpos);
      std::string val = line.substr(strpos+3);
      pars[magMap[mag]] = std::stof(val);
      std::cout<<"setting "<<magMap[mag]<<" magname "<<mag<<" to "<<std::stof(val)<<std::endl;
    }
  }
  else if(useFieldMaps){
    for(int i=0; i<nMagnetPars; i++) fudgeFactor[i] = 1.0;
  
//    fudgeFactor[1] = 1.402333;
//    fudgeFactor[4] = 0.851299;
//    fudgeFactor[5] = 0.830461;
//    fudgeFactor[7] = 1.023666;

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
      nominalPars[i] = pars[i];
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
  }

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

const double* Interface::GetNominalPars(){
  return &nominalPars[0];
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

double Interface::fcn(const double *pars){
  SetInternalPars(pars);
  GenerateInputFile(pars);
  return CalcChisq(pars);
}

double Interface::CalcPrior(const double *pars){
  double chisq = 0;
  for(int i=0; i<nPars; i++){
    chisq += (pars[i]-nominalPars[i])*(pars[i]-nominalPars[i])/(1e-10+(0.1*nominalPars[i])*(0.1*nominalPars[i]));
  }
  return chisq;
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
  std::cout<<"parsing input file "<<baseBeamlineFile<<std::endl;
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
//    //does it have the string bScaling in it?
//    size_t bspos = line.find("bScaling=");
//    size_t k0pos = line.find("B=");
//    size_t k1pos = line.find("k1=");
//    if(bspos == std::string::npos && k0pos==std::string::npos && k1pos==std::string::npos) beamline[beamline.size()-1].append(line);
//    else{
//      std::string label;
//      size_t strpos = 0;
//      if(bspos!=std::string::npos){
//        strpos = bspos;
//        label = "bScaling=";
//      }
//      if(k0pos!=std::string::npos){
//        strpos = k0pos;
//        label = "B=";
//      }
//      if(k1pos!=std::string::npos){
//        strpos = k1pos;
//        label = "k1=";
//      }
//      std::string firstPart = line.substr(0, strpos);
//      firstPart.append(label);
//      beamline[beamline.size()-1].append(firstPart);
//      std::string secondPart = line.substr(strpos);
//      size_t comma = secondPart.find_first_of(",");
//      std::string newSecondPart = secondPart.substr(comma);
//      
//      beamline.push_back(newSecondPart);
//    }
//  }
  if(nPars != beamline.size()-1){
    std::cerr<<"Number of parameters set is not equal to the number in the gmad file supplied "<<__FILE__<<":"<<__LINE__<<std::endl;
    std::cerr<<"Expected "<<nPars<<" but file has "<<beamline.size()-1<<std::endl;
    throw  ;
  }
}

void Interface::CalcBeamPars(){
    TFile *fitFile = new TFile("/home/bdsim_output.root");
    TTree *samplerData = static_cast<TTree*>(fitFile->Get("Event"));
    BDSOutputROOTEventSampler<float> *SSEM1Data = nullptr; //new BDSOutputROOTEventSampler<float>("SSEM1");
    samplerData->SetBranchAddress("SSEM1.", &SSEM1Data);

    for(int i=0; i<samplerData->GetEntries(); i++){
        samplerData->GetEntry(i);
        std::cout<<SSEM1Data->samplerName<<std::endl;
    }
}


double Interface::CalcChisq(const double *pars){

  std::system("bdsim --file=../gmad/optimised.gmad --batch --ngenerate=100 --outfile=/home/bdsim_output --seed=1989 > /dev/null");
  std::system("rebdsimOptics /home/bdsim_output.root /home/bdsim_output_optics.root > /dev/null");

//  CalcBeamPars();

  Optics beamOptics("/home/bdsim_output_optics.root");


  double chisq=0;
  double chisqx = 0;
  double chisqy = 0;
  double chisqwx = 0;
  double chisqwy = 0;

  for(int i=0; i<dat.size(); i++){
    beamOptics.fChain->GetEntry(i+1);
    std::array<double, 4> simulation = {1000.*beamOptics.Mean_x, 1000.*beamOptics.Mean_y, 2000.*beamOptics.Sigma_x, 2000.*beamOptics.Sigma_y};
    std::cout<<"SSEM"<<i+1<<" beam sim postion = \t"<<simulation[0]<<", \t"<<simulation[1]<<" data \t"<<dat[i][0]<<", \t"<<dat[i][1]<<std::endl;
    std::cout<<"SSEM"<<i+1<<" beam sim width   = \t"<<simulation[2]<<", \t"<<simulation[3]<<" data \t"<<dat[i][2]<<", \t"<<dat[i][3]<<std::endl;
    chisqx += (dat[i][0]-simulation[0])*(dat[i][0]-simulation[0])/(0.2*0.2);  //CERN 0.2mm uncert on ssem position x
    chisqy += (dat[i][1]-simulation[1])*(dat[i][1]-simulation[1])/(0.2*0.2);  //position y
    chisqwx += (dat[i][2]-simulation[2])*(dat[i][2]-simulation[2])/(0.2*0.2);  //CERN width with 0.2mm precision x
    chisqwy += (dat[i][3]-simulation[3])*(dat[i][3]-simulation[3])/(0.2*0.2);  //width y
//    std::cout<<"SSEM "<<i+1<<" position cumulative chisq contribution "<<chisq<<std::endl;
  }
//  chisq+=CalcPrior(pars);
  if(fitMode & 0x01) chisq += chisqx;
  if(fitMode & 0x02) chisq += chisqy;
  if(fitMode & 0x04) chisq += chisqwx;
  if(fitMode & 0x08) chisq += chisqwy;

  std::cout<<"returning chisq = "<<chisq<<std::endl;

  std::cout<<"xpos chisq \t ypos chisq \t xwid chisq \t ywid chisq"<<std::endl;
  std::cout<<std::setprecision(4)<<chisqx<<"\t"<<chisqy<<"\t"<<chisqwx<<"\t"<<chisqwy<<std::endl;
  return chisq;
}

