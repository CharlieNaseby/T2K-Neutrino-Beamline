#include "Interface.h"
#include "Optics.h" 

//constructor to set data and the gmad file we'll use as a base to our fit

Interface::Interface(std::string dataFile, std::string baseBeamlineFile, int npars){

  nPars = npars;
  internalPars.resize(npars);
  nominalPars.resize(npars);
  magNames.resize(nPars);
  preFit.resize(nPars);


  TFile *ssemDataFile = new TFile(dataFile.c_str(), "READ");
  TTree *ssemData = dynamic_cast<TTree*> (ssemDataFile->Get("anabeam"));
  double ssemx[19], ssemy[19], ssemwx[19], ssemwy[19], ct[5], magset[31];
  ssemData->SetBranchAddress("magset", magset);
  ssemData->SetBranchAddress("ssemx", ssemx);
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
    if(ct[4] < 100) continue;
    for(int j=0; j<NSSEM; j++){
      ssemxMean[j]+=ssemx[j];
      ssemyMean[j]+=ssemy[j];
      ssemwxMean[j]+=ssemwx[j];
      ssemwyMean[j]+=ssemwy[j];
    }
    nShots++;    
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
  ParseInputFile(baseBeamlineFile);
};

Interface::~Interface(){};

void Interface::SetInitialValues(bool usePrevBestFit, bool useFieldMaps, double *pars){

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

  preFit[0] = magCurrent[0]/100.;
  preFit[1] = magCurrent[1]/100.;
  preFit[2] = magCurrent[2]/100.;
  preFit[3] = magCurrent[4]/100.;
  preFit[4] = magCurrent[5]/100.;
  preFit[5] = magCurrent[6]/100.;
  preFit[6] = magCurrent[7]/100.;
  preFit[7] = magCurrent[8]/100.;
  preFit[8] = magCurrent[9]/100.;
  preFit[9] = magCurrent[10]/100.;
  preFit[10] = magCurrent[11]/100.;

  std::vector<double> fudgeFactor(nPars);
  //magnet scaling factors based on a fit to hard edge run0910580 with prefit using fieldmap vals
  fudgeFactor[0] = 1;
  fudgeFactor[1] = 0.503659;
  fudgeFactor[2] = 0.878798;
  fudgeFactor[3] = 1.11968;
  fudgeFactor[4] = 1.08873;
  fudgeFactor[5] = 1.03649;
  fudgeFactor[6] = 1.4863;
  fudgeFactor[7] = 1.59548;
  fudgeFactor[8] = 1.28218;
  fudgeFactor[9] = 1;
  fudgeFactor[10] = 1;

  for(auto mag : magNames) preFit[magMap[mag]] *= kScaling[mag] * fudgeFactor[magMap[mag]]; //prefit is always the expected value for the parameters based on the currents


  if(usePrevBestFit){
    std::ifstream infile;
    infile.open("bestFitHardEdge.txt");
    std::string line;

    while(std::getline(infile, line)){
      size_t strpos = line.find(" =");
      std::string mag = line.substr(0, strpos);
      std::string val = line.substr(strpos+3);
      pars[magMap[mag]] = std::stof(val);
      std::cout<<"setting "<<magMap[mag]<<" magname "<<mag<<" to "<<std::stof(val)<<std::endl;
    }
  }
  else if(useFieldMaps){
    for(int i=0; i<nPars; i++) fudgeFactor[i] = 1.0;
  
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
    for(int i=0; i<nPars; i++) pars[i] *= fudgeFactor[i];

    //pars[0] = -0.05;  //CERN TODO forcefully set BPV1
    for(int i=0; i<nPars; i++) preFit[i] = pars[i];
  }
  else{
    for(int i=0; i<nPars; i++) pars[i] = preFit[i];
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

void Interface::SetNominalPars(const double *pars){
  for(int i=0; i<nPars; i++) nominalPars[i] = pars[i];
}

const double* Interface::GetNominalPars(){
  return &nominalPars[0];
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
  std::ifstream infile;
  infile.open(baseBeamlineFile);
  beamline.resize(1);
  std::string line;

  while(std::getline(infile, line)){
    //does it have the string bScaling in it?
    size_t bspos = line.find("bScaling=");
    size_t k0pos = line.find("B=");
    size_t k1pos = line.find("k1=");
    line.append("\n");
    if(bspos == std::string::npos && k0pos==std::string::npos && k1pos==std::string::npos) beamline[beamline.size()-1].append(line);
    else{
      std::string label;
      size_t strpos = 0;
      if(bspos!=std::string::npos){
        strpos = bspos;
        label = "bScaling=";
      }
      if(k0pos!=std::string::npos){
        strpos = k0pos;
        label = "B=";
      }
      if(k1pos!=std::string::npos){
        strpos = k1pos;
        label = "k1=";
      }
      std::string firstPart = line.substr(0, strpos);
      firstPart.append(label);
      beamline[beamline.size()-1].append(firstPart);
      std::string secondPart = line.substr(strpos);
      size_t comma = secondPart.find_first_of(",");
      std::string newSecondPart = secondPart.substr(comma);
      
      beamline.push_back(newSecondPart);
    }
  }
  if(nPars != beamline.size()-1){
    std::cerr<<"Number of parameters set is not equal to the number in the gmad file supplied FILE "<<__FILE__<<":"<<__LINE__<<std::endl;
    throw  ;
  }
}

double Interface::CalcChisq(const double *pars){

  std::system("bdsim --file=../gmad/optimised.gmad --batch --ngenerate=50 --outfile=/home/bdsim_output --seed=1989 > /dev/null");
  std::system("rebdsimOptics /home/bdsim_output.root /home/bdsim_output_optics.root > /dev/null");

  Optics beamOptics("/home/bdsim_output_optics.root");

  double chisq=0;
  for(int i=0; i<dat.size(); i++){
    beamOptics.fChain->GetEntry(i+1);
    std::array<double, 4> simulation = {1000.*beamOptics.Mean_x, 1000.*beamOptics.Mean_y, 2000.*beamOptics.Sigma_x, 2000.*beamOptics.Sigma_y};
    std::cout<<"SSEM"<<i+1<<" beam sim postion = "<<simulation[0]<<", "<<simulation[1]<<" data "<<dat[i][0]<<", "<<dat[i][1]<<std::endl;
    std::cout<<"SSEM"<<i+1<<" beam sim width   = "<<simulation[2]<<", "<<simulation[3]<<" data "<<dat[i][2]<<", "<<dat[i][3]<<std::endl;

    if(fitMode & 0x01) for(int j=0; j<2; j++) chisq += (dat[i][j]-simulation[j])*(dat[i][j]-simulation[j])/(0.2);  //CERN 0.2mm uncert on ssem position
    if(fitMode & 0x02) for(int j=2; j<4; j++) chisq += (dat[i][j]-simulation[j])*(dat[i][j]-simulation[j])/(0.2);  //CERN width with 0.2mm precision
   

//    std::cout<<"SSEM "<<i+1<<" position cumulative chisq contribution "<<chisq<<std::endl;
  }
//  chisq+=CalcPrior(pars);
  std::cout<<"returning chisq = "<<chisq<<std::endl;
  return chisq;
}

