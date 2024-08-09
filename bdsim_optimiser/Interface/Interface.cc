#include "Interface.h"
#include "Optics.h" 

std::string reBdsimOptics(std::string inputFileName){

  // emittance on the fly
  bool emittanceOnFly = false;

  std::string outputFileName = RBDS::DefaultOutputName(inputFileName, "_optics");

  DataLoader* dl = nullptr;
  try
    {dl = new DataLoader(inputFileName, false, true);}
  catch (const RBDSException& error)
    {std::cerr << error.what() << std::endl; throw;}
  catch (const std::exception& error)
    {std::cerr << error.what() << std::endl; throw;}

  // beam required to get the mass of the primary particle in EventAnalysis
  Beam*   beam     = dl->GetBeam();
  TChain* beamTree = dl->GetBeamTree();
  BDSOutputROOTEventBeam* outputBeam = beam->beam;
  beamTree->GetEntry(0);
  const std::string& particleName = outputBeam->particle;
  
  TChain* modelTree = dl->GetModelTree();
  if (modelTree->GetEntries() == 0)
    {
      std::cout << "Warning: data file written without Model tree that is required to know the sampler names" << std::endl;
      std::cout << "         only the primary sampler will be analysed if available" << std::endl;
    }

  EventAnalysis* evtAnalysis;
  try
    {
      evtAnalysis = new EventAnalysis(dl->GetEvent(), dl->GetEventTree(),
                                      false, true, false, true, -1, emittanceOnFly, 0, -1, particleName);
      evtAnalysis->Execute();
    }
  catch (const RBDSException& error)
    {std::cerr << error.what() << std::endl; throw;}

  TFile* outputFile = new TFile(outputFileName.c_str(), "RECREATE");

  // add header for file type and version details
  outputFile->cd();
  BDSOutputROOTEventHeader* headerOut = new BDSOutputROOTEventHeader();
  headerOut->Fill(dl->GetFileNames()); // updates time stamp
  headerOut->SetFileType("REBDSIM");
  TTree* headerTree = new TTree("Header", "REBDSIM Header");
  headerTree->Branch("Header.", "BDSOutputROOTEventHeader", headerOut);
  headerTree->Fill();
  headerTree->Write("", TObject::kOverwrite);

  // write merged histograms and optics
  evtAnalysis->Write(outputFile);

  // Don't clone the model tree if only primaries are generated - model not created in BDSIM
  Options* options = dl->GetOptions();
  TChain*  optionsTree = dl->GetOptionsTree();
  BDSOutputROOTEventOptions* ob = options->options;
  optionsTree->GetEntry(0);
  if (!ob->generatePrimariesOnly)
    {
      // clone model tree for nice built in optics plotting
      auto newTree = modelTree->CloneTree();
      newTree->Write("", TObject::kOverwrite);
    }
  
  outputFile->Close();
  delete outputFile;
  std::cout << "Result written to: " << outputFileName << std::endl;
  delete dl;
  delete evtAnalysis;

  return outputFileName;

}

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
    throw;
  }
  internalPars.resize(npars);
  nPars = npars;
}

void Interface::ParamScan(int param, int npoints, double xlo, double xup){
  double step = (xup-xlo)/(double)(npoints-1);
  double store = internalPars[param];
  internalPars[param] = xlo;
  const double tmp[1] = {0.99};
  for(int i=0; i<1; i++){
    std::cout << "Param " << param << "=" << internalPars[param]<< " chisq = " << fcn(internalPars);
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
  fcn(tmp);
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
  std::cout<<"Inside CalcChisq"<<std::endl;
  char *argv[6];
  char execname[256] = "bdsim";
  char path[256] = "--file=../gmad/optimised.gmad";
  char batch[256] = "--batch";
  char ngen[256] = "--ngenerate=200";
  char outfile[256] = "--outfile=/home/bdsim_output";
  char seed[256] = "--seed=1989";
  argv[0] = execname;
  argv[1] = path;
  argv[2] = batch;
  argv[3] = ngen;
  argv[4] = outfile;
  argv[5] = seed;

  std::system("bdsim --file=../gmad/optimised.gmad --batch --ngenerate=200 --outfile=/home/bdsim_output --seed=1989 > /dev/null");
  std::system("rebdsimOptics /home/bdsim_output.root /home/bdsim_output_optics.root > /dev/null");
//  BDSIM *bds = new BDSIM();
//  bds->Initialise(6, argv);
//  std::cout<<"Called new BDSIM"<<std::endl;
//  if (!bds->Initialised())
//    {
//      if (bds->InitialisationResult() == 1)
//        {std::cout << "Intialisation failed" << std::endl; return 1;}
//    }
//  else
//    {bds->BeamOn();}
//  delete bds;

//  std::string opticsFilename = reBdsimOptics("/home/bdsim_output.root");

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

