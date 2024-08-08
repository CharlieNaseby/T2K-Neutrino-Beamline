#include "Interface.h"


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


Interface::Interface(double *data){
    dat.resize(2);
    dat[0] = data[0];
    dat[1] = data[1];
};
Interface::~Interface(){};

double Interface::calc_chisq(const double *pars){

      std::cout<<"inside calc chisq"<<std::endl;
      char *argv[6];
      char path[256] = "--file=/home/extraction_to_arc.gmad";
      char batch[256] = "--batch";
      char ngen[256] = "--ngenerate=200";
      char outfile[256] = "--outfile=/home/bdsim_output";
      char seed[256] = "--seed=1989";
      argv[1] = path;
      argv[2] = batch;
      argv[3] = ngen;
      argv[4] = outfile;
      argv[5] = seed;
      BDSIM *bds = new BDSIM(6, argv);
      if (!bds->Initialised())
        {
          if (bds->InitialisationResult() == 1)
            {std::cout << "Intialisation failed" << std::endl; return 1;}
        }
      else
        {bds->BeamOn();}
      delete bds;

    reBdsimOptics("/home/bdsim_output.root");



    double chisq=0;
    for(int i=0; i<2; i++) chisq += (dat[i]-pars[i])*(dat[i]-pars[i])/pars[i];
    return chisq;
}

