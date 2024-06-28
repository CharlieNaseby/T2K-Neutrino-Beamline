//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr 22 15:45:05 2024 by ROOT version 6.24/00
// from TTree Optics/Optics
// found on file: test.root
//////////////////////////////////////////////////////////

#ifndef Optics_h
#define Optics_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Optics {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        Emitt_x;
   Double_t        Emitt_y;
   Double_t        Alpha_x;
   Double_t        Alpha_y;
   Double_t        Beta_x;
   Double_t        Beta_y;
   Double_t        Gamma_x;
   Double_t        Gamma_y;
   Double_t        Disp_x;
   Double_t        Disp_y;
   Double_t        Disp_xp;
   Double_t        Disp_yp;
   Double_t        Mean_x;
   Double_t        Mean_y;
   Double_t        Mean_xp;
   Double_t        Mean_yp;
   Double_t        Sigma_x;
   Double_t        Sigma_y;
   Double_t        Sigma_xp;
   Double_t        Sigma_yp;
   Double_t        S;
   Double_t        Npart;
   Double_t        Sigma_Emitt_x;
   Double_t        Sigma_Emitt_y;
   Double_t        Sigma_Alpha_x;
   Double_t        Sigma_Alpha_y;
   Double_t        Sigma_Beta_x;
   Double_t        Sigma_Beta_y;
   Double_t        Sigma_Gamma_x;
   Double_t        Sigma_Gamma_y;
   Double_t        Sigma_Disp_x;
   Double_t        Sigma_Disp_y;
   Double_t        Sigma_Disp_xp;
   Double_t        Sigma_Disp_yp;
   Double_t        Sigma_Mean_x;
   Double_t        Sigma_Mean_y;
   Double_t        Sigma_Mean_xp;
   Double_t        Sigma_Mean_yp;
   Double_t        Sigma_Sigma_x;
   Double_t        Sigma_Sigma_y;
   Double_t        Sigma_Sigma_xp;
   Double_t        Sigma_Sigma_yp;
   Double_t        Mean_E;
   Double_t        Mean_t;
   Double_t        Sigma_E;
   Double_t        Sigma_t;
   Double_t        Sigma_Mean_E;
   Double_t        Sigma_Mean_t;
   Double_t        Sigma_Sigma_E;
   Double_t        Sigma_Sigma_t;
   Double_t        xyCorrelationCoefficent;

   // List of branches
   TBranch        *b_Emitt_x;   //!
   TBranch        *b_Emitt_y;   //!
   TBranch        *b_Alpha_x;   //!
   TBranch        *b_Alpha_y;   //!
   TBranch        *b_Beta_x;   //!
   TBranch        *b_Beta_y;   //!
   TBranch        *b_Gamma_x;   //!
   TBranch        *b_Gamma_y;   //!
   TBranch        *b_Disp_x;   //!
   TBranch        *b_Disp_y;   //!
   TBranch        *b_Disp_xp;   //!
   TBranch        *b_Disp_yp;   //!
   TBranch        *b_Mean_x;   //!
   TBranch        *b_Mean_y;   //!
   TBranch        *b_Mean_xp;   //!
   TBranch        *b_Mean_yp;   //!
   TBranch        *b_Sigma_x;   //!
   TBranch        *b_Sigma_y;   //!
   TBranch        *b_Sigma_xp;   //!
   TBranch        *b_Sigma_yp;   //!
   TBranch        *b_S;   //!
   TBranch        *b_Npart;   //!
   TBranch        *b_Sigma_Emitt_x;   //!
   TBranch        *b_Sigma_Emitt_y;   //!
   TBranch        *b_Sigma_Alpha_x;   //!
   TBranch        *b_Sigma_Alpha_y;   //!
   TBranch        *b_Sigma_Beta_x;   //!
   TBranch        *b_Sigma_Beta_y;   //!
   TBranch        *b_Sigma_Gamma_x;   //!
   TBranch        *b_Sigma_Gamma_y;   //!
   TBranch        *b_Sigma_Disp_x;   //!
   TBranch        *b_Sigma_Disp_y;   //!
   TBranch        *b_Sigma_Disp_xp;   //!
   TBranch        *b_Sigma_Disp_yp;   //!
   TBranch        *b_Sigma_Mean_x;   //!
   TBranch        *b_Sigma_Mean_y;   //!
   TBranch        *b_Sigma_Mean_xp;   //!
   TBranch        *b_Sigma_Mean_yp;   //!
   TBranch        *b_Sigma_Sigma_x;   //!
   TBranch        *b_Sigma_Sigma_y;   //!
   TBranch        *b_Sigma_Sigma_xp;   //!
   TBranch        *b_Sigma_Sigma_yp;   //!
   TBranch        *b_Mean_E;   //!
   TBranch        *b_Mean_t;   //!
   TBranch        *b_Sigma_E;   //!
   TBranch        *b_Sigma_t;   //!
   TBranch        *b_Sigma_Mean_E;   //!
   TBranch        *b_Sigma_Mean_t;   //!
   TBranch        *b_Sigma_Sigma_E;   //!
   TBranch        *b_Sigma_Sigma_t;   //!
   TBranch        *b_xyCorrelationCoefficent;   //!

   Optics(char *infilename);
   virtual ~Optics();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     dump_as_csv(std::string filename);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Optics_cxx
Optics::Optics(char *infilename) : fChain(0) 
{
   TFile *infile = new TFile(infilename);
   TTree *tree = (TTree*)infile->Get("Optics");


// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
//   if (tree == 0) {
//      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("test.root");
//      if (!f || !f->IsOpen()) {
//         f = new TFile("test.root");
//      }
//      f->GetObject("Optics",tree);
//
//   }
   Init(tree);
}

Optics::~Optics()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Optics::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Optics::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Optics::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Emitt_x", &Emitt_x, &b_Emitt_x);
   fChain->SetBranchAddress("Emitt_y", &Emitt_y, &b_Emitt_y);
   fChain->SetBranchAddress("Alpha_x", &Alpha_x, &b_Alpha_x);
   fChain->SetBranchAddress("Alpha_y", &Alpha_y, &b_Alpha_y);
   fChain->SetBranchAddress("Beta_x", &Beta_x, &b_Beta_x);
   fChain->SetBranchAddress("Beta_y", &Beta_y, &b_Beta_y);
   fChain->SetBranchAddress("Gamma_x", &Gamma_x, &b_Gamma_x);
   fChain->SetBranchAddress("Gamma_y", &Gamma_y, &b_Gamma_y);
   fChain->SetBranchAddress("Disp_x", &Disp_x, &b_Disp_x);
   fChain->SetBranchAddress("Disp_y", &Disp_y, &b_Disp_y);
   fChain->SetBranchAddress("Disp_xp", &Disp_xp, &b_Disp_xp);
   fChain->SetBranchAddress("Disp_yp", &Disp_yp, &b_Disp_yp);
   fChain->SetBranchAddress("Mean_x", &Mean_x, &b_Mean_x);
   fChain->SetBranchAddress("Mean_y", &Mean_y, &b_Mean_y);
   fChain->SetBranchAddress("Mean_xp", &Mean_xp, &b_Mean_xp);
   fChain->SetBranchAddress("Mean_yp", &Mean_yp, &b_Mean_yp);
   fChain->SetBranchAddress("Sigma_x", &Sigma_x, &b_Sigma_x);
   fChain->SetBranchAddress("Sigma_y", &Sigma_y, &b_Sigma_y);
   fChain->SetBranchAddress("Sigma_xp", &Sigma_xp, &b_Sigma_xp);
   fChain->SetBranchAddress("Sigma_yp", &Sigma_yp, &b_Sigma_yp);
   fChain->SetBranchAddress("S", &S, &b_S);
   fChain->SetBranchAddress("Npart", &Npart, &b_Npart);
   fChain->SetBranchAddress("Sigma_Emitt_x", &Sigma_Emitt_x, &b_Sigma_Emitt_x);
   fChain->SetBranchAddress("Sigma_Emitt_y", &Sigma_Emitt_y, &b_Sigma_Emitt_y);
   fChain->SetBranchAddress("Sigma_Alpha_x", &Sigma_Alpha_x, &b_Sigma_Alpha_x);
   fChain->SetBranchAddress("Sigma_Alpha_y", &Sigma_Alpha_y, &b_Sigma_Alpha_y);
   fChain->SetBranchAddress("Sigma_Beta_x", &Sigma_Beta_x, &b_Sigma_Beta_x);
   fChain->SetBranchAddress("Sigma_Beta_y", &Sigma_Beta_y, &b_Sigma_Beta_y);
   fChain->SetBranchAddress("Sigma_Gamma_x", &Sigma_Gamma_x, &b_Sigma_Gamma_x);
   fChain->SetBranchAddress("Sigma_Gamma_y", &Sigma_Gamma_y, &b_Sigma_Gamma_y);
   fChain->SetBranchAddress("Sigma_Disp_x", &Sigma_Disp_x, &b_Sigma_Disp_x);
   fChain->SetBranchAddress("Sigma_Disp_y", &Sigma_Disp_y, &b_Sigma_Disp_y);
   fChain->SetBranchAddress("Sigma_Disp_xp", &Sigma_Disp_xp, &b_Sigma_Disp_xp);
   fChain->SetBranchAddress("Sigma_Disp_yp", &Sigma_Disp_yp, &b_Sigma_Disp_yp);
   fChain->SetBranchAddress("Sigma_Mean_x", &Sigma_Mean_x, &b_Sigma_Mean_x);
   fChain->SetBranchAddress("Sigma_Mean_y", &Sigma_Mean_y, &b_Sigma_Mean_y);
   fChain->SetBranchAddress("Sigma_Mean_xp", &Sigma_Mean_xp, &b_Sigma_Mean_xp);
   fChain->SetBranchAddress("Sigma_Mean_yp", &Sigma_Mean_yp, &b_Sigma_Mean_yp);
   fChain->SetBranchAddress("Sigma_Sigma_x", &Sigma_Sigma_x, &b_Sigma_Sigma_x);
   fChain->SetBranchAddress("Sigma_Sigma_y", &Sigma_Sigma_y, &b_Sigma_Sigma_y);
   fChain->SetBranchAddress("Sigma_Sigma_xp", &Sigma_Sigma_xp, &b_Sigma_Sigma_xp);
   fChain->SetBranchAddress("Sigma_Sigma_yp", &Sigma_Sigma_yp, &b_Sigma_Sigma_yp);
   fChain->SetBranchAddress("Mean_E", &Mean_E, &b_Mean_E);
   fChain->SetBranchAddress("Mean_t", &Mean_t, &b_Mean_t);
   fChain->SetBranchAddress("Sigma_E", &Sigma_E, &b_Sigma_E);
   fChain->SetBranchAddress("Sigma_t", &Sigma_t, &b_Sigma_t);
   fChain->SetBranchAddress("Sigma_Mean_E", &Sigma_Mean_E, &b_Sigma_Mean_E);
   fChain->SetBranchAddress("Sigma_Mean_t", &Sigma_Mean_t, &b_Sigma_Mean_t);
   fChain->SetBranchAddress("Sigma_Sigma_E", &Sigma_Sigma_E, &b_Sigma_Sigma_E);
   fChain->SetBranchAddress("Sigma_Sigma_t", &Sigma_Sigma_t, &b_Sigma_Sigma_t);
   fChain->SetBranchAddress("xyCorrelationCoefficent", &xyCorrelationCoefficent, &b_xyCorrelationCoefficent);
   Notify();
}

Bool_t Optics::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Optics::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Optics::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Optics_cxx
