#define Optics_cxx
#include "Optics.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <string>

void Optics::dump_as_csv(std::string filename)
{
   Long64_t nentries = fChain->GetEntriesFast();

   std::cout<<"saving..."<<std::endl;

   std::fstream outfile(filename, std::ios::out);
   outfile<<"entry,";
   outfile<<"Emitt_x,";
   outfile<<"Emitt_y,";
   outfile<<"Alpha_x,";
   outfile<<"Alpha_y,";
   outfile<<"Beta_x,";
   outfile<<"Beta_y,";
   outfile<<"Gamma_x,";
   outfile<<"Gamma_y,";
   outfile<<"Disp_x,";
   outfile<<"Disp_y,";
   outfile<<"Disp_xp,";
   outfile<<"Disp_yp,";
   outfile<<"Mean_x,";
   outfile<<"Mean_y,";
   outfile<<"Mean_xp,";
   outfile<<"Mean_yp,";
   outfile<<"Sigma_x,";
   outfile<<"Sigma_y,";
   outfile<<"Sigma_xp,";
   outfile<<"Sigma_yp,";
   outfile<<"S,";
   outfile<<"Npart,";
   outfile<<"Sigma_Emitt_x,";
   outfile<<"Sigma_Emitt_y,";
   outfile<<"Sigma_Alpha_x,";
   outfile<<"Sigma_Alpha_y,";
   outfile<<"Sigma_Beta_x,";
   outfile<<"Sigma_Beta_y,";
   outfile<<"Sigma_Gamma_x,";
   outfile<<"Sigma_Gamma_y,";
   outfile<<"Sigma_Disp_x,";
   outfile<<"Sigma_Disp_y,";
   outfile<<"Sigma_Disp_xp,";
   outfile<<"Sigma_Disp_yp,";
   outfile<<"Sigma_Mean_x,";
   outfile<<"Sigma_Mean_y,";
   outfile<<"Sigma_Mean_xp,";
   outfile<<"Sigma_Mean_yp,";
   outfile<<"Sigma_Sigma_x,";
   outfile<<"Sigma_Sigma_y,";
   outfile<<"Sigma_Sigma_xp,";
   outfile<<"Sigma_Sigma_yp,";
   outfile<<"Mean_E,";
   outfile<<"Mean_t,";
   outfile<<"Sigma_E,";
   outfile<<"Sigma_t,";
   outfile<<"Sigma_Mean_E,";
   outfile<<"Sigma_Mean_t,";
   outfile<<"Sigma_Sigma_E,";
   outfile<<"Sigma_Sigma_t,";
   outfile<<"xyCorrelationCoefficent"<<"\n";


   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      outfile<<jentry<<",";
      outfile<<Emitt_x<<",";
      outfile<<Emitt_y<<",";
      outfile<<Alpha_x<<",";
      outfile<<Alpha_y<<",";
      outfile<<Beta_x<<",";
      outfile<<Beta_y<<",";
      outfile<<Gamma_x<<",";
      outfile<<Gamma_y<<",";
      outfile<<Disp_x<<",";
      outfile<<Disp_y<<",";
      outfile<<Disp_xp<<",";
      outfile<<Disp_yp<<",";
      outfile<<Mean_x<<",";
      outfile<<Mean_y<<",";
      outfile<<Mean_xp<<",";
      outfile<<Mean_yp<<",";
      outfile<<Sigma_x<<",";
      outfile<<Sigma_y<<",";
      outfile<<Sigma_xp<<",";
      outfile<<Sigma_yp<<",";
      outfile<<S<<",";
      outfile<<Npart<<",";
      outfile<<Sigma_Emitt_x<<",";
      outfile<<Sigma_Emitt_y<<",";
      outfile<<Sigma_Alpha_x<<",";
      outfile<<Sigma_Alpha_y<<",";
      outfile<<Sigma_Beta_x<<",";
      outfile<<Sigma_Beta_y<<",";
      outfile<<Sigma_Gamma_x<<",";
      outfile<<Sigma_Gamma_y<<",";
      outfile<<Sigma_Disp_x<<",";
      outfile<<Sigma_Disp_y<<",";
      outfile<<Sigma_Disp_xp<<",";
      outfile<<Sigma_Disp_yp<<",";
      outfile<<Sigma_Mean_x<<",";
      outfile<<Sigma_Mean_y<<",";
      outfile<<Sigma_Mean_xp<<",";
      outfile<<Sigma_Mean_yp<<",";
      outfile<<Sigma_Sigma_x<<",";
      outfile<<Sigma_Sigma_y<<",";
      outfile<<Sigma_Sigma_xp<<",";
      outfile<<Sigma_Sigma_yp<<",";
      outfile<<Mean_E<<",";
      outfile<<Mean_t<<",";
      outfile<<Sigma_E<<",";
      outfile<<Sigma_t<<",";
      outfile<<Sigma_Mean_E<<",";
      outfile<<Sigma_Mean_t<<",";
      outfile<<Sigma_Sigma_E<<",";
      outfile<<Sigma_Sigma_t<<",";
      outfile<<xyCorrelationCoefficent<<"\n";
   }
   outfile.close();
}


void usage(){
    std::cout<<"usage: ./optics_to_csv input.root output.csv"<<std::endl;
}


int main(int argc, char **argv){
    if(argc!=3){
        usage();
        exit(1);
    }
    Optics csv_dumper(argv[1]);
    csv_dumper.dump_as_csv(std::string(argv[2]));
}



