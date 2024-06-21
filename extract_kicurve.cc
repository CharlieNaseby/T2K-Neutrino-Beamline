#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TGraph.h"
#include "TKey.h"


//this executable converts the file kicurve_sad.root available in the root dir of SAD
// to a pandas format usable by the BDSIM config creation python tool

int main(){

    TFile *inf = new TFile("kicurve_sad.root", "READ");

    TIter next(inf->GetListOfKeys());
    TKey *key;
    std::fstream out("kicurve.csv", std::ios::out);
    out << "element, current, kval\n";
    while ((key = (TKey*)next())) {
        TGraph *gr = (TGraph*)key->ReadObj();
        for(int i=0; i<gr->GetN(); i++){
            double x, y;
            gr->GetPoint(i, x, y);
            const char *name = key->GetName();
            out << (name+2) << ", " << x << ", " <<y <<"\n";
        }

    }
    out.close();
}