#pragma once

class Parameter{
private:
    double fitValue, physicalValue;
    std::string permenantState; //free, fixed
    std::string state; //free, fixed, or temporarily fixed for a fit
public:
    std::string name;
    std::string type; //hBend, vBend, quad, beam
    std::string section; //preparation, arc, finalFocus
    int currentIdx; //index of the magnet in the input file magnet current array
    double prefitValue, nominalValue, priorError; //all except fitValue are in real, physical units
    double boundLo = -1e30;
    double boundHi = 1e30;


    Parameter(std::string nameIn, std::string stateIn, std::string parType, std::string sectionIn="NA", int currentIdxIn=-1):
        permenantState(stateIn), name(nameIn), type(parType), section(sectionIn),  currentIdx(currentIdxIn)
        {;}



    std::string getState(){return state;}
    void setState(std::string s){
        if(permenantState == "fixed" && s == "free") std::cout << "Caution attempting to set a permanently fixed parameter to free"<<std::endl;
        else{
          state = s;
        }
    }


    void setPhysicalValue(double physval){physicalValue = physval; fitValue = physval/prefitValue;}
    void setFitValue(double fitval){fitValue = fitval; physicalValue = fitval*prefitValue;}
    double getPhysicalValue(){return physicalValue;}
    double getFitValue(){return fitValue;}
};
