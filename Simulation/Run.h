#ifndef RUN_H
#define RUN_H

#include "Event.h"
#include "yaml-cpp/yaml.h"
#include "TH1D.h"
#include <vector>
#include <string>

class Run: public TObject
{
/*
 *  Class for executing a Run
 *  -------------------------
 *  Parameters:
 *  cfgFileName: string containing yaml file for configurations
 * 
 */
public:
    Run() = default;
    Run(TString cfgFileName);

private:
    void RunConstMult(), RunUniformMult(), RunCustomMult(), CreateDetectors();
    YAML::Node fConfigFile;
    string fOutFileName;
    TTree fTreeGen;
    TTree fTreeRec;

    //Multiplicity info
    unsigned fNEvents,fConstMult;
    string fMultType,fMultFile,fMultHisto;
    vector<unsigned> fMultRange;

    //Vertex dispersion info
    double fSigmaZ,fSigmaX,fSigmaY;

    //Detectors info
    vector<MaterialBudget*> fDetectors;
    vector<bool> fIsDetector;
    vector<double> fRadii, fThickness, fLength;
    vector<string> fMaterial;

    //Reconstruction & simulation options
    bool fDoMultScattering, fDoSmearing, fDoNoise;
    int fMeanNoise;

    bool fVerbose;

    ClassDef(Run, 1)
};


#endif