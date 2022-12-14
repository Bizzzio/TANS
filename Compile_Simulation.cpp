#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"

void Compile_Simulation(TString myopt="fast")
{
    gSystem->AddIncludePath((string("-I ")+gSystem->GetWorkingDirectory()+"/Simulation").c_str());
    TString opt;
    if (myopt.Contains("force"))
        opt="kfg";          //k mantiene il .so, f forza la compilazione (come mettere il ++), g serve se vogliamo usare un debugger
    else
        opt="kg";
    if (myopt.Contains("clean"))
        gSystem->Exec("./clean.sh");
    gSystem->CompileMacro("Simulation/Particle.cpp",opt.Data());
    gSystem->CompileMacro("Simulation/MaterialBudget.cpp",opt.Data());
    gSystem->CompileMacro("Simulation/Detector.cpp",opt.Data());
    gSystem->CompileMacro("Simulation/Event.cpp",opt.Data());
    gSystem->CompileMacro("Simulation/Run.cpp",opt.Data());
    gSystem->CompileMacro("Simulation/Simulation.cpp",opt.Data());
    if (myopt.Contains("run"))
        gROOT->ProcessLine("Simulation()");
}