#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include"Run.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <vector>
#include <string>
 
class Reconstruction: public TObject
{
public:
    Reconstruction();
    ~Reconstruction(){delete fResiduals;}
 
    void VertexReco();
    void MinDca();
    void TrackletsReco();
    std::vector<double> GetTrackletParameters(int j);
    // void VertexReco(int index, double t);
    void MinGlobalDistance();
    MaterialBudget::fPoint FirstOctant(MaterialBudget::fPoint point);
    vector<double> Line(int index, double time);
    vector<double> TimesMinGlobalDistance(int index0);
 
private:
    std::vector<std::vector<MaterialBudget::fPoint>> fIntersections1;
    std::vector<std::vector<MaterialBudget::fPoint>> fIntersections2;
    std::vector<std::vector<Event::fVertMult>> fConfigs;
    std::vector<std::vector<MaterialBudget::fPoint>> fTracklets;
    std::vector<std::vector<double>> fVertexes;
    std::vector<double> fVertexesZ;
    TH1D* fResiduals = new TH1D("Residuals", "Residuals", 500,-0.5,0.5);

    void FillHistoMinDca(TH1D* histo, vector<MaterialBudget::fPoint>& tracklets, vector<double>& vertextemp);
    void FillHistoResiduals();

 
    ClassDef(Reconstruction, 1)
};
 
 
#endif
