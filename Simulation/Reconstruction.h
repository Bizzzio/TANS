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
    std::vector<double> GetTrackletParameters(MaterialBudget::fPoint point1, MaterialBudget::fPoint point2);
    void VertexRecoGeom();
    vector<double> Line(MaterialBudget::fPoint point1, MaterialBudget::fPoint point2, vector<double> velocity, double time);
    void FillHistoVertexGeom(vector<MaterialBudget::fPoint> acceptedtracks, vector<vector<double>> velocities, TH1D* graph);
    double GetTotalDistance(vector<MaterialBudget::fPoint> acceptedtracks, vector<vector<double>> velocities, double time, double mindistance);


 
private:
    std::vector<std::vector<MaterialBudget::fPoint>> fIntersections1;
    std::vector<std::vector<MaterialBudget::fPoint>> fIntersections2;
    std::vector<std::vector<Event::fVertMult>> fConfigs;
    std::vector<std::vector<MaterialBudget::fPoint>> fTracklets;
    std::vector<std::vector<double>> fVertexes;
    std::vector<double> fVertexesZ;
    std::vector<std::vector<double>> fVertexesTrackDca;
    std::vector<double> fVertexesZTrackDca;
    TH1D* fResiduals = new TH1D("Residuals", "Residuals", 500,-5,5);

    void FillHistoMinDca(TH1D* histo, vector<MaterialBudget::fPoint>& tracklets, vector<double>& vertextemp);
    void FillHistoResiduals();

 
    ClassDef(Reconstruction, 1)
};
 
 
#endif
