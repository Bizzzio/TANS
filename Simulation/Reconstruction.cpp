#include "Reconstruction.h"
#include"TDirectory.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <cmath>
 
ClassImp(Reconstruction)
 
Reconstruction::Reconstruction()
{
    // Open the file containing the tree.
    TFile* file = TFile::Open("Simulation/Tree.root");
    TTreeReader RecReader("fTreeRec", file);
   
    // The branch "RecHits detector 1" contains MaterialBudget::fPoint.

    TTreeReaderValue<std::vector<MaterialBudget::fPoint>> intersect1(RecReader, "RecHits_1");
    TTreeReaderValue<std::vector<MaterialBudget::fPoint>> intersect2(RecReader, "RecHits_2");
 
    // Loop over all entries of the TTree or TChain.
    while (RecReader.Next()) {
       // Just access the data as if intersect1 and intersect2 were iterators (note the '*'
       // in front of them):
       fIntersections1.push_back(*intersect1);
       fIntersections2.push_back(*intersect2);
    }

    TTreeReader GenReader("fTreeGen", file);
   
    TTreeReaderValue<std::vector<Event::fVertMult>> configs(GenReader, "Config");
 
    // Loop over all entries of the TTree or TChain.
    while (GenReader.Next()) {
       // Just access the data as if intersect1 and intersect2 were iterators (note the '*'
       // in front of them):
       fConfigs.push_back(*configs);
    }



    //cout << "Intersezioni1: " << fIntersections1.size() << endl;    
    //for (auto i: fIntersections1){
    //    //cout << i.size() << endl;
    //    for (auto j: i){
    //        //j.print();
    //    }
    //    //cout << endl;
    //}
    //cout << endl;
    //cout << "Intersezioni2: " << fIntersections2.size() << endl;
    //for (auto x: fIntersections2){
    //    cout << x.size() << endl;
    //    for (auto y: x){
    //        //y.print();
    //    }
    //    //cout << endl;
    //}
    //cout << endl;
    //cout << "Fine stampa" << endl;
    
    
    FindTracklets();
    cout << "Fine FindTracklets ended" << endl;
    MinDca();
    cout << "Fine MinDca" << endl;
    FillHistoResiduals();
    cout << "FillHIstoResiduals ended" << endl;
    file->Close();
    TFile outfile("outfile.root","recreate");
    fResiduals->Write();
    outfile.Close();
}
 
void Reconstruction::FindTracklets()
{
    //cout << "Entering FindTracklets" << endl;
    double phimax = 0.010; // maximum angle difference for 2 intersections to be a tracklet
    std::vector<MaterialBudget::fPoint> tracklet;
    for (unsigned i=0; i<fIntersections1.size(); ++i)
    {
        for (auto y: fIntersections1[i])
            for (auto j: fIntersections2[i])
                if(abs(y.phi-j.phi)<phimax)
                {
                    tracklet.push_back(y);
                    tracklet.push_back(j);
                }
        if (tracklet.size()>0)
            fTracklets.push_back(tracklet);
        else
            fTracklets.push_back({});
        tracklet.clear();
    }
}
 
void Reconstruction::VertexReco()
{
    cout << "Entering VertexReco" << endl;
    std::vector<double> vertex;
    for(auto i: fTracklets){
        double a,b,c; // parameters for line connecting two intersections of the same tracklet from detector 2 to detector 1
        a = i[0].x - i[1].x;
        b = i[0].y - i[1].y;
        c = i[0].z - i[1].z;
 
        // value of parameter t for intersection between line and orthogonal plane passing for O
        double t = -(a*i[1].x + b*i[1].y + c*i[1].z)/(a*a + b*b + c*c);
 
        // vertex coordinates
        double x, y, z;
        x = i[1].x + a*t;
        y = i[1].y + b*t;
        z = i[1].z + c*t;
 
        vertex.push_back(x);
        vertex.push_back(y);
        vertex.push_back(z);
 
        fVertexes.push_back(vertex);
        vertex.clear();
    }
    for(auto i: fVertexes){
        cout << "Vertex coordinates: ";
        for(auto j: i){
            cout << j << " ";
        }
        cout << endl;
    }
}

void Reconstruction::MinDca()
{
    TH1D* histo = new TH1D("vertex", "vertex", 50, -15,15);
    vector<double> vertexTemp;
    for(auto& i: fTracklets)
    {
        if (i.size()==0)
        {
            fVertexesZ.push_back(nan(""));
            continue;
        }
        FillHistoMinDca(histo, i, vertexTemp);

        int histomax(histo->GetMaximumBin()), count(0);
        double xmin(histo->GetBinLowEdge(histomax)), xmax=xmin+histo->GetBinWidth(histomax), mean(0);
        
        for (auto& i : vertexTemp)
            if (i<xmax && i>xmin)
            {
                mean+=i;
                ++count;
            }
        fVertexesZ.push_back(mean/count);
        //cout<<"Vertex is:"<<fVertexesZ.back();
        histo->Reset();
        vertexTemp.clear();
    }
    delete histo;
}

void Reconstruction::FillHistoMinDca(TH1D* histo, vector<MaterialBudget::fPoint>& tracklets, vector<double>& vertextemp)
{
    for (unsigned j=0;j<tracklets.size();j+=2)
    {
        double a,b,c; // parameters for line connecting two intersections of the same tracklet from detector 2 to detector 1
        a = tracklets[j+1].x - tracklets[j].x;
        b = tracklets[j+1].y - tracklets[j].y;
        c = tracklets[j+1].z - tracklets[j].z;

        double t;
        t = - (a*tracklets[j].x + b*tracklets[j].y) / (a*a+b*b);    //t of closest approach

        vertextemp.push_back(tracklets[j].z + c*t);
        histo->Fill(vertextemp.back());
    }
}

void Reconstruction::FillHistoResiduals()
{
    for(unsigned i=0; i<fConfigs.size(); ++i)
    {
        //cout<<fVertexesZ[i]<<"\t"<<fConfigs[i][0].z<<endl;
        if (!isnan(fVertexesZ[i]))
            fResiduals->Fill(fVertexesZ[i]-fConfigs[i][0].z);
    }
}