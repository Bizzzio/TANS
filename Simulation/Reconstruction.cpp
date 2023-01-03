#include "Reconstruction.h"
#include"TDirectory.h"
#include "TH1D.h"
#include "TCanvas.h"
 
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
    
    
    TrackletsReco();
    cout << "TrackletsReco ended" << endl;
    MinDca();
    cout << "MinDca ended" << endl;
    FillHistoResiduals();
    cout << "FillHIstoResiduals ended" << endl;
    file->Close();
    TFile outfile("outfile.root","recreate");
    fResiduals->Write();
    outfile.Close();
}
 
void Reconstruction::TrackletsReco()
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
        tracklet.clear();
    }
}
 
std::vector<double> Reconstruction::GetTrackletParameters(int j)
{
    std::vector<double> trackletparameters;
    double a,b,c,theta,phi; // parameters for line connecting two intersections of the same tracklet from detector 2 to detector 1
    a = fTracklets[j][0].x - fTracklets[j][1].x;
    b = fTracklets[j][0].y - fTracklets[j][1].y;
    c = fTracklets[j][0].z - fTracklets[j][1].z;
    double norm = TMath::Sqrt(a*a + b*b + c*c);
    double aa, bb, cc;
    aa=a/norm;
    bb=b/norm;
    cc=c/norm;

    theta=TMath::ACos(cc);
    phi=TMath::ATan(bb/aa);
    if (aa<0)
        phi+=TMath::Pi();

    trackletparameters.push_back(a);
    trackletparameters.push_back(b);
    trackletparameters.push_back(c);
    trackletparameters.push_back(theta);
    trackletparameters.push_back(phi);
    return trackletparameters;
}

/*void Reconstruction::VertexReco(int index, double t)
{
    cout << "Entering VertexReco" << endl;
    std::vector<double> vertex;
    double a,b,c; // parameters for line connecting two intersections of the same tracklet from detector 2 to detector 1
    a = fTracklets[index][0].x - fTracklets[index][1].x;
    b = fTracklets[index][0].y - fTracklets[index][1].y;
    c = fTracklets[index][0].z - fTracklets[index][1].z;

    // value of parameter t for intersection between line and orthogonal plane passfTracklets[index]ng for O
    t = -(a*fTracklets[index][1].x + b*fTracklets[index][1].y + c*fTracklets[index][1].z)/(a*a + b*b + c*c);

    // vertex coordfTracklets[index]nates
    double x, y, z;
    x = fTracklets[index][1].x + a*t;
    y = fTracklets[index][1].y + b*t;
    z = fTracklets[index][1].z + c*t;

    vertex.push_back(x);
    vertex.push_back(y);
    vertex.push_back(z);

    fVertexes.push_back(vertex);
    vertex.clear();

    cout << "Vertex coordinates: ";
    for(auto i: vertex){
            cout << i << " ";
        }
    cout << endl;
}*/

void Reconstruction::MinDca()
{
    TH1D* histo = new TH1D("vertex", "vertex", 50, -15,15);
    vector<double> vertexTemp;
    for(auto& i: fTracklets)
    {
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
    for(unsigned i=0; i<fVertexesZ.size(); ++i)
        fResiduals->Fill(fVertexesZ[i]-fConfigs[i][0].z);
}

vector<double> Reconstruction::Line(int index, double time)
{
    double x = fTracklets[index][1].x + GetTrackletParameters(index)[0]*time;
    double y = fTracklets[index][1].y + GetTrackletParameters(index)[1]*time;
    double z = fTracklets[index][1].z + GetTrackletParameters(index)[2]*time;
    vector<double> linecoordinates;
    linecoordinates.push_back(x);
    linecoordinates.push_back(y);
    linecoordinates.push_back(z);
    return linecoordinates;
}

void Reconstruction::MinGlobalDistance()
{
    int index=0;
    int phi=abs(fTracklets[0][0].phi);
    for(unsigned i=1; i<fTracklets.size(); ++i){
        if(abs(fTracklets[i][0].phi)<phi){
            phi=fTracklets[i][0].phi;
            index=i;
        }
    }
    vector<double> distances;
    MaterialBudget::fPoint point0 = FirstOctant(fTracklets[index][1]);
    for(unsigned i=0; i<fTracklets.size(); ++i){
        // if(i==index) distances.push_back(0.);
        MaterialBudget::fPoint point1 = FirstOctant(fTracklets[i][1]);
        double DX = point1.x-point0.x; 
        double DY = point1.y-point0.y;
        double DZ = point1.z-point0.z;
        distances.push_back(TMath::Sqrt(DX*DX + DY*DY + DZ*DZ));
    }
    vector<double> times;
    for(auto i: distances){
        times.push_back(i/(3.*10*10*10*10*10*10*10*10));
    }
    double totaldistance = 0.;
    double t;
    for(unsigned y=0; y<fTracklets.size(); y++){
        for(unsigned j=0; j<fTracklets.size(); j++){
            double dx = Line(y,t+times[y])[0]-Line(j,t+times[j])[0];
            double dy = Line(y,t+times[y])[1]-Line(j,t+times[j])[1];
            double dz = Line(y,t+times[y])[2]-Line(j,t+times[j])[2];
            totaldistance += TMath::Sqrt(dx*dx + dy*dy + dz*dz);
        }
    }
}

vector<double> Reconstruction::TimesMinGlobalDistance(int index0)
{
    double a,b,c; // parameters for line connecting two intersections of the same tracklet from detector 2 to detector 1
    a = fTracklets[index0][0].x - fTracklets[index0][1].x;
    b = fTracklets[index0][0].y - fTracklets[index0][1].y;
    c = fTracklets[index0][0].z - fTracklets[index0][1].z;
    double t;
    t = - (a*fTracklets[index0][1].x + b*fTracklets[index0][1].y) / (a*a+b*b);    //t of closest approach
}

MaterialBudget::fPoint Reconstruction::FirstOctant(MaterialBudget::fPoint point) // turns point into a vector in x,y,z>0 octant with minimum phi
{
    double x = point.x; 
    double y = point.y; 
    double z = point.z;

    if(point.x>0 && point.y>0 && point.z<0) point.z = -z;
    if(point.x>0 && point.y<0 && point.z>0) point.y = -y;
    if(point.x<0 && point.y>0 && point.z>0) point.x = -x;
    if(point.x<0 && point.y<0 && point.z>0) point.x = -x; point.y = -y;
    if(point.x<0 && point.y>0 && point.z<0) point.x = -x; point.z = -z; 
    if(point.x>0 && point.y<0 && point.z<0) point.y = -y; point.z = -z; 
    if(point.x<0 && point.y<0 && point.z<0) point.x = -x; point.y = -y; point.z = -z;

    return point;
}
