#include "Reconstruction.h"
#include"TDirectory.h"
#include "TH1D.h"
#include "TAxis.h"
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
    //MinDca();
    //cout << "MinDca ended" << endl;
    //FillHistoResiduals();
    //cout << "FillHIstoResiduals ended" << endl;
    MinGlobalDistance();
    cout << "MinGlobalDistance ended" << endl;
    file->Close();
    //TFile outfile("outfile.root","recreate");
    //fResiduals->Write();
    //outfile.Close();
    /*TFile outfile("trackletdca.root","recreate");
    ftrackletdca->Write();
    trackletdca.Close();*/
}

/*bool Reconstruction::CheckTracklet(MaterialBudget::fPoint int1, Materialbudget::fPoint int2)
{

}*/

void Reconstruction::TrackletsReco()
{
    //cout << "Entering FindTracklets" << endl;
    double phimax = 0.010; // maximum angle difference for 2 intersections to be a tracklet
    std::vector<MaterialBudget::fPoint> tracklet;
    for (unsigned i=0; i<fIntersections1.size(); ++i)
    {
        int count=0;
        for (auto y: fIntersections1[i]){
            for (auto j: fIntersections2[i]){
                if(abs(y.phi-j.phi)<phimax)
                {
                    //if(CheckTracklet(y,j)){
                        tracklet.push_back(y);
                        tracklet.push_back(j);
                        count++;
                    //}
                }
            }
        }
        cout << "Reconstructed tracklets: " << count << endl;
        if (tracklet.size()>0)
            fTracklets.push_back(tracklet);
        tracklet.clear();
    }
    //std::vector<std::vector<int>> CombinatorialTracklets;
    //for(unsigned event=0; event<fTracklets.size(); event++){
    //    std::vector<int> combtrackindexes; 
    //    for(unsigned j=0; j<fTracklets[event].size(); j+=2){
    //        for(unsigned k=0; k<fTracklets[event].size(); k+=2){
    //            if(k==j) continue;
    //            if(fTracklets[event][j].x==fTracklets[event][k].x && fTracklets[event][j].y==fTracklets[event][k].y 
    //                && fTracklets[event][j].z==fTracklets[event][k].z)
    //                {
    //                    combtrackindexes.push_back(1); // tracklets sharing intersection with first detector
    //                    combtrackindexes.push_back(j);
    //                    combtrackindexes.push_back(k);
    //                }
    //            if(fTracklets[event][j+1].x==fTracklets[event][k+1].x && fTracklets[event][j+1].y==fTracklets[event][k+1].y 
    //                && fTracklets[event][j+1].z==fTracklets[event][k+1].z)
    //                {   
    //                    combtrackindexes.push_back(2); // tracklets sharing intersection with second detector
    //                    combtrackindexes.push_back(j+1);
    //                    combtrackindexes.push_back(k+1);
    //                }
    //        }
    //    }
    //    CombinatorialTracklets.push_back(combtrackindexes);
    //    combtrackindexes.clear();
    //}
    //cout << "Combinatorial tracklets: " << CombinatorialTracklets.size()/3 << endl;
    //for(unsigned i=0; i<fTracklets.size(); i++){
    //    cout << "Event " << i << " tracklets: " << endl;
    //    for(unsigned j=0; j<fTracklets[i].size(); j+=2){
    //        fTracklets[i][j].print();
    //        fTracklets[i][j+1].print();
    //        cout << endl;
    //    }
    //}
}

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
        t = - (a*tracklets[j].x + b*tracklets[j].y) / (a*a+b*b);    // t of closest approach
        cout << "Time of closest approach: " << t << endl;

        vertextemp.push_back(tracklets[j].z + c*t);
        histo->Fill(vertextemp.back());
    }
}

void Reconstruction::FillHistoResiduals()
{
    for(unsigned i=0; i<fVertexesZ.size(); ++i)
        fResiduals->Fill(fVertexesZ[i]-fConfigs[i][0].z);
}

std::vector<double> Reconstruction::GetTrackletParameters(int nevent, int ntracklet)
{
    // returns parameters of tracklet n°ntracklet belonging to the event n°nevent
    std::vector<double> trackletparameters;
    double a,b,c,theta,phi; // parameters for line connecting two intersections of the same tracklet from detector 2 to detector 1
    a = fTracklets[nevent][ntracklet].x - fTracklets[nevent][ntracklet+1].x;
    b = fTracklets[nevent][ntracklet].y - fTracklets[nevent][ntracklet+1].y;
    c = fTracklets[nevent][ntracklet].z - fTracklets[nevent][ntracklet+1].z;
    theta=TMath::ACos(c);
    phi=TMath::ATan(b/a);
    if (a<0)
        phi+=TMath::Pi();

    trackletparameters.push_back(a);
    trackletparameters.push_back(b);
    trackletparameters.push_back(c);
    trackletparameters.push_back(theta);
    trackletparameters.push_back(phi);
    return trackletparameters;
}

vector<double> Reconstruction::Line(int nevent, int index, double time)
{
    // returns point of line starting from intersection with detector 2 at t=time of tracklet n°index belonging to the event n°nevent  
    double x, y, z;
    x = fTracklets[nevent][index+1].x + GetTrackletParameters(nevent, index)[0]*time;
    y = fTracklets[nevent][index+1].y + GetTrackletParameters(nevent, index)[1]*time;
    z = fTracklets[nevent][index+1].z + GetTrackletParameters(nevent, index)[2]*time;
    //cout << x << " " << y << " " << z << endl;
    vector<double> linecoordinates;
    linecoordinates.push_back(x);
    linecoordinates.push_back(y);
    linecoordinates.push_back(z);
    //cout << "Line Ending" << endl;
    return linecoordinates;
}

void Reconstruction::FindShortestTracklet(int nevent, int& index, double& thetamin)
{
    if(fTracklets[nevent][0].z>0){
        thetamin = TMath::Pi()-GetTrackletParameters(nevent,0)[3];
    }
    if(fTracklets[nevent][0].z<=0){
        thetamin = GetTrackletParameters(nevent,0)[3];
    }
    for(unsigned i=2; i<fTracklets[nevent].size(); i+=2){
        double theta1 = GetTrackletParameters(nevent,i)[3];
        if(fTracklets[nevent][i].z<0){
            if(theta1<TMath::Pi()/2){
                if(theta1<thetamin){
                    thetamin=theta1;
                    index=i;
                }
            }
        }
        else{
            if(theta1>=TMath::Pi()/2){
                if(TMath::Pi()-theta1<thetamin){
                    thetamin=theta1;
                    index=i;
                }
            }
        }
    }
}

void Reconstruction::MinGlobalDistance()
{
    for(int nevent=0; nevent<fTracklets.size(); nevent++){
        cout << "Entering MinGlobalDistance" << endl;
            
            // First step: find tracklet with smallest theta
        int index=0;
        double theta;
        FindShortestTracklet(nevent, index, theta);
            // Second step: calculate norm of difference vector between equivalent points belonging to the octant x,y,z>0 of intersections
            // with 2nd detector of all tracklets and tracklet with minimum theta
        std::vector<double> distances;
        MaterialBudget::fPoint point0 = FirstOctant(fTracklets[nevent][index+1]); // intersection in first octant with 2nd detector for tracklet with minimum theta 
        for(unsigned i=0; i<fTracklets[nevent].size(); i+=2){
            MaterialBudget::fPoint point1 = FirstOctant(fTracklets[nevent][i+1]);
            double DX = point1.x-point0.x; 
            double DY = point1.y-point0.y;
            double DZ = point1.z-point0.z;
            distances.push_back(TMath::Sqrt(DX*DX + DY*DY + DZ*DZ));
        }
            // Third step: connecting spacial delays with time delays
        std::vector<double> delays;
        for(auto i: distances){
            delays.push_back(i);
            delays.push_back(i);
        }

        // Fourth step: declare histogram with time on x axis, total distance of all tracklet points on y axis
        int nbins = 1000;
        TH1D* hist = new TH1D("TrackDcaVertex", "TrackDcaVertex", nbins, -2, 2);
        
        // Fifth step: fill histogram
        for(unsigned int i=1; i<=nbins; i++){
            hist->Fill(hist->GetBinCenter(i),FillHistoTrackMinDca(nevent, hist->GetBinCenter(i), delays));
        }

        int binmin = hist->GetMinimumBin();
        double t0 = hist->GetXaxis()->GetBinCenter(binmin);
        cout << "Time minimizing total distance and corresponding bin: " << t0 << " " << binmin << endl;

        for(unsigned i=0; i<fTracklets[nevent].size(); i+=2){
            std::vector<double> trackletvertex = Line(nevent, i, t0+delays[i]);
            fVertexesTrackDca.push_back(trackletvertex);
            fVertexesZTrackDca.push_back(trackletvertex[2]);
        }

        cout << "Reconstructed vertexes: " << endl;
        for(unsigned b=0; b<fVertexesZTrackDca.size(); b++){
            cout << fVertexesZTrackDca[b] << "    ";
        }

        TFile TrackDca("TrackDca.root","recreate");
        hist->Write();
        TrackDca.Close();

        delete hist;
    }
}

double Reconstruction::FillHistoTrackMinDca(int nevent, double time, vector<double> delays)
{
    // fills histogram with total distance between tracklets at time t
    cout << "Entering FillHistoTrackMinDca" << endl;
    
        // Sixth step: loop over all tracklets to find total distance of all tracklets with respect to all tracklets
    double totaldistance=0; // total distance
    
    for(unsigned w=0; w<fTracklets[nevent].size(); w+=2){
        // sum total distance of all tracklets with respect to the w-th tracklet of the nevent-th event
        totaldistance += GetTotalDistance(nevent, w, time+delays[w]);
    }
    cout << "Filling histogram entry for time: " << time << endl;
    return totaldistance;
}

double Reconstruction::GetTotalDistance(int nevent, unsigned trackindex, double time)
{    
        //returns total distance of all tracklets with respect to track0 at t=time
    //cout << "Entering GetTotalDistance for track number " << trackindex << endl;
    std::vector<double> vec1 = Line(nevent,trackindex,time);
    
    double countdistance = 0.;
    for(unsigned z=0; z<fTracklets[nevent].size(); z+=2){
        // loop over all tracklets to sum distance between track0 point and their line point at time t 
        std::vector<double> vec2 = Line(nevent,z,time);
        double dx = vec1[0]-vec2[0];
        double dy = vec1[1]-vec2[1];
        double dz = vec1[2]-vec2[2];
        //cout << dx << " " << dy << " " << dz << endl;
        countdistance += TMath::Sqrt(dx*dx + dy*dy + dz*dz);
    }
    return countdistance;

}

MaterialBudget::fPoint Reconstruction::FirstOctant(MaterialBudget::fPoint point)
{
        // turns generic point into a point in x,y,z>0 octant on the cone defined by the same theta
    cout << "Entering FirstOctant" << endl;
    
    double x = point.x; 
    double y = point.y; 
    double z = point.z;

    if(x>0 && y>0 && z<0) point.z = -z;
    if(x>0 && y<0 && z>0) point.y = -y;
    if(x<0 && y>0 && z>0) point.x = -x;
    if(x<0 && y<0 && z>0) point.x = -x; point.y = -y;
    if(x<0 && y>0 && z<0) point.x = -x; point.z = -z; 
    if(x>0 && y<0 && z<0) point.y = -y; point.z = -z; 
    if(x<0 && y<0 && z<0) point.x = -x; point.y = -y; point.z = -z;

    return point;
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