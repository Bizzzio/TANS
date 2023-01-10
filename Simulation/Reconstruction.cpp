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

   
    TrackletsReco();
    cout << "TrackletsReco ended" << endl;
    //MinDca();
    //cout << "MinDca ended" << endl;
    //cout << "FillHIstoResiduals ended" << endl;
    file->Close();
    VertexRecoGeom();
    cout << "MinGlobalDistance ended" << endl;
    FillHistoResiduals();
    TFile outfile("outfile.root","recreate");
    fResiduals->Write();
    outfile.Close();
    /*TFile outfile("trackletdca.root","recreate");
    ftrackletdca->Write();
    trackletdca.Close();*/
}

void Reconstruction::TrackletsReco()
{
    //cout << "Entering FindTracklets" << endl;
    double phimax = 0.010; // maximum angle difference for 2 intersections to be a tracklet
    std::vector<MaterialBudget::fPoint> tracklet;
    for (unsigned i=0; i<fIntersections1.size(); ++i)
    {
        int count=0;
        int countdiscarded=0;
        for (auto y: fIntersections1[i]){
            for (auto j: fIntersections2[i]){
                if(abs(y.phi-j.phi)<phimax)
                {
                    //if(CheckTracklet(countdiscarded,y,j)){
                        tracklet.push_back(y);
                        tracklet.push_back(j);
                        count++;
                    //}
                }
            }
        }
        //cout << "Reconstructed tracklets: " << count << endl;
        //cout << "Discarded tracklets: " << countdiscarded << endl;
        if (tracklet.size()>0)
            fTracklets.push_back(tracklet);
        tracklet.clear();
    }
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
        //cout << "Time of closest approach: " << t << endl;

        vertextemp.push_back(tracklets[j].z + c*t);
        histo->Fill(vertextemp.back());
    }
}

void Reconstruction::FillHistoResiduals()
{
    for(unsigned i=0; i<fVertexesZ.size(); ++i)
        fResiduals->Fill(fVertexesZ[i]-fConfigs[i][0].z);
}

std::vector<double> Reconstruction::GetTrackletParameters(MaterialBudget::fPoint point1, MaterialBudget::fPoint point2)
{
    // returns coordinate and angle direction of vector of line going from point1 to point2
    std::vector<double> trackletparameters;
    double aa,bb,cc,theta,phi,norm; 
    norm = TMath::Sqrt((point2.x-point1.x)*(point2.x-point1.x)+(point2.y-point1.y)*(point2.y-point1.y)+(point2.z-point1.z)*(point2.z-point1.z));
    aa = (point2.x - point1.x);
    bb = (point2.y - point1.y);
    cc = (point2.z - point1.z);
    theta=TMath::ACos(cc);
    phi=TMath::ATan(bb/aa);
    if (aa<0)
        phi+=TMath::Pi();

    trackletparameters.push_back(aa);
    trackletparameters.push_back(bb);
    trackletparameters.push_back(cc);
    trackletparameters.push_back(theta);
    trackletparameters.push_back(phi);
    return trackletparameters;
}

void Reconstruction::VertexRecoGeom()
{   
    for(int i=0; i<fTracklets.size();i++)
    {
        //if (i%1000 == 0)
            cout<<"Event "<<i<<endl;
  	    double tmin=-2., tmax=0.;
        int nbins=1000;
        TH1D* histo = new TH1D("vertex", "vertex", 50, -15, 15);             //to discard combinatorial bkg
  	    TH1D* hist = new TH1D("TrackDcaVertex", "TrackDcaVertex", nbins, tmin, tmax);       //for time of closest approach
        
        std::vector<double> vertexTemp;
        std::vector<MaterialBudget::fPoint> acceptedtracks;

        FillHistoMinDca(histo, fTracklets[i], vertexTemp);
      	int histomax(histo->GetMaximumBin());
        double xmin(histo->GetBinLowEdge(histomax)), xmax=xmin+histo->GetBinWidth(histomax);
        for (int j=0; j<vertexTemp.size(); j++)
            if (vertexTemp[j]<xmax && vertexTemp[j]>xmin)
            {
            	acceptedtracks.push_back(fTracklets[i][2*j]);
              	acceptedtracks.push_back(fTracklets[i][2*j+1]);
            }
    	std::vector<std::vector<double>> newvelocities;
  		for(int z=0; z<acceptedtracks.size(); z+=2){
    	  newvelocities.push_back(GetTrackletParameters(acceptedtracks[z],acceptedtracks[z+1]));
    	}
      	FillHistoVertexGeom(acceptedtracks, newvelocities, hist);    	
    
        //TFile file("VertexRecoGeom.root", "recreate");
        //hist->Write();
        //histo->Write();
        //file.Close();
        int histomin(hist->GetMinimumBin());
        double vertextime = hist->GetBinCenter(histomin);
        //cout<<"VERTEXTIME:"<< vertextime<<endl;
        for(int y=0; y<acceptedtracks.size(); y+=2){
            std::vector<double> vertexcoordinates = Line(acceptedtracks[y], acceptedtracks[y+1], GetTrackletParameters(acceptedtracks[y],acceptedtracks[y+1]), vertextime);
            fVertexesTrackDca.push_back(vertexcoordinates);
            fVertexesZTrackDca.push_back(vertexcoordinates[2]);
        }
        double Vertexmean=0;
        for (int i=0; i<fVertexesZTrackDca.size(); i++)
        {
            Vertexmean+=fVertexesZTrackDca[i];
        }
        Vertexmean/=fVertexesZTrackDca.size();
        fVertexesZ.push_back(Vertexmean);
        acceptedtracks.clear();
  	    delete histo;
  	    delete hist;    
    }
}		

void Reconstruction::FillHistoVertexGeom(vector<MaterialBudget::fPoint> acceptedtracks, vector<vector<double>> velocities, TH1D* graph)
{
    double mindistance=100000;
  	for(int i=1; i<=graph->GetNbinsX(); i++)
    {
        //if (graph->GetBinCenter(i)<0 && graph->GetBinCenter(i)>-2)
        //    cout << "Filling histogram cell for time " << graph->GetBinCenter(i) << endl;   
        double totaldistance = GetTotalDistance(acceptedtracks, velocities, graph->GetBinCenter(i), mindistance); // total distance among all tracklets
        if (totaldistance<mindistance)
            mindistance=totaldistance;
        graph->Fill(graph->GetBinCenter(i),totaldistance);
    }    
} 

double Reconstruction::GetTotalDistance(vector<MaterialBudget::fPoint> acceptedtracks, vector<vector<double>> velocities, double time, double mindistance)
{    
        //returns total distance of all tracklets with respect to track0 at t=time
    //cout << "Entering GetTotalDistance for track number " << trackindex << endl;
  
  double countdistance=0;
  for(unsigned z=0; z<acceptedtracks.size(); z+=2){
    std::vector<double> vec1 = Line(acceptedtracks[z], acceptedtracks[z+1], velocities[z/2], time);
    if (countdistance>mindistance)
        break;
  	for(unsigned w=z; w<acceptedtracks.size(); w+=2){
      std::vector<double> vec2 = Line(acceptedtracks[w], acceptedtracks[w+1], velocities[w/2], time);
      double dx = vec1[0]-vec2[0];
      double dy = vec1[1]-vec2[1];
      double dz = vec1[2]-vec2[2];
      countdistance += TMath::Sqrt(dx*dx + dy*dy + dz*dz);
      if (countdistance>mindistance)
        break;
    }
  }
  return countdistance;
}

vector<double> Reconstruction::Line(MaterialBudget::fPoint point1, MaterialBudget::fPoint point2, vector<double> velocity, double time)
{
    // returns point of line starting from point1 and passing for point2 at time t 
    double x, y, z;
    x = point1.x + velocity[0]*time;
    y = point1.y + velocity[1]*time;
    z = point1.z + velocity[2]*time;
    //cout << x << " " << y << " " << z << endl;
    //cout << "Line Ending" << endl;
    return {x,y,z};
}