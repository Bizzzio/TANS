#include "Event.h"

//Used to create ACLiC dictionary for vector<Event::fVertMult>
#ifdef __MAKECINT__
#pragma link C++ class vector<Event::fVertMult>+;
#endif

ClassImp(Event)

Event::Event(vector<MaterialBudget*> detectors, unsigned multiplicity, double x, double y, double z, TTree& gentree, TTree& rectree):
fDetectors(detectors)
{

    while (fParticles.size()<multiplicity)          //Creates particles
    {
        Particle* temp = new Particle({x,y,z},gRandom);
        fParticles.push_back(temp);
    }

    //Stores event infos
    fVertMult temp;
    temp.multiplicity = multiplicity;
    temp.x = x;
    temp.y = y;
    temp.z = z;
    fConfig.push_back(temp);
    

    ProcessEvent(gentree, rectree);
    //EventVisual(fParticles);
}

Event::Event(vector<MaterialBudget*> detectors, unsigned multiplicity, double x, double y, double z,vector<double> direction, TTree& gentree, TTree& rectree):
fDetectors(detectors)
{

    while (fParticles.size()<multiplicity)          //Creates particles
    {
        Particle* temp = new Particle({x,y,z},direction);
        fParticles.push_back(temp);
    }

    //Stores event infos
    fVertMult temp;
    temp.multiplicity = multiplicity;
    temp.x = x;
    temp.y = y;
    temp.z = z;
    fConfig.push_back(temp);
    

    ProcessEvent(gentree, rectree);
    //EventVisual(fParticles);
}

Event::Event(vector<MaterialBudget*> detectors, unsigned multiplicity, double x, double y, double z,TH1D* pseudorap, TTree& gentree, TTree& rectree):
fDetectors(detectors)
{

    while (fParticles.size()<multiplicity)          //Creates particles
    {
        Particle* temp = new Particle({x,y,z},pseudorap);
        fParticles.push_back(temp);
    }

    //Stores event infos
    fVertMult temp;
    temp.multiplicity = multiplicity;
    temp.x = x;
    temp.y = y;
    temp.z = z;
    fConfig.push_back(temp);
    

    ProcessEvent(gentree, rectree);
    //EventVisual(fParticles);
}

Event::~Event()
{
    //deallocates memory
    for (auto&i : fParticles)
        delete i;
}

void Event::ProcessEvent(TTree& gentree, TTree& rectree)
{
/*
 *  Function to handle the event
 *  -------------------------
 *  Parameters:
 *  gentree: TTree&
 *      TTree of generated events
 *  rectree: TTree&
 *      TTree of reconstructed events
 */
    gentree.SetBranchAddress("Config", &fConfigprt);
    for(vector<MaterialBudget*>::size_type j = 0; j<fDetectors.size(); j++)
    {
        for (vector<Particle*>::size_type i = 0; i<fParticles.size(); i++)
            fDetectors[j]->Interaction(fParticles[i]);

    }
    gentree.Fill();
    rectree.Fill();
    for (auto& j : fDetectors)
        j->ClearData();
}

void Event::EventVisual(vector<Particle*> particles)
{
/*
 *  Function to visualise the event
 *  -------------------------
 *  Parameters:
 *  particles: vector<Particle*>
 *      vector of generated particles
 */
    gSystem->Load("libGeom");

    TGeoManager *geom = new TGeoManager("world", "the simplest geometry");      //TODO: deallocate memory
    TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
    TGeoMaterial *Fe = new TGeoMaterial("Fe",55.845,26,7.87);
    TGeoMedium   *med = new TGeoMedium("Vacuum",1,mat);
    TGeoMedium *Iron = new TGeoMedium("Iron",0,Fe);
    
    TGeoVolume *top=gGeoManager->MakeBox("Top",med,15.,15.,30.);
    TGeoVolume *Outer=gGeoManager->MakeTubs("Outer",Iron,6.99,7.01,27.,90,360);
    TGeoVolume *Inner=gGeoManager->MakeTubs("Inner",Iron,3.99,4.01,27.,90,360);
    TGeoVolume *BeamPipe=gGeoManager->MakeTubs("BeamPipe",Iron,2.96,3.04,27.,90,360);

    //double phi =32, theta=63;           //gamma = phi(polar)-90, theta(here)=theta(polar)+90
    //double rad_phi= TMath::DegToRad()*phi, rad_theta=TMath::DegToRad()*theta;
    //double len1=8/TMath::Sin(rad_theta)/2;
    std::string Trackstr = "Track";
    unsigned count=0;
    
    top->SetLineColor(kBlue);
    gGeoManager->SetTopVisible(0); // the TOP is invisible
    Outer->SetLineColor(kOrange);
    top->AddNode(Outer,1);
    Inner->SetLineColor(kBlue);
    top->AddNode(Inner,1);
    BeamPipe->SetLineColor(kMagenta);
    BeamPipe->SetFillStyle(1001);
    BeamPipe->SetFillColor(kMagenta);
    //top->SetVisibility(kFALSE);
    top->AddNode(BeamPipe,1);

    for (auto& i : particles)
    {
        vector<double> dir = i->GetPolarDirection(); 
        dir[0] = 180 - dir[0]*TMath::RadToDeg();
        dir[1] = dir[1]*TMath::RadToDeg() - 90;

        double critangle=TMath::ATan(20/27)*TMath::RadToDeg();
        critangle=14;

        if(dir[0] < (180 - critangle)  && dir[0]>critangle)
        {
            double rad_phi= TMath::DegToRad()*dir[1], rad_theta=TMath::DegToRad()*dir[0];
            double len1=8/TMath::Sin(rad_theta)/2;
            TGeoVolume *Track=gGeoManager->MakeTubs((Trackstr+std::to_string(count)).c_str(),Iron,0,0.05,len1,0,360);
            top->AddNode(Track,count, new TGeoCombiTrans(-len1*TMath::Sin(rad_phi)*TMath::Sin(rad_theta),len1*TMath::Cos(rad_phi)*TMath::Sin(rad_theta),-len1*TMath::Cos(rad_theta),new TGeoRotation("ha_8",dir[1],dir[0],0))); // 45,90,0
            ++count;
        }
    }
    
    gGeoManager->SetTopVolume(top);
    gGeoManager->SetVisLevel(3);
    gGeoManager->CloseGeometry();
    top->Draw("ogl");

}