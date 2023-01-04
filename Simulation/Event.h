#ifndef EVENT_H
#define EVENT_H

#include <TSystem.h>
#include <TGeoManager.h>
#include <string>
#include "TStopwatch.h"
#include "Detector.h"

class Event : public TObject
{
/*
 *  Class for executing an event
 *  -------------------------
 *  Parameters:
 *  Detector: vector<MaterialBudget*> 
 *      Contains pointers to the detectors/material budgets used in the experiment
 * 
 *  Multiplicity: unsigned  
 *      Multiplicity of the event
 * 
 *  x,y,z: double
 *      Coordinates of the primary vertex
 * 
 *  gentree,rectree: TTree&
 *      Trees for data output
 * 
 */
public:
    
    //Struct for containing event infos: multiplicity and vertex coordinates
    using fVertMult = struct{
        unsigned int multiplicity; double x, y, z;     //no need for theta in cilindrical coordinates
        void print() {cout<<"Multiplicity': "<< multiplicity <<"\t Vertex coordinates: x: "<< x <<", y: "<< y <<", z: "<< z <<endl;}
    };

    Event(vector<MaterialBudget*> detectors, unsigned multiplicity, double x, double y, double z, TTree& gentree, TTree& rectree);
    Event(vector<MaterialBudget*> detectors, unsigned multiplicity, double x, double y, double z,vector<double> direction, TTree& gentree, TTree& rectree);
    Event(vector<MaterialBudget*> detectors, unsigned multiplicity, double x, double y, double z,TH1D* pseudorap, TTree& gentree, TTree& rectree);
    ~Event();
    void ProcessEvent(TTree&, TTree&);
    void EventVisual(vector<Particle*> particles);

private:
    vector<Particle*>           fParticles;
    vector<MaterialBudget*>     fDetectors;
    vector<fVertMult>           fConfig;
    vector<fVertMult>*          fConfigprt = &fConfig;
    
    ClassDef(Event, 1)

};

#endif