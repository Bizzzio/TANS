#ifndef PARTICLE_H
#define PARTICLE_H

#include "TObject.h"
#include "TMath.h"
#include "Riostream.h"
#include "TRandom.h"
#include "TH1D.h"
#include <vector>

class Particle : public TObject
{
/*
 *  Class for describing particles
 */

public:
    Particle ();                                                                    //default constructor
    Particle (double x, double y, double z, double px, double py, double pz);
    Particle (double* point, double* direction);                                    //array initialisation
    Particle (vector<double> point, vector<double> direction);                      //vector initialisation
    Particle (vector<double> point, TRandom* rndmptr);                              //random initialisation
    Particle (vector<double> point, TH1D* theta);                              

    vector<double> GetDirection()       const             {return {fPx,fPy,fPz};}
    vector<double> GetPolarDirection()  const             {return {fTheta,fPhi};}
    vector<double> GetPoint()           const             {return {fX,fY,fZ};}

    Particle& SetDirection(vector<double> direction)            {return SetDirection(direction[0],direction[1],direction[2]);}
    Particle& SetDirection(double *direction)                   {return SetDirection(direction[0],direction[1],direction[2]);}
    Particle& SetDirection(double px, double py, double pz);
    Particle& SetPolarDirection(double theta, double phi)       {return SetDirection(TMath::Sin(theta)*TMath::Cos(phi), TMath::Sin(theta)*TMath::Sin(phi),TMath::Cos(theta));}
    Particle& SetPoint(vector<double> point)                    {return SetPoint(point[0],point[1],point[2]);}
    Particle& SetPoint(double *point)                           {return SetPoint(point[0],point[1],point[2]);}
    Particle& SetPoint(double x, double y, double z)            {fX=x; fY=y; fZ=z; return *this;}


private:
    double fX,fY,fZ;
    double fPx,fPy,fPz;
    double fTheta, fPhi;
    double ComputePhi(double x, double y);

    ClassDef(Particle, 1)
};

#endif