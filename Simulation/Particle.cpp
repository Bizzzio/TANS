#include "Particle.h"

ClassImp(Particle)

Particle::Particle() :
fX(0.), fY(0.), fZ(0.), fPx(1.), fPy(0.), fPz(0.), fTheta(0.), fPhi(0.) {}

Particle::Particle(double x, double y, double z, double px, double py, double pz)
{
    SetDirection(px,py,pz).SetPoint(x,y,z);
}

Particle::Particle (double* point, double* direction)
{
    SetDirection(direction[0],direction[1],direction[2]).SetPoint(point[0],point[1],point[2]);
}  

Particle::Particle (vector<double> point, vector<double> direction)
{
    SetDirection(direction[0],direction[1],direction[2]).SetPoint(point[0],point[1],point[2]);
}  

Particle::Particle (vector<double> point, TRandom* rndmptr)
{
    double direction[3];
    rndmptr->Sphere(direction[0],direction[1],direction[2],1);
    SetDirection(direction[0],direction[1],direction[2]).SetPoint(point[0],point[1],point[2]);
}

Particle::Particle (vector<double> point, TH1D* theta)
{
    SetPolarDirection(theta->GetRandom(),2*TMath::Pi()*gRandom->Rndm()).SetPoint(point[0],point[1],point[2]);
}

Particle& Particle::SetDirection(double px, double py, double pz)
{
    double norm= px*px + py*py + pz*pz;
    if (TMath::Abs(norm - 1) > 1e-6 && norm > 1e-6)
    {
        norm=TMath::Sqrt(norm);
        cout<< "\033[93mWarning: trying to set direction to ("<<px<<", "<<py<<", "<<pz<<"), which is not normalized, setting the normalised vector ";
        px/=norm; py/=norm; pz/=norm;
        cout<<"("<<px<<", "<<py<<", "<<pz<<") instead \033[0m \n";
    }
    else if (norm < 1e-6)
    {
        cout<< "\033[93mWarning: trying to set direction to null vector, setting the normalised vector to (1, 0, 0) instead\033[0m \n";
        px=1; py=0; pz=0;
    }

    fPx=px; fPy=py; fPz=pz;
    fTheta=TMath::ACos(fPz);
    fPhi=ComputePhi(fPx,fPy);

    return *this;
}

double Particle::ComputePhi(double x, double y)
{
/*
 *  Helper function that evauates the azimuth
 *  -------------------------
 *  Parameters:
 *  x,y: double
 *      coordinates of a point
 * 
 *  -------------------------
 *  Returns:
 *      double containing azimuth of the point of coordinates (x,y)
 * 
 */
    double phi = TMath::ATan(y / x);
    if (x<0)
        phi+=TMath::Pi();
    return phi;
}