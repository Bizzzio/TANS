#ifndef DETECTOR_H
#define DETECTOR_H

#include "Particle.h"
#include "MaterialBudget.h"

class Detector : public MaterialBudget
{
/*
 *  Class for describing a cylindrical detector
 *  -------------------------
 *  Parameters:
 *  thickness: double 
 *      Thickness of the layer
 * 
 *  radius: double 
 *      Radius of the layer
 * 
 *  lenght: double 
 *      Lenght of the layer
 * 
 *  density: double 
 *      Density of the layer
 * 
 *  z: int 
 *      Atomic number of the layer's material
 * 
 *  a: int 
 *      Mass number of the layer's material
 * 
 *  material: string
 *      String containing material symbol
 *      Currently only "Si" and "Be" are supported
 *  
 *  multscat: bool
 *      Enables multiple scattering
 * 
 *  smearing: bool
 *      Activates smearing
 * 
 *  noise: bool
 *      Activates noise
 * 
 */
public:

    Detector();
    Detector(double thickness, double radius, double length, string material, bool multscat, bool smearing, bool noise, int meannoise);

    virtual bool IsDetector() const                             {return true;}
    virtual void Interaction(Particle* particle);

    Detector& SetStatus(bool smearing, bool noise)              {fSmearing = smearing; fNoise = noise; return *this;}
    Detector& SetStatus(vector<bool> status);

    virtual void ClearData()                                    {fTrueHit.clear(); fRecoHit.clear(); if(fNoise) Noise();}
    virtual void FillData(Particle*);
    virtual void SetBranchAddress(TTree& gentree,TTree& rectree, unsigned countdet);
    MaterialBudget::fPoint GetSmearedIntersection(const MaterialBudget::fPoint& intersection);
    virtual MaterialBudget::fPoint GetLastIntersection(const Particle* particle)  {return fTrueHit.back();}
    vector<bool> GetStatus() const                              {return {fMultScat, fSmearing, fNoise};}
    vector<MaterialBudget::fPoint> GetTrueHits() const          {return fTrueHit;}
    vector<MaterialBudget::fPoint> GetRecoHits() const          {return fRecoHit;}
    
private:
    void Noise();
    
    bool fSmearing, fNoise;
    int fMeanNoise;
    vector<MaterialBudget::fPoint> fTrueHit, fRecoHit;
    vector<MaterialBudget::fPoint>* fTrueHitPtr = &fTrueHit;
    vector<MaterialBudget::fPoint>* fRecoHitPtr = &fRecoHit;

    double fMuAngular=0, fSigmaAngular=3e-3;    // Parameters for smearing (in cm)
    double fMuZ=0, fSigmaZ=1.2e-2;              // Parameters for smearing (in cm)

};

#endif