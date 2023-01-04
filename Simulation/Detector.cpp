#include "Detector.h"

Detector::Detector() : MaterialBudget() {};

Detector::Detector(double thickness, double radius, double length, string material, bool multscat=1, bool smearing=1, bool noise=0, int meannoise=0)
: MaterialBudget(thickness, radius, length, material, multscat), fMeanNoise(meannoise)
{
    SetStatus(smearing, noise);
    if(fNoise)
        Noise();
}

Detector& Detector::SetStatus(vector<bool> status)
{
    if(status.size() == 2){  
        fSmearing = status[0];
        fNoise = status[1];
    }
    else{
        cout << "Invalid features for detector, switching on smearing and noise phenomena" << endl;
        fSmearing = true;
        fNoise = true;
    }
    return *this;
}

void Detector::Interaction(Particle* particle)
{
/*
 *  Function that handles the interaction between the material budget and a particle
 *  -------------------------
 *  Parameters:
 *  part: Particle*
 *      Pointer to particle
 * 
 */
    FillData(particle); 
    if(fMultScat)
        MultScattering(particle);
}

void Detector::FillData(Particle* particle)
{
/*
 *  Function that fills the vectors containing the intersections
 *  -------------------------
 *  Parameters:
 *  part: Particle*
 *      Pointer to particle
 * 
 */
    MaterialBudget::fPoint intersection = GetIntersection(particle);
    if (intersection.isIntersection){
        fTrueHit.push_back(intersection);
        if(fSmearing){
            MaterialBudget::fPoint smearedintersection = GetSmearedIntersection(intersection);
            if(smearedintersection.isIntersection){
                fRecoHit.push_back(smearedintersection);
            }
        }
    }
}

MaterialBudget::fPoint Detector::GetSmearedIntersection(const MaterialBudget::fPoint& intersection)
{
/*
 *  Function that returns the smeared intersection with a particle
 *  -------------------------
 *  Parameters:
 *  intersection: const MaterialBudget::fPoint&
 *      intersection point
 * 
 *  -------------------------
 *  Returns:
 *      MaterialBudget::fPoint containing smeared intersection
 * 
 */
    MaterialBudget::fPoint SmearedIntersection;
    
    double ar=gRandom->Gaus(fMuAngular,fSigmaAngular);
    double zr=gRandom->Gaus(fMuZ,fSigmaZ);

    //smearing
    SmearedIntersection.phi = intersection.phi + ar/fRadius;
    SmearedIntersection.z = intersection.z + zr;

    if (SmearedIntersection.z < -fLength/2 || SmearedIntersection.z > fLength/2){
        SmearedIntersection.isIntersection = false;
    }
    else{
        SmearedIntersection.isIntersection=true;
        SmearedIntersection.x = fRadius*TMath::Cos(SmearedIntersection.phi);
        SmearedIntersection.y = fRadius*TMath::Sin(SmearedIntersection.phi);
    }
    return SmearedIntersection;
}

void Detector::SetBranchAddress(TTree& gentree,TTree& rectree, unsigned countdet)
{
/*
 *  Function that sets the address of the branch of the output tree to the
 *  vectors containing the intersections
 *  -------------------------
 *  Parameters:
 *  gentree,rectree: TTree&
 *      Trees for data output
 * 
 *  countdet: unsigned
 *      running count of detectors
 * 
 */
    gentree.SetBranchAddress((string("GenHits_")+std::to_string(countdet)).c_str(), &fTrueHitPtr);
    if (fSmearing)
        rectree.SetBranchAddress((string("RecHits_")+std::to_string(countdet)).c_str(), &fRecoHitPtr);
    else    
        rectree.SetBranchAddress((string("RecHits_")+std::to_string(countdet)).c_str(), &fTrueHitPtr);

}

void Detector::Noise()
{
    MaterialBudget::fPoint noise;
    noise.isIntersection = true;
    static double  pi=TMath::Pi(), Twopi=2*pi;
    int noisehits=gRandom->Poisson(fMeanNoise);
    for (int i=0; i<noisehits;++i)
    {
        noise.z=(gRandom->Rndm()-0.5)*fLength;
        double phi=gRandom->Rndm()*Twopi;
        noise.phi=phi-pi/2;
        noise.x=fRadius*TMath::Cos(phi);
        noise.y=fRadius*TMath::Sin(phi);
        fRecoHit.push_back(noise);
    }
}