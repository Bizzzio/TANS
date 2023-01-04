#include "Run.h"
#include "Reconstruction.h"

ClassImp(Run)

Run::Run(TString cfgFileName) : 
fConfigFile(YAML::LoadFile(cfgFileName.Data())),
fOutFileName(fConfigFile["OutFileName"].as<std::string>()),
fTreeGen("fTreeGen","fTreeGen"),
fTreeRec("fTreeRec","fTreeRec"),
fNEvents(fConfigFile["NEvents"].as<unsigned>()),

//Multiplicity settings
fMultType(fConfigFile["Multiplicity"]["MultType"].as<std::string>()),
fConstMult(fConfigFile["Multiplicity"]["MultConst"].IsNull() ? 0 : fConfigFile["Multiplicity"]["MultConst"].as<unsigned>()),
fMultFile(fConfigFile["Multiplicity"]["MultFile"].as<std::string>()),
fMultHistoName(fConfigFile["Multiplicity"]["MultHistoName"].as<std::string>()),
fMultRange(fConfigFile["Multiplicity"]["MultRange"].IsNull() ?  std::vector<unsigned>{} : fConfigFile["Multiplicity"]["MultRange"].as<std::vector<unsigned> >()),
//Angular distribution settings
fDistrType(fConfigFile["AngularDistr"]["DistrType"].as<std::string>()),
fDistrFile(fConfigFile["AngularDistr"]["DistrFile"].as<std::string>()),
fDistrHisto(fConfigFile["AngularDistr"]["DistrHisto"].as<std::string>()),
fDistrConst(fConfigFile["AngularDistr"]["DistrConst"].as<std::vector<double> >()),

//Vertex dispersion settings
fSigmaX(fConfigFile["SigmaX"].as<double>()),
fSigmaY(fConfigFile["SigmaY"].as<double>()),
fSigmaZ(fConfigFile["SigmaZ"].as<double>()),

//Material budgets settings
fIsDetector(fConfigFile["IsDetector"].as<std::vector<bool> >()),
fRadii(fConfigFile["Radii"].as<std::vector<double> >()),
fThickness(fConfigFile["Thickness"].as<std::vector<double> >()),
fLength(fConfigFile["Lenght"].as<std::vector<double> >()),
fMaterial(fConfigFile["Material"].as<std::vector<string> >()),
fDoMultScattering(fConfigFile["DoMultScattering"].as<bool>()),
fDoSmearing(fConfigFile["DoSmearing"].as<bool>()),
fDoNoise(fConfigFile["DoNoise"].as<bool>()),
fMeanNoise(fConfigFile["MeanNoise"].as<int>()),

fVerbose(fConfigFile["Verbose"].as<bool>())
{
    TStopwatch SimulationTime, WritingTime, ReconstructionTime;
    SimulationTime.Start();

    std::vector<Event::fVertMult>* config = nullptr;        //used to define Branch type
    fTreeGen.Branch("Config",&config);                      //Branch for multipliciy and vertex position

    CreateDetectors();
    if (fMultType.find("kConst") != std::string::npos)
        fMultFunction=&Run::GetConstMult;
    else if (fMultType.find("kUniform") != std::string::npos)
        fMultFunction=&Run::GetUniformMult;
    else if (fMultType.find("kCustom") != std::string::npos)
    {
        TFile* infile = TFile::Open(fMultFile.c_str());
        fMultHisto = (TH1D*)infile->Get(fMultHistoName.c_str());
        fMultHisto -> SetDirectory(0);
        infile->Close();
        fMultFunction=&Run::GetCustomMult;
    }           //TODO: add default case

    if (fDistrType.find("kConst") != std::string::npos)
        RunConstDistr();
    else if (fDistrType.find("kUniform") != std::string::npos)
        RunUniformDistr();
    else if (fDistrType.find("kCustom") != std::string::npos)      
        RunCustomDistr();

    SimulationTime.Stop();
    
    WritingTime.Start();
    TFile hfile(fOutFileName.c_str(),"recreate");
    fTreeGen.Write();
    fTreeRec.Write();
    hfile.Close();
    WritingTime.Stop();
    ReconstructionTime.Start();
    Reconstruction analysis;
    ReconstructionTime.Stop();
    cout<<"Simulation time:"<<endl;
    SimulationTime.Print("u");
    cout<<"Writing time:"<<endl;
    WritingTime.Print("u");
    cout<<"Reconstruction time:"<<endl;
    ReconstructionTime.Print("u");
}

inline unsigned Run::GetCustomMult()
{
    unsigned mult;
    if (fMultRange.size()>0)
    {
        do {
            mult = fMultHisto->GetRandom();
        }while (mult < fMultRange[0] || mult > fMultRange[1]);
    }
    else
        mult = fMultHisto->GetRandom();
    return mult;
}

void Run::RunConstDistr()
{
/*
 *  Function that runs fNEvents with constant multiplicity
 */
    if (fDistrConst.size() != 3)
    {
        cout<<"\033[93mWarning: DistrConst is not a vector (dim!=3), setting it to (1,0,0)\033[0m \n";
        fDistrConst={1,0,0};
    }
    for (unsigned i=1; i<=fNEvents; i++)
    {
        if (i%10000==0 && fVerbose)
            cout<<"Processing event "<<i<<endl;
        Event(fDetectors,(this->*fMultFunction)(), gRandom->Gaus(0,fSigmaX),gRandom->Gaus(0,fSigmaY),gRandom->Gaus(0,fSigmaZ),fDistrConst,fTreeGen,fTreeRec);     
    }
}

void Run::RunUniformDistr()          //TODO: check case "MultRange" is not defined
{
/*
 *  Function that runs fNEvents with multiplicity distributed evenly in the range (fMultRange[0],fMultRange[1])
 */
    for (unsigned i=1; i<=fNEvents; i++)
    {
        if (i%10000==0 && fVerbose)
            cout<<"Processing event "<<i<<endl;
        Event(fDetectors,(this->*fMultFunction)(), gRandom->Gaus(0,fSigmaX),gRandom->Gaus(0,fSigmaY),gRandom->Gaus(0,fSigmaZ),fTreeGen,fTreeRec);     
    }
}

void Run::RunCustomDistr()
{
/*
 *  Function that runs fNEvents with multiplicity distributed according to an histogram
 *  If fMultRange is define, multiplicity is restricted to that range
 */
    TFile* infile = TFile::Open(fDistrFile.c_str());
    TH1D* histopseudorap = (TH1D*)infile->Get(fDistrHisto.c_str());
    histopseudorap->SetDirectory(0);
    infile->Close();
    double mintheta=2*TMath::ATan(exp(-(histopseudorap->GetBinLowEdge(histopseudorap->GetNbinsX()) + histopseudorap->GetBinWidth(histopseudorap->GetNbinsX()))));
    cout<<"Mintetha "<<mintheta<<endl;
    double maxtheta=2*TMath::ATan(exp(-(histopseudorap->GetBinLowEdge(1))));  
    cout<<"Maxtetha "<<maxtheta<<endl;
    cout<<2*TMath::ATan(exp(-(histopseudorap->GetBinLowEdge(1))))<<endl;
    TH1D* histotheta = new TH1D("Theta","Theta",histopseudorap->GetNbinsX(),mintheta,maxtheta);
    for (int i=1; i<=histopseudorap->GetNbinsX();i++)
    {
        double theta=histotheta->GetBinCenter(i);
        double eta=-TMath::Log(TMath::Tan(theta/2));
        double etadensity = histopseudorap->GetBinContent(histopseudorap->FindBin(eta));
        double derivative=TMath::Abs(1/TMath::Cos(theta/2)/TMath::Sin(theta/2));
        histotheta->Fill(histotheta->GetBinCenter(i),etadensity*derivative);
    }
    for (unsigned i=0; i<fNEvents; i++)
    {
        if (i%10000==0 && i!=0 && fVerbose)
            cout<<"Processing event "<<i<<endl;
        Event(fDetectors,(this->*fMultFunction)(), gRandom->Gaus(0,fSigmaX),gRandom->Gaus(0,fSigmaY),gRandom->Gaus(0,fSigmaZ),histotheta,fTreeGen,fTreeRec);     
    }
    delete histotheta;
}


void Run::CreateDetectors()     //TODO: add default cases
{
/*
 *  Function that generates the array of detectors based on the ConfigFile info
 */
    std::vector<MaterialBudget::fPoint>*    genhits     = nullptr;          //used to define Branch type
    std::vector<MaterialBudget::fPoint>*    rechits     = nullptr;          //used to define Branch type

    unsigned DetCount = 0;                                                  //counts number of detectors
    for (vector<bool>::size_type i=0; i<fIsDetector.size(); ++i)
    {
        if (fIsDetector[i])
        {
            if (fVerbose)
                cout<<"Creating Detector"<<endl;
            fDetectors.push_back(new Detector(fThickness[i],fRadii[i],fLength[i],fMaterial[i],fDoMultScattering,fDoSmearing,fDoNoise,fMeanNoise));
            ++DetCount;
            fTreeGen.Branch((string("GenHits_")+std::to_string(DetCount)).c_str(),&genhits);
            fTreeRec.Branch((string("RecHits_")+std::to_string(DetCount)).c_str(),&rechits);
            fDetectors[i]->SetBranchAddress(fTreeGen,fTreeRec,DetCount);    //SetsAddress of corresponding branches
        }
        else
        {
            if (fVerbose)
                cout<<"Creating MaterialBudget"<<endl;
            fDetectors.push_back(new MaterialBudget(fThickness[i],fRadii[i],fLength[i],fMaterial[i],fDoMultScattering));
        }
    }

    //Sorts detectors/materialbudgets based on the radii
    std::sort(fDetectors.begin(),fDetectors.end(),[] (MaterialBudget const *const a, MaterialBudget const *const b) { return a->GetRadius() < b->GetRadius(); });
}
