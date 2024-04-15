#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stddef.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <signal.h>
#include <TROOT.h>
#include <TFile.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TF1.h>
#include <TGraph2D.h>
//#include <vector>
//#include <TLorentzVector.h>
#include <TVector3.h>
#include <TRandom3.h>
#include "JGenBeamEnergy.h"   //JLA  including new jlab routine ...
#include <TGenPhaseSpace.h>     //  ... replacing TGenPhaseSpace
#include "JGenFermiMomentum.h"
#include <TInterpreter.h>
//for HDDM format
#include "HDDM/hddm_s.hpp"
#include "particleType.h"

#define mass_protPDG 0.93827201
#define mass_neutPDG 0.939565378
#define mass_kaonplusPDG 0.493677
#define mass_kaonminusPDG 0.493677
#define mass_kaonzeroPDG 0.497614
#define mass_photonPDG 0.0
#define mass_cascPDG 1.31486
#define mass_cascminusPDG 1.32171
#define mass_lambdPDG 1.115683
#define mass_sigzeroPDG 1.192642
#define mass_sigplusPDG 1.18937
#define mass_sigminusPDG 1.19745
#define mass_piplusPDG 0.139570
#define mass_piminusPDG 0.139570
#define mass_pizeroPDG 0.1349766

#define G3ID_prot 14
#define G3ID_neut 13
#define G3ID_kaonplus 11
#define G3ID_kaonminus 12
#define G3ID_kaonlong 10
#define G3ID_kaonshort 16
#define G3ID_photon 1
#define G3ID_casc 22
#define G3ID_cascminus 23
#define G3ID_lambd 18
#define G3ID_sigzero 20
#define G3ID_sigplus 19
#define G3ID_sigminus 21
#define G3ID_piplus 8
#define G3ID_piminus 9
#define G3ID_pizero 7

#define PDGID_photon 22
#define PDGID_kaonlong 130
#define PDGID_kaonshort 310
#define PDGID_kaonplus 321
#define PDGID_kaonminus -321
#define PDGID_prot 2212
#define PDGID_neut 2112
#define PDGID_lambd 3122
#define PDGID_sigplus 3222
#define PDGID_sigzero 3212
#define PDGID_sigminus 3112
#define PDGID_casc 3322
#define PDGID_cascminus 3312
#define PDGID_piplus 211
#define PDGID_piminus -211
#define PDGID_pizero 111

//const double TARGET_Z_CENTER = 65.;   // GlueX standard
//const double TARGET_Z_LENGTH = 30.;   // GlueX standard
const double TARGET_Z_CENTER = 60.;  // cm
const double TARGET_Z_LENGTH = 40.;  // cm
const double TARGET_RADIUS = 3.;     // cm
const double KPT_Z_LENGTH = 40.;     // cm
const double DISTANCE_KPT_SPECTROMETER = 24.*100.;   // cm
const double SPEED_OF_LIGHT_CM_NS = 29.9792;   // cm/ns

using namespace std;
void PrintUsage (char *processName);
int main (int argc, char **argv){
    
    gROOT->ProcessLine(".L Loader.C+");
    
    extern char *optarg; //to be used by getopt to return optionargument
    char *ProgName = argv[0];
    char *ROOTFILE = NULL;
    char *REACTION = NULL;
    char *OBS = NULL;
    int WillBeRootOutput = 0;
    int PrintEvents=0;
    int c;
    int max = 1;
    bool SIMULATE_BEAM_TIMING = true;
    bool SIMULATE_KPT_THICKNESS = true;
    bool SIMULATE_TARGET_SIZE = true;
    
    JGenBeamEnergy *beamE=new JGenBeamEnergy(kaon, histo, 0.0, 8.0); //JLA  default plain distribution
    if (argc == 1){
        PrintUsage (ProgName);
        exit (0);
    }
    while ((c = getopt (argc, argv, "hE:F:M:R:P:S:t:c")) != -1){ //returns -1 if no more options are present
        switch (c){
            case 'h':
                PrintUsage (ProgName);
                break;
            case 'M':
                max = atoi (optarg);
                break;
            case 'F':
                WillBeRootOutput = 1;
                ROOTFILE = optarg;
                break;
            case 'R':
                REACTION = optarg;
                break;
            case 'P':
                PrintEvents = 1;
                break;
            case 'E':
                beamE = new JGenBeamEnergy(optarg);   //JLA added option -E to generate energy distribution
                break;
            case 'S':
                OBS = optarg;   //added option -S to generate polarised cross sections
                break;
            case 't':
                SIMULATE_BEAM_TIMING = false;
                break;
            case 'c':
                SIMULATE_KPT_THICKNESS = false;
                break;
            case 'g':
                SIMULATE_TARGET_SIZE = false;
                break;
            default:
                fprintf (stderr, "Unrecognized argument: [-%c]\n\n", c);
                PrintUsage (ProgName);
                exit (1);
                break;
        }
    }
    enum keywordReaction_t {kl1,kl2,kl3,kl4,kl5,kl6,kln1,kln2,kln3,kln4,kln5,kln6,kln7,kln8,g1,g2,g3,n1,n2,n3,n4,n5,n6}; //reaction keyword
    char* keywordReaction [] = {"kl1","kl2","kl3","kl4","kl5","kl6","kln1","kln2","kln3","kln4","kln5","kln6","kln7","kln8","g1","g2","g3","n1","n2","n3","n4","n5","n6", NULL }; //reaction keyword
    int   ReactionKey=0;
    while (REACTION && keywordReaction [ReactionKey] && strcasecmp (keywordReaction[ReactionKey], REACTION)) {
        ReactionKey++;
    }
    if (!REACTION || !keywordReaction[ReactionKey]){
        cout<<"Reaction not found. Generating default."<<endl;
        ReactionKey=0;
    }
    
    enum keywordObs_t {nsol,sol2,sol4}; //reaction keyword
    char* keywordObs [] = {"nsol","sol2","sol4", NULL }; //reaction keyword
    int   ObsKey=0;
    while (OBS && keywordObs [ObsKey] && strcasecmp (keywordObs[ObsKey], OBS)) {
        ObsKey++;
    }
    if (!OBS || !keywordObs[ObsKey]){
        cout<<"Solution not found. Generating phase space."<<endl;
        ObsKey=0;
    }
    
    
    //Variables for Tree
    int num_tracks;
    vector<int> pdg_ID, geant_ID, charge;
    vector<TLorentzVector> part4Vect;
    vector<TVector3> vertex;
    double cs, poly, costhetaK, costhetaPi, Wval;
    /////////
    bool cascadegen=false;
    
    double beamMass;
    double *FSmasses, *FSmasses1; //final state masses
    TLorentzVector W, target4vec, temp4vect, pizero4vec;
    TLorentzVector *temp4vectpoint=NULL;
    TVector3 tempVert, tempVert2, beta, yaxis, zaxis; //vertex position
    double weight;
    double lifetime;
    double piontheta;
    double gammaval;
    zaxis.SetXYZ(0,0,1);
    
    geant_ID.clear();
    charge.clear();
    part4Vect.clear();
    vertex.clear();
    pdg_ID.clear();
    
    Particle_t targetType = Proton;
    Particle_t beamType;
    int fermi=0;
    
    
    int usesol=0;
    TGraph2D *SolXS=NULL;
    TGraph2D *SolPY=NULL;
    Double_t w, t, x,y;
    Int_t nlines = 0;
    ifstream in1;

    double maxxs=0.0, maxpy=0.0;
    switch ((keywordObs_t) ObsKey) {
        case nsol:
            usesol=0;
            cout<<"Generating phase space"<<endl;
            break;
        case sol2:
            cout<<"Generating sol2"<<endl;
            usesol=1;
            in1.open("SOL2--KLn-KpXim.dat");
            SolXS = new TGraph2D();
            SolPY = new TGraph2D();
            SolXS->SetName("Sol2_XS");
            SolPY->SetName("Sol2_Py");
            while (1) {
                in1 >> w>> t>> x >> y;
                if (!in1.good()) break;
                if (nlines < 5) printf("w=%8f, t=%8f, x=%8f, y=%8f\n",w,t,x,y);
                SolXS->SetPoint(nlines, w, t, x);
                SolPY->SetPoint(nlines, w, t, y/x);
                nlines++;
            }
            printf(" found %d points\n",nlines);
            in1.close();
            maxxs=SolXS->GetMaximum();
            maxpy=SolPY->GetMaximum();
            maxxs=0.013;
            break;
        case sol4:
            cout<<"Generating sol4"<<endl;
            usesol=1;
            in1.open("SOL4--KLn-KpXim.dat");
            SolXS = new TGraph2D();
            SolPY = new TGraph2D();
            SolXS->SetName("Sol2_XS");
            SolPY->SetName("Sol2_Py");
            while (1) {
                in1 >> w>> t>> x >> y;
                if (!in1.good()) break;
                if (nlines < 5) printf("w=%8f, t=%8f, x=%8f, y=%8f\n",w,t,x,y);
                SolXS->SetPoint(nlines, w, t, x);
                SolPY->SetPoint(nlines, w, t, y/x);
                nlines++;
            }
            printf(" found %d points\n",nlines);
            in1.close();
            maxxs=SolXS->GetMaximum();
            maxpy=SolPY->GetMaximum();
            maxxs=0.045;
            break;
        default:
            fprintf (stderr, "Unrecognized argument: [-%c]\n\n. Rinning nosolution", c);
            usesol=0;
            break;
    }
    
    
    TFile *RootOut;
    TTree *mytree;
    string mystring;//=ROOTFILE;
    string mystring2;//=mystring.substr(0, mystring.length()-4);
    mystring2.append("hddm");
    
    
    mystring=ROOTFILE;
    mystring2=mystring.substr(0, mystring.length()-4);
    mystring2.append("hddm");
    
    if (WillBeRootOutput){
        cout<<"Saving events in file "<<mystring<<" and "<<mystring2<<endl;
    }
    char *Reactionstring=keywordReaction[ReactionKey];

    RootOut=new TFile(ROOTFILE, "recreate", "");    // This is ROOT output file
    mytree = new TTree ("mytree", "A TTree object");
    mytree->Branch("num_tracks",&num_tracks);
    mytree->Branch("part4Vect",&part4Vect);
    mytree->Branch("geant_ID",&geant_ID);
    mytree->Branch("pdg_ID",&pdg_ID);
    mytree->Branch("charge",&charge);
    mytree->Branch("vertex",&vertex);
    mytree->Branch("Reaction",Reactionstring,"Reactionstring/C", 1024);
    
    
    switch ((keywordReaction_t) ReactionKey) {
        case kl1:
            num_tracks=4;
            beamType = KLong;
            beamMass=mass_kaonzeroPDG;
            geant_ID.push_back(G3ID_kaonlong);
            geant_ID.push_back(G3ID_prot);
            geant_ID.push_back(G3ID_kaonplus);
            geant_ID.push_back(G3ID_neut);
            pdg_ID.push_back(PDGID_kaonlong);
            pdg_ID.push_back(PDGID_prot);
            pdg_ID.push_back(PDGID_kaonplus);
            pdg_ID.push_back(PDGID_neut);
            charge.push_back(0);
            charge.push_back(1);
            charge.push_back(1);
            charge.push_back(0);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_kaonplusPDG;
            FSmasses[1]=mass_neutPDG;
            break;
        case kl2:
            num_tracks=4;
            beamType = KLong;
            beamMass=mass_kaonzeroPDG;
            geant_ID.push_back(G3ID_kaonlong);
            geant_ID.push_back(G3ID_prot);
            geant_ID.push_back(G3ID_kaonshort);
            geant_ID.push_back(G3ID_prot);
            pdg_ID.push_back(PDGID_kaonlong);
            pdg_ID.push_back(PDGID_prot);
            pdg_ID.push_back(PDGID_kaonshort);
            pdg_ID.push_back(PDGID_prot);
            charge.push_back(0);
            charge.push_back(1);
            charge.push_back(0);
            charge.push_back(1);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_kaonzeroPDG;
            FSmasses[1]=mass_protPDG;
            break;
        case kl3:
            num_tracks=4;
            beamType = KLong;
            beamMass=mass_kaonzeroPDG;
            geant_ID.push_back(G3ID_kaonlong);
            geant_ID.push_back(G3ID_prot);
            geant_ID.push_back(G3ID_kaonplus);
            geant_ID.push_back(G3ID_casc);
            pdg_ID.push_back(PDGID_kaonlong);
            pdg_ID.push_back(PDGID_prot);
            pdg_ID.push_back(PDGID_kaonplus);
            pdg_ID.push_back(PDGID_casc);
            charge.push_back(0);
            charge.push_back(1);
            charge.push_back(1);
            charge.push_back(0);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_kaonplusPDG;
            FSmasses[1]=mass_cascPDG;
            break;
        case kl4:
            num_tracks=4;
            beamMass=mass_kaonzeroPDG;
            geant_ID.push_back(G3ID_kaonlong);
            geant_ID.push_back(G3ID_prot);
            geant_ID.push_back(G3ID_piplus);
            geant_ID.push_back(G3ID_lambd);
            pdg_ID.push_back(PDGID_kaonlong);
            pdg_ID.push_back(PDGID_prot);
            pdg_ID.push_back(PDGID_piplus);
            pdg_ID.push_back(PDGID_lambd);
            charge.push_back(0);
            charge.push_back(1);
            charge.push_back(1);
            charge.push_back(0);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_piplusPDG;
            FSmasses[1]=mass_lambdPDG;
            break;
        case kl5:
            num_tracks=4;
            beamMass=mass_kaonzeroPDG;
            geant_ID.push_back(G3ID_kaonlong);
            geant_ID.push_back(G3ID_prot);
            geant_ID.push_back(G3ID_pizero);
            geant_ID.push_back(G3ID_sigplus);
            pdg_ID.push_back(PDGID_kaonlong);
            pdg_ID.push_back(PDGID_prot);
            pdg_ID.push_back(PDGID_pizero);
            pdg_ID.push_back(PDGID_sigplus);
            charge.push_back(0);
            charge.push_back(1);
            charge.push_back(0);
            charge.push_back(1);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_pizeroPDG;
            FSmasses[1]=mass_sigplusPDG;
            break;
        case kl6:
            num_tracks=4;
            beamMass=mass_kaonzeroPDG;
            geant_ID.push_back(G3ID_kaonlong);
            geant_ID.push_back(G3ID_prot);
            geant_ID.push_back(G3ID_piplus);
            geant_ID.push_back(G3ID_sigzero);
            pdg_ID.push_back(PDGID_kaonlong);
            pdg_ID.push_back(PDGID_prot);
            pdg_ID.push_back(PDGID_piplus);
            pdg_ID.push_back(PDGID_sigzero);
            charge.push_back(0);
            charge.push_back(1);
            charge.push_back(1);
            charge.push_back(0);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_piplusPDG;
            FSmasses[1]=mass_sigzeroPDG;
            break;
        case kln1:
            num_tracks=4;
            beamType = KLong;
            beamMass=mass_kaonzeroPDG;
            geant_ID.push_back(G3ID_kaonlong);
            geant_ID.push_back(G3ID_neut);
            geant_ID.push_back(G3ID_kaonminus);
            geant_ID.push_back(G3ID_prot);
            pdg_ID.push_back(PDGID_kaonlong);
            pdg_ID.push_back(PDGID_neut);
            pdg_ID.push_back(PDGID_kaonminus);
            pdg_ID.push_back(PDGID_prot);
            charge.push_back(0);
            charge.push_back(0);
            charge.push_back(-1);
            charge.push_back(1);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_kaonminusPDG;
            FSmasses[1]=mass_protPDG;
            fermi=1;
            break;
        case kln2:
            num_tracks=4;
            beamType = KLong;
            beamMass=mass_kaonzeroPDG;
            geant_ID.push_back(G3ID_kaonlong);
            geant_ID.push_back(G3ID_neut);
            geant_ID.push_back(G3ID_kaonshort);
            geant_ID.push_back(G3ID_neut);
            pdg_ID.push_back(PDGID_kaonlong);
            pdg_ID.push_back(PDGID_neut);
            pdg_ID.push_back(PDGID_kaonshort);
            pdg_ID.push_back(PDGID_neut);
            charge.push_back(0);
            charge.push_back(0);
            charge.push_back(0);
            charge.push_back(0);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_kaonzeroPDG;
            FSmasses[1]=mass_neutPDG;
            fermi=1;
            break;
        case kln3:
            num_tracks=4;
            beamType = KLong;
            beamMass=mass_kaonzeroPDG;
            geant_ID.push_back(G3ID_kaonlong);
            geant_ID.push_back(G3ID_neut);
            geant_ID.push_back(G3ID_kaonplus);
            charge.push_back(0);
            charge.push_back(0);
            charge.push_back(1);
            pdg_ID.push_back(PDGID_kaonlong);
            pdg_ID.push_back(PDGID_neut);
            pdg_ID.push_back(PDGID_kaonplus);
            if (usesol==1){
                geant_ID.push_back(G3ID_piminus);
                geant_ID.push_back(G3ID_lambd);
                charge.push_back(-1);
                charge.push_back(0);
                pdg_ID.push_back(PDGID_pizero);
                pdg_ID.push_back(PDGID_lambd);
                mytree->Branch("cs",&cs);
                mytree->Branch("poly",&poly);
                mytree->Branch("costhetaK",&costhetaK);
                mytree->Branch("W",&Wval);
                mytree->Branch("costhetaPi",&costhetaPi);
            }
            else {
                geant_ID.push_back(G3ID_cascminus);
                charge.push_back(-1);
                pdg_ID.push_back(PDGID_cascminus);
            }
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_kaonplusPDG;
            FSmasses[1]=mass_cascminusPDG;
            fermi=1;
            if (usesol==1){
                FSmasses1=new double[2];
                FSmasses1[0]=mass_piminusPDG;
                FSmasses1[1]=mass_lambdPDG;
            }
            cascadegen=true;
            break;
        case kln4:
            num_tracks=4;
            beamMass=mass_kaonzeroPDG;
            geant_ID.push_back(G3ID_kaonlong);
            geant_ID.push_back(G3ID_neut);
            geant_ID.push_back(G3ID_pizero);
            geant_ID.push_back(G3ID_lambd);
            pdg_ID.push_back(PDGID_kaonlong);
            pdg_ID.push_back(PDGID_neut);
            pdg_ID.push_back(PDGID_pizero);
            pdg_ID.push_back(PDGID_lambd);
            charge.push_back(0);
            charge.push_back(0);
            charge.push_back(0);
            charge.push_back(0);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_pizeroPDG;
            FSmasses[1]=mass_lambdPDG;
            fermi=1;
            break;
        case kln5:
            num_tracks=4;
            beamMass=mass_kaonzeroPDG;
            geant_ID.push_back(G3ID_kaonlong);
            geant_ID.push_back(G3ID_neut);
            geant_ID.push_back(G3ID_pizero);
            geant_ID.push_back(G3ID_sigzero);
            pdg_ID.push_back(PDGID_kaonlong);
            pdg_ID.push_back(PDGID_neut);
            pdg_ID.push_back(PDGID_pizero);
            pdg_ID.push_back(PDGID_sigzero);
            charge.push_back(0);
            charge.push_back(0);
            charge.push_back(0);
            charge.push_back(0);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_pizeroPDG;
            FSmasses[1]=mass_sigzeroPDG;
            fermi=1;
            break;
        case kln6:
            num_tracks=4;
            beamMass=mass_kaonzeroPDG;
            geant_ID.push_back(G3ID_kaonlong);
            geant_ID.push_back(G3ID_neut);
            geant_ID.push_back(G3ID_piminus);
            geant_ID.push_back(G3ID_sigplus);
            pdg_ID.push_back(PDGID_kaonlong);
            pdg_ID.push_back(PDGID_neut);
            pdg_ID.push_back(PDGID_piminus);
            pdg_ID.push_back(PDGID_sigplus);
            charge.push_back(0);
            charge.push_back(0);
            charge.push_back(-1);
            charge.push_back(1);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_piminusPDG;
            FSmasses[1]=mass_sigplusPDG;
            fermi=1;
            break;
        case kln7:
            num_tracks=4;
            beamMass=mass_kaonzeroPDG;
            geant_ID.push_back(G3ID_kaonlong);
            geant_ID.push_back(G3ID_neut);
            geant_ID.push_back(G3ID_piplus);
            geant_ID.push_back(G3ID_sigminus);
            pdg_ID.push_back(PDGID_kaonlong);
            pdg_ID.push_back(PDGID_neut);
            pdg_ID.push_back(PDGID_piplus);
            pdg_ID.push_back(PDGID_sigminus);
            charge.push_back(0);
            charge.push_back(0);
            charge.push_back(1);
            charge.push_back(-1);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_piplusPDG;
            FSmasses[1]=mass_sigminusPDG;
            fermi=1;
            break;
        case kln8:
            num_tracks=4;
            beamType = KLong;
            beamMass=mass_kaonzeroPDG;
            geant_ID.push_back(G3ID_kaonlong);
            geant_ID.push_back(G3ID_neut);
            geant_ID.push_back(G3ID_kaonshort);
            geant_ID.push_back(G3ID_casc);
            pdg_ID.push_back(PDGID_kaonlong);
            pdg_ID.push_back(PDGID_neut);
            pdg_ID.push_back(PDGID_kaonshort);
            pdg_ID.push_back(PDGID_casc);
            charge.push_back(0);
            charge.push_back(0);
            charge.push_back(0);
            charge.push_back(0);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_kaonzeroPDG;
            FSmasses[1]=mass_cascPDG;
            fermi=1;
            break;
        case g1:
            num_tracks=4;
            beamType = Gamma;
            beamMass=mass_photonPDG;
            geant_ID.push_back(G3ID_photon);
            geant_ID.push_back(G3ID_prot);
            geant_ID.push_back(G3ID_kaonplus);
            geant_ID.push_back(G3ID_lambd);
            pdg_ID.push_back(PDGID_photon);
            pdg_ID.push_back(PDGID_prot);
            pdg_ID.push_back(PDGID_kaonplus);
            pdg_ID.push_back(PDGID_lambd);
            charge.push_back(0);
            charge.push_back(1);
            charge.push_back(1);
            charge.push_back(0);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_kaonplusPDG;
            FSmasses[1]=mass_lambdPDG;
            break;
        case g2:
            num_tracks=4;
            beamType = Gamma;
            beamMass=mass_photonPDG;
            geant_ID.push_back(G3ID_photon);
            geant_ID.push_back(G3ID_prot);
            geant_ID.push_back(G3ID_kaonplus);
            geant_ID.push_back(G3ID_sigzero);
            pdg_ID.push_back(PDGID_photon);
            pdg_ID.push_back(PDGID_prot);
            pdg_ID.push_back(PDGID_kaonplus);
            pdg_ID.push_back(PDGID_sigzero);
            charge.push_back(0);
            charge.push_back(1);
            charge.push_back(1);
            charge.push_back(0);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_kaonplusPDG;
            FSmasses[1]=mass_sigzeroPDG;
            break;
        case g3:
            num_tracks=4;
            beamType = Gamma;
            beamMass=mass_photonPDG;
            geant_ID.push_back(G3ID_photon);
            geant_ID.push_back(G3ID_prot);
            geant_ID.push_back(G3ID_kaonshort);
            geant_ID.push_back(G3ID_sigplus);
            pdg_ID.push_back(PDGID_photon);
            pdg_ID.push_back(PDGID_prot);
            pdg_ID.push_back(PDGID_kaonshort);
            pdg_ID.push_back(PDGID_sigplus);
            charge.push_back(0);
            charge.push_back(1);
            charge.push_back(0);
            charge.push_back(1);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_kaonzeroPDG;
            FSmasses[1]=mass_sigplusPDG;
            break;
        case n1:
            num_tracks=5;
            beamType = Neutron;
            beamMass=mass_neutPDG;
            geant_ID.push_back(G3ID_neut);
            geant_ID.push_back(G3ID_prot);
            geant_ID.push_back(G3ID_kaonplus);
            geant_ID.push_back(G3ID_lambd);
            geant_ID.push_back(G3ID_neut);
            pdg_ID.push_back(PDGID_neut);
            pdg_ID.push_back(PDGID_prot);
            pdg_ID.push_back(PDGID_kaonplus);
            pdg_ID.push_back(PDGID_lambd);
            pdg_ID.push_back(PDGID_neut);
            charge.push_back(0);
            charge.push_back(1);
            charge.push_back(1);
            charge.push_back(0);
            charge.push_back(0);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_kaonplusPDG;
            FSmasses[1]=mass_lambdPDG;
            FSmasses[2]=mass_neutPDG;
            break;
        case n2:
            num_tracks=5;
            beamType = Neutron;
            beamMass=mass_neutPDG;
            geant_ID.push_back(G3ID_neut);
            geant_ID.push_back(G3ID_prot);
            geant_ID.push_back(G3ID_kaonplus);
            geant_ID.push_back(G3ID_sigzero);
            geant_ID.push_back(G3ID_neut);
            pdg_ID.push_back(PDGID_neut);
            pdg_ID.push_back(PDGID_prot);
            pdg_ID.push_back(PDGID_kaonplus);
            pdg_ID.push_back(PDGID_sigzero);
            pdg_ID.push_back(PDGID_neut);
            charge.push_back(0);
            charge.push_back(1);
            charge.push_back(1);
            charge.push_back(0);
            charge.push_back(0);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_kaonplusPDG;
            FSmasses[1]=mass_sigzeroPDG;
            FSmasses[2]=mass_neutPDG;
            break;
        case n3:
            num_tracks=5;
            beamType = Neutron;
            beamMass=mass_neutPDG;
            geant_ID.push_back(G3ID_neut);
            geant_ID.push_back(G3ID_prot);
            geant_ID.push_back(G3ID_kaonshort);
            geant_ID.push_back(G3ID_sigplus);
            geant_ID.push_back(G3ID_neut);
            pdg_ID.push_back(PDGID_neut);
            pdg_ID.push_back(PDGID_prot);
            pdg_ID.push_back(PDGID_kaonshort);
            pdg_ID.push_back(PDGID_sigplus);
            pdg_ID.push_back(PDGID_neut);
            charge.push_back(0);
            charge.push_back(1);
            charge.push_back(0);
            charge.push_back(1);
            charge.push_back(0);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_kaonzeroPDG;
            FSmasses[1]=mass_sigplusPDG;
            FSmasses[2]=mass_neutPDG;
            break;
        case n4:
            num_tracks=5;
            beamType = Neutron;
            beamMass=mass_neutPDG;
            geant_ID.push_back(G3ID_neut);
            geant_ID.push_back(G3ID_prot);
            geant_ID.push_back(G3ID_kaonshort);
            geant_ID.push_back(G3ID_lambd);
            geant_ID.push_back(G3ID_prot);
            pdg_ID.push_back(PDGID_neut);
            pdg_ID.push_back(PDGID_prot);
            pdg_ID.push_back(PDGID_kaonshort);
            pdg_ID.push_back(PDGID_lambd);
            pdg_ID.push_back(PDGID_prot);
            charge.push_back(0);
            charge.push_back(1);
            charge.push_back(0);
            charge.push_back(0);
            charge.push_back(1);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_kaonzeroPDG;
            FSmasses[1]=mass_lambdPDG;
            FSmasses[2]=mass_protPDG;
            break;
        case n5:
            num_tracks=5;
            beamType = Neutron;
            beamMass=mass_neutPDG;
            geant_ID.push_back(G3ID_neut);
            geant_ID.push_back(G3ID_prot);
            geant_ID.push_back(G3ID_kaonshort);
            geant_ID.push_back(G3ID_sigzero);
            geant_ID.push_back(G3ID_prot);
            pdg_ID.push_back(PDGID_neut);
            pdg_ID.push_back(PDGID_prot);
            pdg_ID.push_back(PDGID_kaonshort);
            pdg_ID.push_back(PDGID_sigzero);
            pdg_ID.push_back(PDGID_prot);
            charge.push_back(0);
            charge.push_back(1);
            charge.push_back(0);
            charge.push_back(0);
            charge.push_back(1);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_kaonzeroPDG;
            FSmasses[1]=mass_sigzeroPDG;
            FSmasses[2]=mass_protPDG;
            break;
        case n6:
            num_tracks=5;
            beamType = Neutron;
            beamMass=mass_neutPDG;
            geant_ID.push_back(G3ID_neut);
            geant_ID.push_back(G3ID_prot);
            geant_ID.push_back(G3ID_neut);
            geant_ID.push_back(G3ID_neut);
            geant_ID.push_back(G3ID_piplus);
            pdg_ID.push_back(PDGID_neut);
            pdg_ID.push_back(PDGID_prot);
            pdg_ID.push_back(PDGID_neut);
            pdg_ID.push_back(PDGID_neut);
            pdg_ID.push_back(PDGID_piplus);
            charge.push_back(0);
            charge.push_back(1);
            charge.push_back(0);
            charge.push_back(0);
            charge.push_back(1);
            FSmasses=new double[num_tracks-2];
            FSmasses[0]=mass_neutPDG;
            FSmasses[1]=mass_neutPDG;
            FSmasses[2]=mass_piplusPDG;
            break;
    }
    
    beamE->Generate(); // Generate beam energy
    temp4vect=beamE->GetP4();
    if (fabs(temp4vect.M()-beamMass)>0.0001){
        cout<<"Selected Beam (option -E) and reaction (option -R) do not match"<<endl;
        return 1;
    }
    target4vec.SetXYZT(0,0,0,mass_protPDG);
    TGenPhaseSpace eventPS, eventPS1;
    TF1 *lifetimef1=NULL;
    if (cascadegen && usesol==1){
        lifetimef1=new TF1("lifetimef1","TMath::Exp(-x/([0]*[1]))",0, 5); //in ns tau= 2.9E-10
        lifetimef1->SetParameter(0, 0.1639);
    }
    double runtemp=0;
    
    std::ofstream *outfile = new ofstream(mystring2.c_str());
    
    if (WillBeRootOutput && !outfile->is_open()) {
        std::cerr << "Unable to open output hddm file \"" << mystring2
        << "\" for writing." << std::endl;
        exit(-3);
    }
    
    hddm_s::ostream *outstream = new hddm_s::ostream(*outfile);
    
    
    if (!ROOTFILE){
        cout<<"File not Specified. Printing generated events on screen"<<endl;
        PrintEvents=1;
    }
    
    int numberofloops=0;
    int Nevents = 0;
    TRandom3 randomNum;
    int runNumber=9000; //Run number for recosnrtuction
    //get max weight
    TLorentzVector deuteron, starget;
    deuteron.SetXYZT(0,0,0,1.87561);
    
    if (usesol==1 && !cascadegen){
        cout<<"No polarised cross-section for reaction chosen. Generating Phase space."<<endl;
    }
    
    cout << "Generating the events..." <<endl;
    for (Nevents = 0; Nevents <max ;){
        numberofloops=0;
        if (fermi==1){
            double ran1 = 0.5*randomNum.Rndm();
            double ran2 = 11.0*randomNum.Rndm();
            if (JGenFermiMomentum::Instance().Spectral(ran1)<ran2)
                continue;
            double fcos = (2*randomNum.Rndm())-1;
            double fphi = 2*3.141592653*randomNum.Rndm();
            double fpx = ran1*sqrt(1-fcos*fcos)*cos(fphi);
            double fpy = ran1*sqrt(1-fcos*fcos)*sin(fphi);
            double fpz = ran1*fcos;
            double fEp = sqrt(fpx*fpx+fpy*fpy+fpz*fpz+0.93827231*0.93827231);
            double fEn = sqrt(fpx*fpx+fpy*fpy+fpz*fpz+0.939565378*0.939565378);
            starget.SetXYZT (-fpx, -fpy, -fpz, fEp);
            target4vec=deuteron-starget;
        }
        
        part4Vect.clear();
        vertex.clear();
        tempVert.SetXYZ(0, 0, 65);
        if(SIMULATE_TARGET_SIZE) {
        	double target_r = TARGET_RADIUS * sqrt(randomNum.Uniform());
        	double target_phi = randomNum.Uniform(0.,TMath::TwoPi());
        	tempVert.SetXYZ(target_r*cos(target_phi),target_r*sin(target_phi),TARGET_Z_CENTER+randomNum.Uniform(-TARGET_Z_LENGTH/2.,TARGET_Z_LENGTH/2.));
        }
        if (Nevents % 100 == 0){
            fprintf (stderr, "%d\r", Nevents);
            fflush (stderr);
        }
        
        beamE->Generate(); // Generate beam energy
        W=beamE->GetP4()+target4vec;
        
        if (eventPS.SetDecay(W, num_tracks-2, FSmasses)){
            
            part4Vect.push_back(beamE->GetP4());
            part4Vect.push_back(target4vec);
            vertex.push_back(tempVert);
            vertex.push_back(tempVert);
            
            weight = eventPS.Generate ();
            
            //cout<<"Max weight: "<<eventPS.GetWtMax()<<endl;
            if ((num_tracks>4) && (randomNum.Uniform(0,eventPS.GetWtMax())>weight)) // to produce flat distributions for events with more than 2 FS particles.
                continue;
            
            if (usesol==1 && cascadegen){
                TLorentzVector cms, Kaoncms;
                double costheta, trueW;
                double evalxs=0.0;
                double evalpy=0.0;
                temp4vectpoint=eventPS.GetDecay(0);
                temp4vect = *temp4vectpoint;
                //Kaon
                part4Vect.push_back(temp4vect);
                vertex.push_back(tempVert);
                
                cms=W;
                trueW=W.M();
                Wval=trueW;
                Kaoncms=temp4vect;
                Kaoncms.Boost(-cms.BoostVector());
                costheta=TMath::Cos(Kaoncms.Theta());
                costhetaK=costheta;
                evalxs=SolXS->Interpolate(trueW,costheta);
                cs=evalxs;
                runtemp=randomNum.Uniform(0,maxxs);
                if (runtemp>evalxs){
                    continue;
                }
                W=*eventPS.GetDecay(1);
                if (lifetimef1 && eventPS1.SetDecay(W, 2, FSmasses1)){
                    weight = eventPS1.Generate();
                    gammaval=eventPS.GetDecay(1)->Gamma();
                    lifetimef1->SetParameter(1,gammaval);
                    lifetime=lifetimef1->GetRandom(0,5.0);
                    beta=eventPS.GetDecay(1)->BoostVector();
                    tempVert2=tempVert+29.9792*lifetime*beta;
                    pizero4vec=*eventPS1.GetDecay(0);
                    pizero4vec.Boost(-beta);
                    yaxis=eventPS.GetDecay(1)->Vect().Cross(zaxis);
                    piontheta=TMath::Cos(pizero4vec.Angle(yaxis));
                    runtemp=randomNum.Uniform(0,1.50);
                    evalpy=SolPY->Interpolate(trueW,costheta);
                    poly=evalpy;
                    if(runtemp>(1+evalpy*0.458*piontheta)){
                        continue;
                    }
                    for (int fspartl=0;fspartl<2; fspartl++){
                        temp4vectpoint=eventPS1.GetDecay(fspartl);
                        temp4vect = *temp4vectpoint;
                        part4Vect.push_back(temp4vect);
                        vertex.push_back(tempVert2);
                    }
                    if (fermi==1){
                        part4Vect.push_back(starget);
                        geant_ID.push_back(G3ID_prot);
                        pdg_ID.push_back(PDGID_prot);
                        charge.push_back(1);
                    }
                    
                    if(PrintEvents){
                        cout<<"("<<part4Vect.at(0).M()<<","<<part4Vect.at(0).Px()<<","<<part4Vect.at(0).Py()<<","<<part4Vect.at(0).Pz()<<") ("
                        <<part4Vect.at(1).M()<<","<<part4Vect.at(1).Px()<<","<<part4Vect.at(1).Py()<<","<<part4Vect.at(1).Pz()<<") -> ";
                        for (int fspartl=2;fspartl<part4Vect.size(); fspartl++){
                            cout<<"("<<part4Vect.at(fspartl).M()<<","<<part4Vect.at(fspartl).Px()<<","<<part4Vect.at(fspartl).Py()<<","<<part4Vect.at(fspartl).Pz()<<") ";
                        }
                        cout<<endl;
                    }
                    
                    // Start a new event
                    hddm_s::HDDM record;
                    hddm_s::PhysicsEventList pes = record.addPhysicsEvents();
                    pes().setRunNo(runNumber);
                    pes().setEventNo(Nevents); //event number
                    hddm_s::ReactionList rs = pes().addReactions();
                    hddm_s::TargetList ts = rs().addTargets();
                    ts().setType(targetType);
                    hddm_s::PropertiesList tpros = ts().addPropertiesList();
                    tpros().setCharge(ParticleCharge(targetType));
                    tpros().setMass(ParticleMass(targetType));
                    hddm_s::MomentumList tmoms = ts().addMomenta();
                    
                    
                    
					// for now, assume that t=0 happens when the expected time of the photon bunch reaches z=0
					// (the front face of the spectrometer)
					// currently assuming that the distance between the center of the KPT target and the front
					// of the spectrometer is 24 m, and the KPT has a length of 40 cm
					// so set the event time based on that, assuming the kaons are moving purely in the z-direction
					// note that standard GlueX units are cm and ns
					double event_time = - DISTANCE_KPT_SPECTROMETER / SPEED_OF_LIGHT_CM_NS;
					if(SIMULATE_BEAM_TIMING) {
						double z_production = 0.;
						if(SIMULATE_KPT_THICKNESS) 
							z_production = randomNum.Uniform(KPT_Z_LENGTH) - KPT_Z_LENGTH/2.;
						double beam_velocity = (part4Vect.at(0).P() * SPEED_OF_LIGHT_CM_NS ) / part4Vect.at(0).E();
						event_time += ( DISTANCE_KPT_SPECTROMETER + z_production + vertex[0].Z() ) / beam_velocity;  
					}
                    
                    tmoms().setPx(0);
                    tmoms().setPy(0);
                    tmoms().setPz(0);
                    tmoms().setE(ParticleMass(targetType));
                    hddm_s::BeamList bs = rs().addBeams();
                    bs().setType(beamType);
                    
                    hddm_s::PropertiesList bpros = bs().addPropertiesList();
                    bpros().setCharge(ParticleCharge(beamType));
                    bpros().setMass(ParticleMass(beamType));
                    hddm_s::MomentumList bmoms = bs().addMomenta();
                    bmoms().setPx(part4Vect.at(0).Px());
                    bmoms().setPy(part4Vect.at(0).Py());
                    bmoms().setPz(part4Vect.at(0).Pz());
                    bmoms().setE(part4Vect.at(0).E());
                    
                    hddm_s::VertexList vs = rs().addVertices();
                    hddm_s::OriginList os = vs().addOrigins();
                    hddm_s::ProductList ps = vs().addProducts(part4Vect.size()-2);
                    
                    os().setT(event_time);
                    os().setVx(vertex.at(0).X());
                    os().setVy(vertex.at(0).Y());
                    os().setVz(vertex.at(0).Z());
                    for (int i=2; i < part4Vect.size(); i++) {
                        ps(i-2).setType((Particle_t) geant_ID[i]);
                        ps(i-2).setPdgtype(pdg_ID[i]);
                        ps(i-2).setId(i-1);         // unique value for this particle within the event
                        ps(i-2).setParentid(0);     // All internally generated particles have no parent
                        ps(i-2).setMech(0);       //   maybe this should be set to something?
                        
                        hddm_s::MomentumList pmoms = ps(i-2).addMomenta();
                        pmoms().setPx(part4Vect.at(i).Px());
                        pmoms().setPy(part4Vect.at(i).Py());
                        pmoms().setPz(part4Vect.at(i).Pz());
                        pmoms().setE(part4Vect.at(i).E());
                    }
                    if (WillBeRootOutput){
                        mytree->Fill ();
                        *outstream << record;
                    }
                    Nevents++;
                }
            }
            else {
                for (int fspartl=0;fspartl<num_tracks-2; fspartl++){
                    temp4vectpoint=eventPS.GetDecay(fspartl);
                    temp4vect = *temp4vectpoint;
                    part4Vect.push_back(temp4vect);
                    vertex.push_back(tempVert);
                }
                
                if(PrintEvents){
                    cout<<"("<<part4Vect.at(0).M()<<","<<part4Vect.at(0).Px()<<","<<part4Vect.at(0).Py()<<","<<part4Vect.at(0).Pz()<<") ("
                    <<part4Vect.at(1).M()<<","<<part4Vect.at(1).Px()<<","<<part4Vect.at(1).Py()<<","<<part4Vect.at(1).Pz()<<") -> ";
                    for (int fspartl=2;fspartl<num_tracks; fspartl++){
                        cout<<"("<<part4Vect.at(fspartl).M()<<","<<part4Vect.at(fspartl).Px()<<","<<part4Vect.at(fspartl).Py()<<","<<part4Vect.at(fspartl).Pz()<<") ";
                    }
                    cout<<endl;
                }
                
                // Start a new event
                hddm_s::HDDM record;
                
                hddm_s::PhysicsEventList pes = record.addPhysicsEvents();
                
                pes().setRunNo(runNumber);
                pes().setEventNo(Nevents); //event number
                
                hddm_s::ReactionList rs = pes().addReactions();
                hddm_s::TargetList ts = rs().addTargets();
                ts().setType(targetType);
                hddm_s::PropertiesList tpros = ts().addPropertiesList();
                tpros().setCharge(ParticleCharge(targetType));
                tpros().setMass(ParticleMass(targetType));
                hddm_s::MomentumList tmoms = ts().addMomenta();
                
				// for now, assume that t=0 happens when the expected time of the photon bunch reaches z=0
				// (the front face of the spectrometer)
				// currently assuming that the distance between the center of the KPT target and the front
				// of the spectrometer is 24 m, and the KPT has a length of 40 cm
				// so set the event time based on that, assuming the kaons are moving purely in the z-direction
				// note that standard GlueX units are cm and ns
				double event_time = - DISTANCE_KPT_SPECTROMETER / SPEED_OF_LIGHT_CM_NS;
				if(SIMULATE_BEAM_TIMING) {
					double z_production = 0.;
					if(SIMULATE_KPT_THICKNESS) 
						z_production = randomNum.Uniform(KPT_Z_LENGTH) - KPT_Z_LENGTH/2.;
					double beam_velocity = (part4Vect.at(0).P() * SPEED_OF_LIGHT_CM_NS ) / part4Vect.at(0).E();
					event_time += ( DISTANCE_KPT_SPECTROMETER + z_production + vertex[0].Z() ) / beam_velocity;  
				}
				
                tmoms().setPx(0);
                tmoms().setPy(0);
                tmoms().setPz(0);
                tmoms().setE(ParticleMass(targetType));
                hddm_s::BeamList bs = rs().addBeams();
                bs().setType(beamType);
                
                hddm_s::PropertiesList bpros = bs().addPropertiesList();
                bpros().setCharge(ParticleCharge(beamType));
                bpros().setMass(ParticleMass(beamType));
                hddm_s::MomentumList bmoms = bs().addMomenta();
                bmoms().setPx(part4Vect.at(0).Px());
                bmoms().setPy(part4Vect.at(0).Py());
                bmoms().setPz(part4Vect.at(0).Pz());
                bmoms().setE(part4Vect.at(0).E());
                
                hddm_s::VertexList vs = rs().addVertices();
                hddm_s::OriginList os = vs().addOrigins();
                hddm_s::ProductList ps = vs().addProducts(num_tracks-2);
                
                os().setT(event_time);
                os().setVx(vertex.at(0).X());
                os().setVy(vertex.at(0).Y());
                os().setVz(vertex.at(0).Z());
                
                
                for (int i=2; i < num_tracks; i++) {
                    
                    ps(i-2).setType((Particle_t) geant_ID[i]);
                    ps(i-2).setPdgtype(pdg_ID[i]);
                    ps(i-2).setId(i-1);         // unique value for this particle within the event
                    ps(i-2).setParentid(0);     // All internally generated particles have no parent
                    ps(i-2).setMech(0);       //   maybe this should be set to something?
                    
                    hddm_s::MomentumList pmoms = ps(i-2).addMomenta();
                    
                    pmoms().setPx(part4Vect.at(i).Px());
                    pmoms().setPy(part4Vect.at(i).Py());
                    pmoms().setPz(part4Vect.at(i).Pz());
                    pmoms().setE(part4Vect.at(i).E());
                    
                }
                
                if (WillBeRootOutput){
                    mytree->Fill ();
                    *outstream << record;
                }
                Nevents++;
            }
        }
        
        numberofloops++;
        if (numberofloops>10000){
            cout<<"Reaction Threshold above energies selected. Change Beam energy"<<endl;
            return 1;
        }
    }
    if (WillBeRootOutput){
        RootOut->Write ();
        delete outfile;
        
    }
    cout << "Number of events processed: " << Nevents << endl;
    return 0;
}


void PrintUsage (char *processName){
    cout << "\nUsage: " << processName << " [options] -FoutputRootFile\n\n";
    cout << "\nExample for 1 million events saved in root file generated.root for reaction Klong p --> Ks p\n\t and kaon beam kinetic energies between 1 and 4 GeV sampled from histogram :\n\t ~>" <<processName<< "  -M1000000 -Fgenerated.root -Ekaon:histo:1.0:4.0 -Rkl2\n\n";
    cout << "\toptions:\n";
    cout << "\t-h\t\t\tThis information\n";
    cout << "\t-F <filename>\t\tdirect ROOT <filename>. See below for root format\n";
    cout << "\t-R <reaction code>\tReaction to Generate. See below for codes\n";
    cout << "\t-M[#]\t\t\tProcess # number of events\n";
    cout << "\t-P\t\t\tPrint generated events\n";
    cout << "\t-t\t\t\tDisable simulated beam timing\n";
    cout << "\t-t\t\t\tDisable simulated production target (KPT) distribution\n";
    cout << "\t-g\t\t\tDisable simulated cyrotarget (LH2/LD2) distribution\n";
    cout << "\t-E <expression>\t\tSee below Values of E is kinetic energy of beam. If no E is specified, Kaon beam is assumed\n\t\t\t\t and it is sampled from BeamProfile_kaons.root \n\n";
    cout << "\t-S <Solution>\t\tSelected solution for polarised cross section. If no S is specified, phase space is generated.\n\t\t\t\tNow only reaction kln3 is supported. \n\n";
    cout << "\t <reaction code> \n";
    cout << "\tkl1\t Klong p --> K+ n \n";
    cout << "\tkl2\t Klong p --> Ks p \n";
    cout << "\tkl3\t Klong p --> K+ Xi \n";
    cout << "\tkl4\t Klong p --> pi+ Lambda \n";
    cout << "\tkl5\t Klong p --> pi0 Sigma+ \n";
    cout << "\tkl6\t Klong p --> pi+ Sigma \n";
    cout << "\tkln1\t Klong n --> K- p \n";
    cout << "\tkln2\t Klong n --> Ks n \n";
    cout << "\tkln3\t Klong n --> K+ Xi- \n";
    cout << "\tkln4\t Klong n --> pi0 Lambda \n";
    cout << "\tkln5\t Klong n --> pi0 Sigma \n";
    cout << "\tkln6\t Klong n --> pi- Sigma+ \n";
    cout << "\tkln7\t Klong n --> pi+ Sigma- \n";
    cout << "\tkln8\t Klong n --> Ks Xi \n";
    cout << "\tg1\t g p --> K+ Lambda \n";
    cout << "\tg2\t g p --> K+ Sigma \n";
    cout << "\tg3\t g p --> Ks Sigma+ \n";
    cout << "\tn1\t n p --> K+ Lambda n \n";
    cout << "\tn2\t n p --> K+ Sigma n \n";
    cout << "\tn3\t n p --> Ks Sigma+ n\n";
    cout << "\tn4\t n p --> Ks Lambda p\n";
    cout << "\tn5\t n p --> Ks Sigma p\n";
    cout << "\tn6\t n p --> n n pi+\n\n";
    cout << "\t <Solution> \n";
    cout << "\tnsol\t No solution; phase space is generated\n";
    cout << "\tsol2\t Solution 2 generated\n";
    cout << "\tsol4\t Solution 4 generated\n";
    cout << "\t<expression> (no whitespace): \n";
    cout << "\t[particle]:[distribution]:[x]:[y] \tDistribution (mono, plain, histo) for particle (kaon, neutron, photon) between [x] and [y]\n";
    cout << "\t              or only [x] for monoenergetic beams. For histo, the beam profile is sampled from  BeamProfile_particle.root\n\n";
    cout << "\tSome examples for option -E: \n";
    cout << "\t[x]           \t\t\tMonoenergetic kaon beam E = [x] GeV\n";
    cout << "\t[x]:[y]       \t\t\tPlain kaon distribution from [x] to [y] GeV\n";
    cout << "\tmono:[x]      \t\t\tMonoenergetic kaon beam as above\n";
    cout << "\tplain:[x]:[y] \t\t\tPlain kaon distribution as above\n";
    cout << "\thisto \t\t\t\tDistribution according to BeamProfile_kaon.root file\n";
    cout << "\thisto:[x] \t\t\tDistribution according to BeamProfile_kaon.root file. [x] is ignored\n";
    cout << "\thisto:[x]:[y] \t\t\tDistribution according to BeamProfile_kaon.root between [x] and [y]\n";
    cout << "\t[particle] \t\t\tParticle (kaon, neutron, photon) beam distribution from BeamProfile_[particle].root\n";
    cout << "\t[particle]:[x] \t\t\tParticle (kaon, neutron, photon) monoenergetic beam E=[x] GeV\n";
    cout << "\t[particle]:[x]:[y] \t\tParticle (kaon, neutron, photon) plain distribution from [x] to [y] GeV\n";
    cout << "\t[particle]:mono:[x] \t\tParticle (kaon, neutron, photon) monoenergetic beam E=[x] GeV\n";
    cout << "\t[particle]:plain:[x]:[y] \tParticle (kaon, neutron, photon) plain distribution from [x] to [y] GeV\n";
    cout << "\t[particle]:histo \t\tDistribution according to BeamProfile_particle.root file\n";
    cout << "\t[particle]:histo:[x] \t\tDistribution according to BeamProfile_particle.root file. [x] is ignored\n";
    cout << "\t[particle]:histo:[x]:[y] \tDistribution according to BeamProfile_particle.root file between [x] and [y]\n\n";
    cout << "\tRoot file contains tree named mytree and the following:\n";
    cout << "\tnum_tracks\t int number of tracks (beam+target+generated particles)\n";
    cout << "\tpart4Vect\t vector<TLorentzVector> vector of 4vectors of tracks (beam=part4Vect[0]; target=part4Vect[1]; genpart=part4Vect[...])\n";
    cout << "\tgeant_ID\t vector<int> geant_ID (beam+target+generated particles)\n";
    cout << "\tpdg_ID\t\t vector<int> pdg_ID (beam+target+genpart)\n";
    cout << "\tcharge\t\t vector<int> charge (beam+target+genpart)\n";
    cout << "\tvertex\t\t vector<TVector3> vertex of each track (beam=vertex[0]; target=vertex[1]; genpart=vertex[...])\n";
    exit (0);
}
