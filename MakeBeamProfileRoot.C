#include <TGraph.h>
#include <TH1D.h>
#include "TFile.h"
#include "TF1.h"

using namespace std;

void MakeBeamProfileRoot(){
    
    int whichroot=0;// 0 to read in KL beam profile data file and create a hist; 1 to create the photon beam profile from function
    
    string loc_NewFileName="BeamProfile.root";
    TFile *locFile = new TFile(loc_NewFileName.c_str(), "RECREATE");
    
    TF1 *fox3 = new TF1("fox3","1.0/x*(18.3484*(1+(1-x)*(1-x))-17.6817*(1-x)*(2.0/3.0))", 0, 12);
    
     TF1 *fo24p = new TF1("fo24p","-0.0511313*exp(-0.5*(x-1)*(x-1))+0.180109*exp(-0.5*(x-2)*(x-2))-0.0535854*exp(-0.5*(x-3)*(x-3))+0.217572*exp(-0.5*(x-4)*(x-4))+0.0026012*exp(-0.5*(x-5)*(x-5))+0.0792481*exp(-0.5*(x-6)*(x-6))+0.000730144*exp(-0.5*(x-7)*(x-7))+0.016732*exp(-0.5*(x-8)*(x-8))+0.000844976*exp(-0.5*(x-9)*(x-9))+0.00178876*exp(-0.5*(x-10)*(x-10))-0.000623221*exp(-0.5*(x-11)*(x-11))",0,12); //Function from Mikhail for the new KL beam profile, x is momentum.
    
    int numbbins=1;
    if (whichroot==0){
        numbbins=10000;
        TH1D *BeamProfile1=new TH1D("BeamProfile1","K_{L} beam profile; E_{K_{L}} (GeV); counts", numbbins,0, 10);
        TH1D *BeamProfile=new TH1D("BeamProfile","K_{L} beam profile; E_{K_{L}} (GeV); counts", numbbins,0, 10);
    }
    else if (whichroot==1) {
        numbbins=12000;
        TH1D *BeamProfile=new TH1D("BeamProfile","#gamma beam profile; E_{#gamma} (GeV); counts", numbbins,0, 12);
    }
    Double_t xval, yvalp, yvalEk;
    
    if (whichroot==0){
        const char nick[]="./kl_mom_prof.dat";
        TGraph *KLbeam=new TGraph(nick,"%lg %lg"," \t");
        
        for (int i=0; i<numbbins; i++){
            xval=10.0/numbbins*i+10.0/(2*numbbins);
            BeamProfile1->Fill(xval, KLbeam->Eval(xval));
            xval=BeamProfile->GetBinCenter(i);
            yvalp=fo24p->Eval(xval);
            yvalEk=TMath::Sqrt(yvalp*yvalp+0.497611*0.497611)-0.497611;
            BeamProfile->SetBinContent(i, yvalEk/0.048);
            cout<<xval<<" "<<yvalEk<<endl;

        }
    }
    else if (whichroot==1){
        for (int i=0; i<numbbins; i++){
            xval=BeamProfile->GetBinCenter(i);
            BeamProfile->SetBinContent(i, fox3->Eval(xval/12.6));
        }
    }
    BeamProfile1->Draw();
    BeamProfile->Draw("same");
    
    //fox3->Draw("same");
    locFile->Write();
    
}
