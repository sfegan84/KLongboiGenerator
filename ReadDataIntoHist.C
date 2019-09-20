
#include "Riostream.h"
void ReadDataIntoHist() {
    //  Read data from an ascii file and create a root file with an histogram and an ntuple.
    //   see a variant of this macro in basic2.C

    
    
    // read file $ROOTSYS/tutorials/tree/basic.dat
    // this file has 3 columns of float data
    TString dir = gSystem->UnixPathName(__FILE__);
    dir.ReplaceAll("ReadDataIntoHist.C","");
    dir.ReplaceAll("/./","/");
    ifstream in1, in2;
    in1.open(Form("%sSOL2--KLn-KpXim.dat",dir.Data()));
    in2.open(Form("%sSOL4--KLn-KpXim.dat",dir.Data()));

    Double_t w, t, x,y;
    Int_t nlines = 0;
    TFile *f = new TFile("XS_Py_KLn_KpXim.root","RECREATE");
    
    
    TGraph2D *g2x = new TGraph2D();
    TGraph2D *g2y = new TGraph2D();
    TGraph2D *g4x = new TGraph2D();
    TGraph2D *g4y = new TGraph2D();
    g2x->SetName("Sol2_XS");
    g4x->SetName("Sol4_XS");
    g2y->SetName("Sol2_Py");
    g4y->SetName("Sol4_Py");

    while (1) {
        in1 >> w>> t>> x >> y;
        if (!in1.good()) break;
        if (nlines < 5) printf("w=%8f, t=%8f, x=%8f, y=%8f\n",w,t,x,y);
        g2x->SetPoint(nlines, w, t, x);
        g2y->SetPoint(nlines, w, t, y/x);
        nlines++;
    }
    printf(" found %d points\n",nlines);
    in1.close();
    nlines = 0;
    while (1) {
        in2 >> w>> t>> x >> y;
        if (!in2.good()) break;
        if (nlines < 5) printf("w=%8f, t=%8f, x=%8f, y=%8f\n",w,t,x,y);
        g4x->SetPoint(nlines, w, t, x);
        g4y->SetPoint(nlines, w, t, y/x);
        nlines++;
    }
    printf(" found %d points\n",nlines);
    in2.close();

    g2x->Write();
    g4x->Write();
    g2y->Write();
    g4y->Write();
    
    f->Close();
}

