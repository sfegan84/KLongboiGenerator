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
#include <TLorentzVector.h>

using namespace std;

class PrintEvents{

 public:

  PrintEvents()=default;
  
  void Write(vector<TLorentzVector> * part4Vect);

  void WriteLund(vector<TLorentzVector> * part4Vect, vector<int> * pdg_ID, vector<TVector3> * vertex);

  void WriteHEPmc(vector<TLorentzVector> * part4Vect);

};
