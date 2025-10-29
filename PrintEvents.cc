#include "PrintEvents.h"



void PrintEvents::Write(vector<TLorentzVector> * part4Vect) {

  //vector<TLorentzVector> part4Vect;
  
  //initial state particles (beam and target)
  cout<<"("<<part4Vect->at(0).M()<<","<<part4Vect->at(0).Px()<<","<<part4Vect->at(0).Py()<<","<<part4Vect->at(0).Pz()<<") ("
  <<part4Vect->at(1).M()<<","<<part4Vect->at(1).Px()<<","<<part4Vect->at(1).Py()<<","<<part4Vect->at(1).Pz()<<") -> ";

  //final state particles
  for (int fspartl=2;(int)fspartl<part4Vect->size(); fspartl++){
    cout<<"("<<part4Vect->at(fspartl).M()<<","<<part4Vect->at(fspartl).Px()<<","<<part4Vect->at(fspartl).Py()<<","<<part4Vect->at(fspartl).Pz()<<") ";
  }
  cout<<endl;

}


void PrintEvents::WriteLund(vector<TLorentzVector> * part4Vect, vector<int> * pdg_ID, vector<TVector3> * vertex) {
  //Implementation of functionality to write a Lund output file. Using CLAS12 conventions. UD == used designed, users can assign any meaning to them

  //Event Header
  // 1 number of particles
  // 2 Mass number of target (UD)
  // 3 Atomic number of target (UD)
  // 4 Target polarisation (UD)
  // 5 z component of first particle spin
  // 6 Beam type electron=11 photon=22 (UD)
  // 7 Beam energy GeV (UD)
  // 8 Interacted nucleon id 2212 or 2112 (UD)
  // 9 Process ID (UD)
  // 10 Event weight (UD)
  cout << "Lund format output selected" << endl;

  cout << "THIS IS A TEST: Lund file will look like this (redirect to a file stream)" << endl;

  cout << (int)part4Vect->size() << "\t"
       << right << "1"
       << right << "1"
       << right << "0"
       << right << "0"
       << right << "130"
       << right << part4Vect->at(0).E()
       << right << pdg_ID->at(1)
       << right << "0"
       << right << "1";

  //Particle Loop
  //1 index
  //2 Lifetime (ns) (UD)
  //3 type (only 1 is propagated in Geant4)
  //4 particle id
  //5 Index of parent (UD)
  //6 Index of first daughter (UD)
  //7 momentum x (GeV)
  //8 momentum y (GeV)
  //9 momentum z (GeV)
  //10 Energy of particle (GeV) (UD)
  //11 Mass of particle (GeV) (UD)
  //12 Vertex x (cm)
  //13 Vertex y (cm)
  //14 Vertex z (cm)

   for (int partl=0;(int)partl<part4Vect->size(); partl++){
        cout   << right << partl
	       << right << "0"
                << setw(6) << right << "1"
                << setw(14) << right << pdg_ID->at(partl)
                << setw(6) << right << "0"
                << setw(6) << right << "0"
	        << setw(14) << right << part4Vect->at(partl).Px()
                << setw(14) << right << part4Vect->at(partl).Py()
                << setw(14) << right << part4Vect->at(partl).Pz()
                << setw(14) << right << part4Vect->at(partl).E()
                << setw(14) << right << part4Vect->at(partl).M()
	       << setw(6) << right << vertex->at(partl).X()
	       << setw(6) << right << vertex->at(partl).Y()
	       << setw(6) << right << vertex->at(partl).Z()
                << endl;
   }

}


void PrintEvents::WriteHEPmc(vector<TLorentzVector> * part4Vect) {
  //Implementation of functionality to write a HEPMC output file
  cout << "HEPMC format output selected" << endl;

  cout << "THIS IS A TEST: HEPMC file will look like this (redirect to a file stream)" << endl;

}
