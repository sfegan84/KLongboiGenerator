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
       << right << "1" << "\t"
       << right << "1" << "\t"
       << right << "0" << "\t"
       << right << "0" << "\t"
       << right << "130" << "\t"
       << right << part4Vect->at(0).E() << "\t"
       << right << pdg_ID->at(1) << "\t"
       << right << "0" << "\t"
       << right << "1" << endl;

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
     cout   << right << partl << "\t"
	       << right << "0"
                << setw(6) << right << "1"
                << setw(14) << right << pdg_ID->at(partl)
                << setw(6) << right << "0"
                << setw(6) << right << "0"
	        << setw(14) << right << part4Vect->at(partl).Px()
                << setw(14) << right << part4Vect->at(partl).Py()
                << setw(14) << right << part4Vect->at(partl).Pz()
                << setw(14) << right << part4Vect->at(partl).E()
                << setw(14) << right << part4Vect->at(partl).M() << "\t"
	       << setw(6) << right << vertex->at(partl).X() << "\t"
	       << setw(6) << right << vertex->at(partl).Y() << "\t"
	       << setw(6) << right << vertex->at(partl).Z()
                << endl;
   }

}


void PrintEvents::WriteHEPmc(vector<TLorentzVector> * part4Vect, vector<int> * pdg_ID, vector<TVector3> * vertex) {
  //Implementation of functionality to write a HEPMC output file
  cout << "HEPMC format output selected" << endl;

  cout << "THIS IS A TEST: HEPMC file will look like this (redirect to a file stream)" << endl;



  cout << "HepMC::Version 3.02.03" << endl;
  cout << "HepMC::Asciiv3-START_EVENT_LISTING" << endl;

  cout << "E " << part4Vect->size() << " " << part4Vect->size() << endl;
  cout << "U GEV CM" << endl;

     for (int partl=0;(int)partl<2; partl++){
       cout << setw(2) << "P"
	    << setw(2) << right << partl+1
	    << setw(2) << right << "0"
	    << setw(6) << right << pdg_ID->at(partl)
	    << setw(12) << right << part4Vect->at(partl).Px()
	    << setw(12) << right << part4Vect->at(partl).Py()
	    << setw(12) << right << part4Vect->at(partl).Pz()
	    << setw(12) << right << part4Vect->at(partl).E()
	    << setw(12) << right << part4Vect->at(partl).M() << "\t";
	 if(partl ==0){
	   cout << setw(2) << right << "4"
	            << endl;
	     }
	 else{
	   cout << setw(2) << right << "1"
	            << endl;
	     }	   

     }
     for (int partl=2;(int)partl<part4Vect->size(); partl++){
       cout << setw(2) << "P"
	    << setw(2) << right << partl+1
	    << setw(2) << right << "0"
	    << setw(6) << right << pdg_ID->at(partl)
	    << setw(12) << right << part4Vect->at(partl).Px()
	    << setw(12) << right << part4Vect->at(partl).Py()
	    << setw(12) << right << part4Vect->at(partl).Pz()
	    << setw(12) << right << part4Vect->at(partl).E()
	    << setw(12) << right << part4Vect->at(partl).M() << "\t"
	    << setw(2) << right << "1"
	    << endl;

       cout << setw(2) << "V"
	    << setw(2) << right << partl+1
	    << setw(2) << right << "0\t"
	    << setw(6) << right << vertex->at(partl).X() << "\t"
	    << setw(6) << right << vertex->at(partl).Y() << "\t"
	    << setw(6) << right << vertex->at(partl).Z() << "\t"
	    << setw(6) << right << part4Vect->at(partl).T()
	    << endl;
   }
  

  cout << "HepMC::Asciiv3-END_EVENT_LISTING" << endl;

}
