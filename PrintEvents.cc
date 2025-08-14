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

