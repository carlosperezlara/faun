#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>

#include "FDetectorMPC.h"
#include "FAPD.h"

#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

FDetectorMPC::FDetectorMPC() {
}

FDetectorMPC::~FDetectorMPC() {
}

void FDetectorMPC::Dump() {
  std::cout << "fCrystals[" << fCrystals.size() << "]" << std::endl;
}

TList* FDetectorMPC::Init() {
  TList *output = new TList();
  output->SetName( "MPC" );
  int x[23] = { -1, 0, +1, -2, -1, 0, +1, +2, -2, -1, 0, +1, +2, -2, -1, 0, +1, +2, -2, -1, 0, +1, +2 };
  int y[23] = { +2, +2, +2, +1, +1, +1, +1, +1, 0, 0, 0, 0, 0, -1, -1, -1, -1, -1, +2, +2, +2, +2, +2 };
  for(int i=0; i!=23; ++i) {
    fCrystalX.push_back(x[i]);
    fCrystalY.push_back(y[i]);
    fCrystals.push_back(new FAPD());
  }
  InputGains("MPC.calib");
  if(1) {
    fQA_APDS     = new TH2D( "MPC_APDS",    ";Crystal;slice",  23,-0.5,22.5,1024,-0.5,1023.5);
    fQA_SIGNALS  = new TH2D( "MPC_SIGNALS", ";Crystal;Sgn",    23,-0.5,22.5,300,-99.5,2399.5);
    fQA_ENERGIES = new TH2D( "MPC_ENERGIES",";Crystal;Ene",    23,-0.5,22.5,300,-4.5,49.5);
    fQA_ENERGY   = new TH1D( "MPC_ENERGY",  ";Energy",        300,-4.5,49.5);
    fQA_CENTROID = new TH2D( "MPC_CENTROID",";<x>;<y>",       69,-3,+3,69,-3,+3);
    fQA_CMN      = new TProfile( "MPC_CMN", ";slice;<CMN>",  1024,-0.5,1023.5);
    output->Add( fQA_APDS );
    output->Add( fQA_SIGNALS );
    output->Add( fQA_ENERGIES );
    output->Add( fQA_ENERGY );
    output->Add( fQA_CENTROID );
    output->Add( fQA_CMN );
  }
  return output;
}

void FDetectorMPC::DoQA() {
  if(fQA_APDS)
    for(unsigned int i=0; i!=fCrystals.size(); ++i)
      for(unsigned int j=0; j!=1024; ++j)
	fQA_APDS->Fill(i,j,fCrystals[i]->GetDataSlice(j));
  if(fQA_SIGNALS)
    for(unsigned int i=0; i!=fCrystals.size(); ++i)
      fQA_SIGNALS->Fill(i,fCrystals[i]->Signal());
  if(fQA_ENERGIES)
    for(unsigned int i=0; i!=fCrystals.size(); ++i)
      fQA_ENERGIES->Fill(i,fCrystals[i]->Energy());
  if(fQA_ENERGY)
      fQA_ENERGY->Fill(fEnergy);
  if(fQA_CENTROID)
      fQA_CENTROID->Fill(CentroidX(),CentroidY());
  if(fQA_CMN)
    for(unsigned int i=0; i!=fCommonNoise.size(); ++i)
      fQA_CMN->Fill(i,fCommonNoise[i]);
}

void FDetectorMPC::EstimateCommonNoise() {
  fCommonNoise.clear();
  std::vector<double> mylist;
  for(unsigned int j=0; j!=1024; ++j) {
    mylist.clear();
    for(std::vector<FAPD*>::size_type i = 0; i!=fCrystals.size(); ++i) {
      double test = fCrystals[i]->GetDataSlice(j);
      if(test<2015) mylist.push_back( test );
    }
    double noise = 0;
    if(mylist.size()>0) {
      sort(mylist.begin(),mylist.end());
      noise = mylist[mylist.size()/2];
    }
    fCommonNoise.push_back( noise );
  }
  return;
}

void FDetectorMPC::Reset() {
  for(std::vector<FAPD*>::size_type i = 0; i!=fCrystals.size(); ++i) {
    fCrystals[i]->Reset();
  }
}

void FDetectorMPC::ReadEnergy() {
  EstimateCommonNoise();
  fEnergy = 0;
  for(std::vector<FAPD*>::size_type i = 0; i!=fCrystals.size(); ++i) {
    fCrystals[i]->SetRange(0,75,200,675);
    fCrystals[i]->ComputeSignal(fCommonNoise);
    fEnergy += fCrystals[i]->Energy();
  }
}

double FDetectorMPC::CentroidX() {
  std::vector<FAPD*>::size_type size = fCrystals.size();
  if(size<1) return -999.;
  double x = 0;
  double en = 0;
  for(std::vector<FAPD*>::size_type i = 0; i!=size; ++i)
    if(fCrystals[i]->Energy()>1.0) {
      x += fCrystals[i]->Energy()*fCrystalX[i];
      en += fCrystals[i]->Energy();
    }
  if(en<1.0) return -999.;
  return x/en;
}

double FDetectorMPC::CentroidY() {
  std::vector<FAPD*>::size_type size = fCrystals.size();
  if(size<1) return -999.;
  double x = 0;
  double en = 0;
  for(std::vector<FAPD*>::size_type i = 0; i!=size; ++i)
    if(fCrystals[i]->Energy()>1.0) {
      x += fCrystals[i]->Energy()*fCrystalY[i];
      en += fCrystals[i]->Energy();
    }
  if(en<1.0) return -999.;
  return x/en;
}

void FDetectorMPC::InputGains(char *file) {
  std::ifstream calib;
  calib.open(file);
  double tmp;
  for(std::vector<FAPD*>::size_type i=0; i!=fCrystals.size(); ++i)
    fCrystals[i]->SetGain( 1 );
  for(std::vector<FAPD*>::size_type i=0; i!=fCrystals.size(); ++i) {
    calib >> tmp;
    fCrystals[i]->SetGain( tmp );
    printf("Gain read for DetectorMPC || crystal %lu <=> %f\n",i,tmp);
    if(!calib.good()) break;
  }
  calib.close();
  return;
}
