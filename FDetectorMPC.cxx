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
#include "TF1.h"

FDetectorMPC::FDetectorMPC() {
  fDOQA = true;
}

FDetectorMPC::~FDetectorMPC() {
}

void FDetectorMPC::Dump() {
  for(unsigned int i=0; i!=fCrystals.size(); ++i)
    std::cout << "MPC Detector | APD " << i << " | [" << fCrystals[i]->GetGain() << "]" << std::endl;
}

TList* FDetectorMPC::Init() {
  InputGains();
  TList *output = new TList();
  output->SetName( "MPC" );
  if(fDOQA) {
    for(int i=0; i!=23; ++i) {
      fQA_S[i]  = new TH2D( Form("MPC_CRYSTAL%d",i), ";ts;adc",1024,-0.5,1023.5,4106,-0.5,4105.5);
      fQA_CMN0_S[i]  = new TH2D( Form("MPC_CMN0_CRYSTAL%d",i), ";ts;adc-cpn",1024,-0.5,1023.5,4106,-0.5,4105.5);
      output->Add( fQA_S[i] );
      output->Add( fQA_CMN0_S[i] );
    }
    fQA_CMN0_BL = new TH2D( "MPC_CMN0_BL",";crystal;baseline",23,-0.5,22.5,200,1900.5,2100.5);
    output->Add( fQA_CMN0_BL );
    fQA_SIGNALS  = new TH2D( "MPC_SIGNALS", ";Crystal;Sgn",      23,-0.5,22.5,300,-99.5,2399.5);
    fQA_ENERGIES = new TH2D( "MPC_ENERGIES",";Crystal;Ene",      23,-0.5,22.5,300,-4.5,49.5);
    fQA_ENERGY   = new TH2D( "MPC_ENERGY",  ";crystal;Energy",   23,-0.5,22.5,300,-4.5,49.5);
    fQA_CENTROID = new TH2D( "MPC_CENTROID",";<x>;<y>",          15,-6.25,+6.25,15,-6.25,+6.25);
    fQA_CMN0      = new TProfile( "MPC_CMN0", ";slice;<CMN0>", 1024,-0.5,1023.5,"s");
    fQA_CMN0_N    = new TH2D( "MPC_CMN0_N",   ";slice;N_{vec}",1024,-0.5,1023.5,25,-0.5,24.5);
    fQA_CMN1      = new TH1D( "MPC_CMN1",     ";CMN1",          200,-500,+500);
    fQA_CMN1_N    = new TH1D( "MPC_CMN1_N",   ";N_{vec}",        24,-0.5,23.5);
    output->Add( fQA_SIGNALS );
    output->Add( fQA_ENERGIES );
    output->Add( fQA_ENERGY );
    output->Add( fQA_CENTROID );
    output->Add( fQA_CMN0 );
    output->Add( fQA_CMN0_N );
    output->Add( fQA_CMN1 );
    output->Add( fQA_CMN1_N );
  }
  return output;
}

void FDetectorMPC::DoQA() {
  if(fDOQA)
    for(unsigned int i=0; i!=fCrystals.size(); ++i)
      for(unsigned int j=0; j!=1024; ++j) {
	fQA_S[i]->Fill(j,fCrystals[i]->GetDataSlice(j));
	fQA_CMN0_S[i]->Fill(j,fCrystals[i]->GetDataSlice(j)-fCommonNoise0[j]);
      }
  if(fQA_SIGNALS)
    for(unsigned int i=0; i!=fCrystals.size(); ++i)
      fQA_SIGNALS->Fill(i,fCrystals[i]->Signal());
  if(fQA_ENERGIES)
    for(unsigned int i=0; i!=fCrystals.size(); ++i)
      fQA_ENERGIES->Fill(i,fCrystals[i]->Energy());
  if(fQA_ENERGY) {
    double sum=0.0;
    for(unsigned int i=0; i!=fCrystals.size(); ++i) {
      sum += fCrystals[i]->Energy();
      fQA_ENERGY->Fill(i,sum);
    }
  }
  if(fQA_CENTROID)
      fQA_CENTROID->Fill(CentroidX(),CentroidY());
  if(fQA_CMN0)
    for(unsigned int i=0; i!=fCommonNoise0.size(); ++i)
      fQA_CMN0->Fill(i,fCommonNoise0[i]);
  if(fQA_CMN1)
    fQA_CMN1->Fill(fCommonNoise1);
}

void FDetectorMPC::EstimateCommonNoise0() {
  fCommonNoise0.clear();
  // estimate baseline for each crystal in this event;
  float baseline[23];
  TH1D *bl = new TH1D("bl","bl",100,-0.5,99.5);
  TF1 *fbl = new TF1("fbl","[0]");
  for(int i=0; i!=23; ++i) {
    for(int j=0; j!=100; ++j)
      bl->SetBinContent( i+1, j+1, fCrystals[i]->GetDataSlice(j) );
    fbl->SetParameter(0,1950);
    bl->Fit( "fbl", "QWW" );
    baseline[i] = fbl->GetParameter(0);
    fQA_CMN0_BL->Fill(i,baseline[i]);
    bl->Reset();
  }
  delete bl;
  delete fbl;

  int cut[23] = { 2020, 2035, 2088, 2031, 2035, 2035, 2057, 2028, 1999, 2023,
		  2026, 2055, 2040, 2017, 2031, 2026, 2074, 2077, 2027, 2145,
		  2068, 2093, 2078 };
  std::vector<double> mylist;
  for(int j=0; j!=1024; ++j) {
    mylist.clear();
    for(std::vector<FAPD*>::size_type i = 0; i!=fCrystals.size(); ++i) {
      double test = fCrystals[i]->GetDataSlice(j) - baseline[i];
      mylist.push_back( test );
    }
    double noise = 0;
    unsigned int size = mylist.size();
    if(fQA_CMN0_N)
      fQA_CMN0_N->Fill(j,size);
    if(size>0) {
      sort(mylist.begin(),mylist.end());
      noise = mylist[3];//23-4x4
    }
    fCommonNoise0.push_back( noise );
  }
  return;
}

void FDetectorMPC::EstimateCommonNoise1() {
  fCommonNoise1 = 0;
  std::vector<double> mylist;
  for(std::vector<FAPD*>::size_type i = 0; i!=fCrystals.size(); ++i) {
    double test = fCrystals[i]->Signal();
    if(test<50) mylist.push_back( test );
  }
  unsigned int size = mylist.size();
  if(fQA_CMN1_N)
    fQA_CMN1_N->Fill(size);
  if(size>0) {
    sort(mylist.begin(),mylist.end());
    fCommonNoise1 = mylist[size/2];
  }
  return;
}

void FDetectorMPC::Reset() {
  fCorrupt = false;
  for(std::vector<FAPD*>::size_type i = 0; i!=fCrystals.size(); ++i) {
    fCrystals[i]->Reset();
  }
}

void FDetectorMPC::ReadEnergy() {
  // each crystal has a trace amplitude that goes from 0 to 4096
  // the baseline is at around 2000
  // however there are events where the signal is at 0
  // this procedure removes those
  for(std::vector<FAPD*>::size_type i = 0; i!=fCrystals.size(); ++i) {
    float sumJ=0;
    for(unsigned int j=0; j!=1024; ++j)
      sumJ +=fCrystals[i]->GetDataSlice(j) / 1024.0;
    if( sumJ < 100 ) {
      fCorrupt = true;
      return;
    }
  }
  //
  // there is noise from the readout common to each time slice
  // from the variation 
  EstimateCommonNoise0();
  for(std::vector<FAPD*>::size_type i = 0; i!=fCrystals.size(); ++i) {
    fCrystals[i]->SetRange(0,150,230,675);
    double chk = fCrystals[i]->ComputeSignal(fCommonNoise0);
    //if(chk>15000) fCorrupt = true; // kill early signals
  }
  EstimateCommonNoise1();
  fCommonNoise1 = 0;
  fEnergy = 0;
  for(std::vector<FAPD*>::size_type i = 0; i!=fCrystals.size(); ++i) {
    fCrystals[i]->ComputeSignal(fCommonNoise0,fCommonNoise1);
    fEnergy += fCrystals[i]->Energy();
  }
}

double FDetectorMPC::CentroidX() {
  std::vector<FAPD*>::size_type size = fCrystals.size();
  if(size<1) return -999.;
  double x = 0;
  double en = 0;
  for(std::vector<FAPD*>::size_type i = 0; i!=size; ++i)
    if(fCrystals[i]->Energy()>0.5) {
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

void FDetectorMPC::InputGains() {
  // Geometry
  double x[23] = {-2.5, 0.0,+2.5, -5.0,-2.5, 0.0,+2.5,+5.0, -5.0,-2.5, 0.0,+2.5,+5.0, -5.0,-2.5, 0.0,+2.5,+5.0, -5.0,-2.5, 0.0,+2.5,+5.0};
  double y[23] = {+5.0,+5.0,+5.0, +2.5,+2.5,+2.5,+2.5,+2.5,  0.0, 0.0, 0.0, 0.0, 0.0, -2.5,-2.5,-2.5,-2.5,-2.5, -5.0,-5.0,-5.0,-5.0,-5.0};
  for(int i=0; i!=23; ++i) {
    fCrystalX.push_back(x[i]);
    fCrystalY.push_back(y[i]);
    fCrystals.push_back(new FAPD());
  }

  // Gains
  std::ifstream calib("MPC.gains");
  double tmp;
  for(std::vector<FAPD*>::size_type i=0; i!=fCrystals.size(); ++i)
    fCrystals[i]->SetGain( 1 );
  for(std::vector<FAPD*>::size_type i=0; i!=fCrystals.size(); ++i) {
    calib >> tmp;
    if(!calib.good()) break;
    fCrystals[i]->SetGain( tmp );
  }
  calib.close();
  return;
}
