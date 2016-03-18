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
      fQA_F100[i]  = new TH2D( Form("MPC_F100_%d",i), ";ev;adc",100,-0.5,99.5,1024,-0.5,1023.5);
      fQA_S[i]  = new TH2D( Form("MPC_CRYSTAL%d",i), ";ts;adc",1024,-0.5,1023.5,4106,-0.5,4105.5);
      fQA_CMN0_S[i]  = new TH2D( Form("MPC_CMN0_CRYSTAL%d",i), ";ts;adc-cpn",1024,-0.5,1023.5,4106,-0.5,4105.5);
      output->Add( fQA_F100[i] );
      output->Add( fQA_S[i] );
      output->Add( fQA_CMN0_S[i] );
    }
    fQA_BL = new TH2D( "MPC_BL",";crystal;baseline",23,-0.5,22.5,200,1849.5,2050.5);
    output->Add( fQA_BL );
    fQA_SIGNALS  = new TH2D( "MPC_SIGNALS", ";Crystal;Sgn",      23,-0.5,22.5,460,-100,2200);
    output->Add( fQA_SIGNALS );
    fQA_ENERGIES = new TH2D( "MPC_ENERGIES",";Crystal;Ene",      23,-0.5,22.5,300,-4.5,49.5);
    output->Add( fQA_ENERGIES );
    fQA_ENERGY   = new TH2D( "MPC_ENERGY",  ";crystal;Energy",   23,-0.5,22.5,300,-4.5,49.5);
    fQA_CENTROID = new TH2D( "MPC_CENTROID",";<x>;<y>",          15,-6.25,+6.25,15,-6.25,+6.25);
    fQA_CMN0      = new TH2D( "MPC_CMN0", ";slice;<CMN0>", 1024,-0.5,1023.5,400,-150,150);
    fQA_CMN0_N    = new TH1D( "MPC_CMN0_N",   ";crystalindex",23,-0.5,22.5);
    fQA_CMN1      = new TH1D( "MPC_CMN1",     ";CMN1",          200,-500,+500);
    fQA_CMN1_N    = new TH1D( "MPC_CMN1_N",   ";N_{vec}",        24,-0.5,23.5);
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
  static int proEvents=0;
  if(fDOQA)
    for(unsigned int i=0; i!=fCrystals.size(); ++i) {
      float weight=0;
      for(unsigned int j=0; j!=1024; ++j) {
	weight = fCrystals[i]->GetDataSlice(j);
	fQA_F100[i]->Fill(proEvents,j,weight);
	fQA_S[i]->Fill(j,fCrystals[i]->GetDataSlice(j));
	fQA_CMN0_S[i]->Fill(j,fCrystals[i]->GetDataSlice(j)-fCommonNoise0[j]);
      }
      fQA_SIGNALS->Fill(i,fCrystals[i]->Signal());
      fQA_BL->Fill(i,fCrystals[i]->Baseline());
      fQA_ENERGIES->Fill(i,fCrystals[i]->Energy());
    }
  if(fQA_ENERGY) {
    double sum=0.0;
    for(unsigned int i=0; i!=fCrystals.size(); ++i) {
      sum += fCrystals[ fOrd[i] ]->Energy();
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

  proEvents++;
}

int FDetectorMPC::getneig(int nei[23], int sorted) {
  //   00 01 02
  //03 04 05 06 07
  //08 09 10 11 12
  //13 14 15 16 17
  //18 19 20 21 22
  std::vector<int> list;
  // n0 n1 n2
  // n3 XX n4
  // n5 n6 n7
  for(int tt=0; tt!=sorted; ++tt) {
    int t = fOrd[tt];
    //n0
    if(t==5||t==6||t==7) list.push_back( t-5 );
    if(t==9||t==10||t==11||t==12||
       t==14||t==15||t==16||t==17||
       t==19||t==20||t==21||t==22) list.push_back( t-6 );
    //n1
    if(t==4||t==5||t==6) list.push_back( t-4 );
    if(t>7) list.push_back( t-5 );
    //n2
    if(t==3||t==4||t==5) list.push_back( t-3 );
    if(t==8||t==9||t==10||t==11||
       t==13||t==14||t==15||t==16||
       t==18||t==19||t==20||t==21) list.push_back( t-4 );
    //n3
    if(t==1||t==2||
       t==4||t==5||t==6||t==7||
       t==9||t==10||t==11||t==12||
       t==14||t==15||t==16||t==17||
       t==19||t==20||t==21||t==22)
      list.push_back( t-1 );
    //n4
    if(t==0||t==1||
       t==3||t==4||t==5||t==6||
       t==8||t==9||t==10||t==11||
       t==13||t==14||t==15||t==16||
       t==18||t==19||t==20||t==21)
      list.push_back( t+1 );
    //n5
    if(t==0||t==1||t==2) list.push_back( t+3 );
    if(t==4||t==5||t==6||t==7||
       t==9||t==10||t==11||t==12||
       t==14||t==15||t==16||t==17) list.push_back( t+4 );
    //n6
    if(t==0||t==1||t==2) list.push_back( t+4 );
    if(t>2&&t<18) list.push_back( t+5 );
    //n7
    if(t==0||t==1||t==2) list.push_back( t+5 );
    if(t==3||t==4||t==5||t==6||
       t==8||t==9||t==10||t==11||
       t==13||t==14||t==15||t==16) list.push_back( t+6 );
  }
  /*
  std::cout << list.size() << "::";
  for(int j=0; j!=list.size(); ++j)
    std::cout << list[j] << "|";
  std::cout << std::endl;
  */
  std::vector<int>::iterator it;
  sort(list.begin(),list.end());
  it = unique(list.begin(),list.end());
  list.resize( distance(list.begin(),it) );
  /*
  std::cout << list.size() << "::";
  for(int j=0; j!=list.size(); ++j)
    std::cout << list[j] << "|";
  std::cout << std::endl;
  */
  for(int i=0; i!=sorted; ++i)
    for(int j=0; j!=list.size(); ++j)
      if(fOrd[i]==list[j]) list[j] = 99;
  sort(list.begin(),list.end());
  it = unique(list.begin(),list.end());
  list.resize( distance(list.begin(),it) );
  /*
  std::cout << list.size() << "::";
  for(int j=0; j!=list.size(); ++j)
    std::cout << list[j] << "|";
  std::cout << std::endl;
  */
  int nn=0;
  for(int j=0; j!=list.size(); ++j) 
    if(list[j]!=99)
      nei[nn++] = list[j];
  return nn;
}

void FDetectorMPC::Sort() {
  float signal[23];
  for(int i=0; i!=23; ++i)
    signal[i] = fCrystals[i]->Signal();

  //sorting crystals by beam proximity
  int neig[23];
  for(int i=1; i!=23; ++i)
    if(signal[i]>signal[fOrd[0]])
      fOrd[0] = i;
  int nord=1, nneig;
  while( nord<23 ) {
    /*
    std::cout << "NUMBER OF SORTED:" << nord << std::endl;
    for(int ii=0; ii!=nord; ++ii)
      std::cout << fOrd[ii] << "|";
    std::cout << std::endl;
    */
    nneig = getneig(neig,nord);
    /*
    std::cout << " NUMBER OF NEIGHBOURS:" << nneig << std::endl;
    for(int ii=0; ii!=nneig; ++ii)
      std::cout << neig[ii] << "|";
    std::cout << std::endl;
    */
    if(nneig>0) {
      fOrd[nord] = neig[0];
      float maxSgn = signal[ fOrd[nord] ];
      for(int i=1; i!=nneig; ++i) {
	int pk = neig[i];
	if(signal[pk]>maxSgn) {
	  maxSgn = signal[pk];
	  fOrd[ nord ] = pk;
	}
      }
      ++nord;
    }
  }
}

void FDetectorMPC::EstimateCommonNoise0() {
  fCommonNoise0.clear();
  std::vector<float> mylist;
  for(int j=0; j!=fCrystals[0]->NumberOfSlices(); ++j) {
    mylist.clear();
    float delta = 0;
    int ncfcpn = 7;
    for(int fl=23-ncfcpn; fl!=23; ++fl)
      mylist.push_back( fCrystals[ fOrd[fl] ]->GetDataSlice(j) - fCrystals[ fOrd[fl] ]->Baseline() );
    fCommonNoise0.push_back( mylist[ ncfcpn/2 ] );
  }
  return;
}

void FDetectorMPC::EstimateCommonNoise1() {
  fCommonNoise1 = 0;
  float delta=0;
  int ncfcn=7;
  std::vector<float> mylist;
  for(int fl=23-ncfcn; fl!=23; ++fl)
    mylist.push_back( fCrystals[ fOrd[fl] ]->Signal() );
  fCommonNoise1 = mylist[ncfcn/2];
  return;
}

void FDetectorMPC::Reset() {
  fCorrupt = false;
  for(std::vector<FAPD*>::size_type i = 0; i!=fCrystals.size(); ++i) {
    fCrystals[i]->Reset();
  }
}

void FDetectorMPC::ComputeSignals() {
  for(std::vector<FAPD*>::size_type i = 0; i!=fCrystals.size(); ++i) {
    fCrystals[i]->ComputeBaseline( fCommonNoise0 );
    //std::cout << i;
    fCrystals[i]->ComputeSignal( fCommonNoise0, fCommonNoise1 );
  }
}

void FDetectorMPC::ReadEnergy() {
  // each crystal has a trace amplitude that goes from 0 to 4096
  // the baseline is at around 2000

  // this procedure finds outliers and expectorates event
  for(std::vector<FAPD*>::size_type i = 0; i!=fCrystals.size(); ++i) {
    float sumJ=0;
    for(unsigned int j=0; j!=1024; ++j)
      if( fCrystals[i]->GetDataSlice(j)<100 || fCrystals[i]->GetDataSlice(j)>4090 ) {
	fCorrupt = true;
	return;
      }
  }

  fCommonNoise0.clear();
  fCommonNoise1 = 0;
  ComputeSignals(); // estimates signal

  // there is noise from the readout common to each time slice
  // the noise is removed by averaging the signal of the eight
  // furthest crystals from the interaction region
  //std::cout << "******\n";
  const int numberOfIterations = 3;
  for(int ii=0; ii!=numberOfIterations; ++ii) {
    Sort(); // sort
    //for(std::vector<FAPD*>::size_type i = 0; i!=fCrystals.size(); ++i)
    //  printf("|%.1f",fCrystals[ fOrd[i] ]->Signal());
    //std::cout << std::endl;
    EstimateCommonNoise0(); // estimate CPN
    ComputeSignals();
    EstimateCommonNoise1(); // estimate CN
    ComputeSignals();
  }
  if(fCommonNoise0.size()==0)
    for(int i=0; i!=1024; ++i)
      fCommonNoise0.push_back(0);
  //for(std::vector<FAPD*>::size_type i = 0; i!=fCrystals.size(); ++i)
  //  printf("|%.1f",fCrystals[ fOrd[i] ]->Signal());
  //std::cout << std::endl;

  fEnergy = 0;
  for(std::vector<FAPD*>::size_type i = 0; i!=fCrystals.size(); ++i)
    fEnergy += fCrystals[i]->Energy();
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
    fCrystals[i]->SetRange(0,150,230,630);
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
