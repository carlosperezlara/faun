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
    fQA_FSIGNALS  = new TH2D( "MPC_FSIGNALS", ";Crystal;Sgn",      23,-0.5,22.5,300,-99.5,2399.5);
    output->Add( fQA_FSIGNALS );
    fQA_CMN0_BL = new TH2D( "MPC_CMN0_BL",";crystal;baseline",23,-0.5,22.5,300,1800.5,2100.5);
    output->Add( fQA_CMN0_BL );
    fQA_SIGNALS  = new TH2D( "MPC_SIGNALS", ";Crystal;Sgn",      23,-0.5,22.5,300,-99.5,2399.5);
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
  if(fDOQA)
    for(unsigned int i=0; i!=fCrystals.size(); ++i) {
      for(unsigned int j=0; j!=1024; ++j) {
	fQA_S[i]->Fill(j,fCrystals[i]->GetDataSlice(j));
	fQA_CMN0_S[i]->Fill(j,fCrystals[i]->GetDataSlice(j)-fCommonNoise0[j]);
      }
      fQA_SIGNALS->Fill(i,fCrystals[i]->Signal());
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

void FDetectorMPC::setup() {
  // estimate baseline for each crystal in this event;
  TH1D *bl = new TH1D("bl","bl",100,-0.5,99.5);
  TF1 *fbl = new TF1("fbl","[0]");
  for(int i=0; i!=23; ++i) {
    for(int j=0; j!=100; ++j)
      bl->SetBinContent( j+1, fCrystals[i]->GetDataSlice(j) );
    fbl->SetParameter(0,1950);
    bl->Fit( "fbl", "QLN" );
    fBaseline[i] = fbl->GetParameter(0);
    fQA_CMN0_BL->Fill(i,fBaseline[i]);
    bl->Reset();
  }
  delete bl;
  delete fbl;
  /*
  std::cout << "BASELINE:" << std::endl;
  for(int ii=0; ii!=23; ++ii)
    std::cout << fBaseline[ii] << "|";
  std::cout << std::endl;
  */
  //fast estimation of the signal for each crystal in this event
  float signal[23];
  for(int i=0; i!=23; ++i) {
    signal[i] = 0;
    for(int j=250; j!=300; ++j)
      signal[i] += (fCrystals[i]->GetDataSlice(j)-fBaseline[i])/50.;
    fQA_SIGNALS->Fill(i,signal[i]);
  }
  /*
  std::cout << "SIGNAL:" << std::endl;
  for(int ii=0; ii!=23; ++ii)
    std::cout << signal[ii] << "|";
  std::cout << std::endl;
  */
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
  std::vector<double> mylist;
  for(int j=0; j!=1024; ++j) {
    mylist.clear();
    float delta = 0;
    //delta += fCrystals[ fOrd[20] ]->GetDataSlice(j) - fBaseline[ fOrd[20] ];
    //fQA_CMN0_N->Fill( fOrd[20] );
    //fCommonNoise0.push_back( delta );
    int ncfcpn = 8;
    for(int fl=23-ncfcpn; fl!=23; ++fl) {
      delta += fCrystals[ fOrd[fl] ]->GetDataSlice(j) - fBaseline[ fOrd[fl] ];
      fQA_CMN0_N->Fill( fOrd[fl] );
    }
    fCommonNoise0.push_back( delta/ncfcpn );
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
  // computes baseline and sort
  setup();
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
