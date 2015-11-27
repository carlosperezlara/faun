#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>

#include "FDetectorEX.h"
#include "FMP.h"

#include "TList.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TH3D.h"

FDetectorEX::FDetectorEX() {
  fDOQA = true;
}

FDetectorEX::~FDetectorEX() {
}

double FDetectorEX::Energy(int lyr) {
  if(lyr<0||lyr>8) return fEnergy;
  double en = 0;
  for(int i=128*(lyr); i!=128*(lyr+1); ++i)
    if(fStateMap[i]==0)
      en += fMinipads[i]->Energy();
  return en;
}

void FDetectorEX::Fill(std::vector<short> packet, short cid[2], short clk[2]) {
  if(packet.size()!=2048) {
    fCorrupt = true;
    return;
  }
  fIsPresent = true;

  for(int i=0; i!=1024; ++i)
    fMinipads[i]->Fill(packet[fLowMap[i]],packet[fHighMap[i]]);
  fCIDFor = cid[0];
  fCIDBck = cid[1];
  fClkFor = clk[0];
  fClkBck = clk[1];
}

TList* FDetectorEX::Init(){ //FDetector has a funtion named Init which returns TList
  InputGains();
  TList *output = new TList();
  output->SetName( "EX" );
  if(fDOQA) {
    fQA_CellID = new TH2D("EX_CellID","EX_CellID",64,-0.5,63.5,64,-0.5,63.5);
    output->Add( fQA_CellID );
    fQA_Clock = new TH2D("EX_Clock","EX_Clock",8,-0.5,7.5,8,-0.5,7.5);
    output->Add( fQA_Clock );

    fQA_RawMap = new TH2D("EX_RawMap","EX_RawMap;idx;adc",2048,-0.5,2047.5,256,-0.5,255.5);
    output->Add( fQA_RawMap );
    fQA_PedMap = new TH2D("EX_PedMap","EX_PedMap;idx;adc",2048,-0.5,2047.5,20,-10.5,+9.5);
    output->Add( fQA_PedMap );
    fQA_RawPedMap = new TH2D("EX_RawPedMap","EX_RawPedMap;idx;adc",2048,-0.5,2047.5,260,-9.5,250.5);
    output->Add( fQA_RawPedMap );
    fQA_MIPMap = new TH2D("EX_MIPMap","EX_MIPMap;idx;adc",2048,-0.5,2047.5,260,-9.5,250.5);
    output->Add( fQA_MIPMap );
    fQA_EneMap = new TH2D("EX_EneMap","EX_EneMap;idx;GeV",2048,-0.5,2047.5,260,0.0,0.005);
    output->Add( fQA_EneMap );
    fQA_HL = new TH3D("EX_HL","EX_HL;adc;adc;idx",260,-9.5,250.5,260,-9.5,250.5,1024,-0.5,1023.5);
    output->Add( fQA_HL );
    fQA_Layer0 = new TH3D("EX_Layer0","EX_Layer0;x;y;lyr",32,-0.5,31.5,4,-0.5,3.5,4,-0.5,3.5);
    output->Add( fQA_Layer0 );
    fQA_Layer1 = new TH3D("EX_Layer1","EX_Layer1;x;y;lyr",4,-0.5,3.5,32,-0.5,31.5,4,-0.5,3.5);
    output->Add( fQA_Layer1 );

    fQA_CommonNoise = new TProfile("EX_CommonNoise","EX_CommonNoise;idx;<CMN>",2048,-0.5,2047.5,"s");
    fQA_CommonNoise->SetMarkerStyle(20);
    output->Add( fQA_CommonNoise );
    fQA_NCN = new TH2D("EX_NCommonNoise","EX_NCommonNoise;idx rw;no pads",64,-0.5,63.5,33,-0.5,32.5);
    output->Add( fQA_NCN );

    fQA_EneLyr = new TH2D("EX_EneLyr","EX_EneLyr;GeV",8,-0.5,7.5,1000,-0.001,0.999);
    output->Add( fQA_EneLyr );
    fQA_Energy = new TH1D("EX_Energy","EX_Energy;GeV",1000,-0.001,0.299);
    output->Add( fQA_Energy );
    fQA_Centroid = new TH3D("EX_Centroid","EX_Centroid;x;y;lyr",96,-3.2,+3.2,96,-3.2,+3.2,8,-0.5,7.5);
    output->Add( fQA_Centroid );
  }
  return output;
}

void FDetectorEX::DoQA(){
  FADC *tmp;
  if(fQA_CellID) {
    fQA_CellID->Fill(fCIDFor,fCIDBck);
  }
  if(fQA_Clock) {
    fQA_Clock->Fill(fClkFor,fClkBck);
  }

  if(fQA_RawMap) {
    for(std::vector<FMP*>::size_type i=0; i!=fMinipads.size(); ++i) {
      tmp = fMinipads[i]->GetHigh();
      fQA_RawMap->Fill( i, tmp->Get() );
      tmp = fMinipads[i]->GetLow();
      fQA_RawMap->Fill( i+1024, tmp->Get() );
    }
  }
  if(fQA_PedMap) {
    for(std::vector<FMP*>::size_type i=0; i!=fMinipads.size(); ++i) {
      tmp = fMinipads[i]->GetHigh();
      fQA_PedMap->Fill( i, tmp->Get()-tmp->GetPMean() );
      tmp = fMinipads[i]->GetLow();
      fQA_PedMap->Fill( i+1024, tmp->Get()-tmp->GetPMean() );
    }
  }
  if(fQA_RawPedMap) {
    for(int i=0; i!=1024; ++i) {
      if(fStateMap[i]) continue;
      fQA_RawPedMap->Fill( i, fMinipads[i]->GetHigh()->RawPed() );
      fQA_RawPedMap->Fill( i+1024, fMinipads[i]->GetLow()->RawPed() );
    }
  }
  if(fQA_MIPMap) {
    for(int i=0; i!=1024; ++i) {
      if(fStateMap[i]) continue;
      fQA_MIPMap->Fill( i, fMinipads[i]->GetHigh()->RawPed()-fCommonNoise0[i] );
      fQA_MIPMap->Fill( i+1024, fMinipads[i]->GetLow()->RawPed()-fCommonNoise1[i] );
      fQA_EneMap->Fill( i, fMinipads[i]->GetHigh()->Energy() );
      fQA_EneMap->Fill( i+1024, fMinipads[i]->GetLow()->Energy() );
      //printf( "%d %f %f\n", i, fMinipads[i]->GetHigh()->Energy(), fMinipads[i]->GetLow()->Energy() );
      fQA_HL->Fill( fMinipads[i]->GetHigh()->RawPed()-fCommonNoise0[i], fMinipads[i]->GetLow()->RawPed()-fCommonNoise1[i], i );
    }
  }
  if(fQA_Layer0) {
    for(std::vector<FMP*>::size_type i=0; i!=fMinipads.size(); ++i) {
      if(fStateMap[i]) continue;
      tmp = fMinipads[i]->GetHigh();
      double w = tmp->RawPed();
      if( w < 5*tmp->GetPSigma() ) continue;
      int lyr = i/128;
      int mm = (i%128)/32;
      int mM = (i%128)%32;
      if(lyr%2) fQA_Layer1->Fill( mm, mM, (lyr-1)/2, w );
      else fQA_Layer0->Fill( mM, mm, lyr/2, w );
    }
  }

  if(fQA_EneLyr) {
    for(int i=0; i!=8; ++i)
      fQA_EneLyr->Fill(i,Energy(i));
  }
  if(fQA_Energy) {
    fQA_Energy->Fill(fEnergy);
  }
  if(fQA_Centroid) {
    for(int i=0; i!=8; ++i)
      fQA_Centroid->Fill( CentroidX(i), CentroidY(i), i );
  }
}


void FDetectorEX::Dump() {
  for(std::vector<FMP*>::size_type i=0; i!=fMinipads.size(); ++i) {
    std::cout << "EX Detector | " << i << " | [" << fHighMap[i] << ";" << fLowMap[i] << "] | [";
    std::cout << fMinipads[i]->GetHigh()->GetPMean() << ";" << fMinipads[i]->GetLow()->GetPMean() << "] | [";
    std::cout << fMinipads[i]->GetHigh()->GetPSigma() << ";" << fMinipads[i]->GetLow()->GetPSigma() << "] | [";
    std::cout << fMinipads[i]->GetHigh()->GetGain() << ";" << fMinipads[i]->GetLow()->GetGain() << "]." << std::endl;
  }
}

void FDetectorEX::ReadEnergy(){
  if(fCIDFor>46) fCorrupt = true;
  if(fCIDBck>46) fCorrupt = true;
  EstimateCommonNoise();
  fEnergy = 0;
  for(std::vector<FMP*>::size_type i=0; i!=fMinipads.size(); ++i)
    if(fStateMap[i]==0)
      fEnergy += fMinipads[i]->Energy();
}

void FDetectorEX::EstimateCommonNoise() {
  fCommonNoise0.clear();
  fCommonNoise1.clear();
  std::vector<double> cmh;
  std::vector<double> cml;
  FADC *adc;
  double tmp;
  for( int lyr=0; lyr!=8; ++lyr ) {
    for( int row=0; row!=4; ++row ) {
      cmh.clear();
      cml.clear();
      for( int cha=0; cha!=32; ++cha ) {
	int idx = 128*lyr + row*32 +cha;
	if(fStateMap[idx]) continue;
	adc = fMinipads[idx]->GetHigh();
	tmp = adc->Get() - adc->GetPMean();
	if( tmp < 5*adc->GetPSigma() ) cmh.push_back( tmp );
	adc = fMinipads[idx]->GetLow();
	tmp = adc->Get() - adc->GetPMean();
	if( tmp < 5*adc->GetPSigma() ) cml.push_back( tmp );
      }
      if(cmh.size()==0) cmh.push_back(0);
      if(cml.size()==0) cml.push_back(0);
      fQA_NCN->Fill(lyr*4+row,cmh.size());
      fQA_NCN->Fill(32+lyr*4+row,cml.size());
      sort(cmh.begin(),cmh.end());
      sort(cml.begin(),cml.end());
      for(int cha=row*32+lyr*128; cha!=(row+1)*32+lyr*128; ++cha)
	fCommonNoise0.push_back( cmh[cmh.size()/2] );
      for(int cha=row*32+lyr*128; cha!=(row+1)*32+lyr*128; ++cha)
	fCommonNoise1.push_back( cml[cml.size()/2] );
    }
  }

  for(int i=0; i!=1024; ++i) {
    fQA_CommonNoise->Fill( i, fCommonNoise0[i] );
    fQA_CommonNoise->Fill( i+1024, fCommonNoise1[i] );
    fMinipads[i]->GetHigh()->SetCMN( fCommonNoise0[i] );
    fMinipads[i]->GetLow()->SetCMN( fCommonNoise1[i] );
  }
}

void FDetectorEX::Reset(){
  fCorrupt = false;
  fIsPresent = false;
  for(int i=0; i!=1024; ++i)
    fMinipads[i]->Reset();
}

void FDetectorEX::InputGains(){
  // Geometry
  for(int idx=0; idx!=1024; ++idx) {
    if( (idx/128)%2==0 ) {
      fMinipadX.push_back( ((idx%128)%32)*0.21 -3.255 );
      fMinipadY.push_back( ((idx%128)/32)*1.52 -3.040 );
    } else {
      fMinipadX.push_back( ((idx%128)/32)*1.52 -3.040 );
      fMinipadY.push_back( ((idx%128)%32)*0.21 -3.255 );
    }
    fMinipadZ.push_back( (idx/128)*1.0 );
    fMinipads.push_back( new FMP() );
  }

  // Map
  fLowMap.clear();
  fHighMap.clear();
  std::ifstream lomap("EX.lmap");
  std::ifstream himap("EX.hmap");
  int a, b;
  for(std::vector<FMP*>::size_type i=0; i!=fMinipads.size(); ++i) {
    himap >> a;
    lomap >> b;
    if(!lomap.good()||!himap.good()) break;
    fHighMap.push_back( a );
    fLowMap.push_back( b );
  }
  lomap.close();
  himap.close();

  double tmp1, tmp2, tmp3, tmp4;
  for(std::vector<FMP*>::size_type i=0; i!=fMinipads.size(); ++i) {
    fMinipads[i]->SetGains( 1, 1 );
    fMinipads[i]->SetPedestals( 0, 0, 0, 0 );
  }

  // State map
  int scan=0;
  std::ifstream stmap("EX.statemap");
  for(std::vector<FMP*>::size_type i=0; i!=fMinipads.size(); ++i) {
    stmap >> a >> b;
    if(!stmap.good()) break;
    for(;scan!=a;++scan) fStateMap.push_back(0);
    fStateMap.push_back( b );
  }
  stmap.close();
  for(;scan!=1024;++scan) fStateMap.push_back(0);

  // Pedestals
  std::ifstream hiped("EX.hped");
  std::ifstream loped("EX.lped");
  for(std::vector<FMP*>::size_type i=0; i!=fMinipads.size(); ++i) {
    hiped >> tmp1 >> tmp2;
    loped >> tmp3 >> tmp4;
    if(!hiped.good()||!loped.good()) break;
    fMinipads[i]->SetPedestals( tmp1, tmp2, tmp3, tmp4 );
  }
  lomap.close();
  himap.close();

  // Gains
  std::ifstream hgain("EX.hgains");
  std::ifstream lgain("EX.lgains");
  for(std::vector<FMP*>::size_type i=0; i!=fMinipads.size(); ++i) {
    hgain >> tmp1;
    lgain >> tmp2;
    if(!hgain.good()||!lgain.good()) break;
    fMinipads[i]->SetGains( tmp1, tmp2 );
  }
  hgain.close();
  lgain.close();
}

double FDetectorEX::CentroidX(int lyr){
  double sumX = 0;
  double sumE = 0;
  for(int i=lyr*128; i!=(lyr+1)*128; ++i) {
    if(fStateMap[i]) continue;
    sumE += fMinipads[i]->Energy();
    sumX += fMinipadX[i]*fMinipads[i]->Energy();
  }
  if(sumE<0.1) return -999.;
  return sumX/sumE;
}

double FDetectorEX::CentroidY(int lyr){
  double sumY = 0;
  double sumE = 0;
  for(int i=lyr*128; i!=(lyr+1)*128; ++i) {
    if(fStateMap[i]) continue;
    sumE += fMinipads[i]->Energy();
    sumY += fMinipadY[i]*fMinipads[i]->Energy();
  }
  if(sumE<0.1) return -999.;
  return sumY/sumE;
}



