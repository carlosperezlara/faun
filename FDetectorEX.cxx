#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>

#include "FDetectorEX.h"
#include "FMP.h"

#include "TList.h"

FDetectorEX::FDetectorEX() {
}

FDetectorEX::~FDetectorEX() {
}

void FDetectorEX::Fill(std::vector<short> packet, short cid[2], short clk[2]) {
  int low[1024];
  int high[1024];
  for(int i=0; i!=1024; ++i)
    fMinipads[i]->Fill(packet[low[i]],packet[high[i]]);
  fCIDFor = cid[0];
  fCIDBck = cid[1];
  fClkFor = clk[0];
  fClkBck = clk[1];
}

TList* FDetectorEX::Init(){ //FDetector has a funtion named Init which returns TList
  TList *output = new TList();
  output->SetName( "EX" );
  for(int idx=0; idx!=1024; ++idx) {
    if( (idx/128)%2==0 ) {
      fMinipadX.push_back( ((idx%128)%32)*0.21 -32.55 );
      fMinipadY.push_back( ((idx%128)/32)*1.52 -5.32  );
    } else {
      fMinipadX.push_back( ((idx%128)/32)*1.52 -5.32  );
      fMinipadY.push_back( ((idx%128)%32)*0.21 -32.55 );
    }
    fMinipadZ.push_back( (idx/128)*1.0 );
    fMinipads.push_back( new FMP() );
  }
  InputGains("EX.calib");
  return output;
}

void FDetectorEX::Dump() {
  std::cout << "fMinipads[" << fMinipads.size() << "]" << std::endl;
}

void FDetectorEX::ReadEnergy(){
}

void FDetectorEX::Reset(){
  fCorrupt = false;
  for(int i=0; i!=1024; ++i)
    fMinipads[i]->Reset();
}

void FDetectorEX::DoQA(){
}

void FDetectorEX::InputGains(char*){
  for(int i=0; i!=1024; ++i) {
    double low = 0;
    double high = 0;
    fMinipads[i]->SetGains(low,high);
  }
}

double FDetectorEX::CentroidX(){
  return 0;
}

double FDetectorEX::CentroidY(){
  return 0;
}



