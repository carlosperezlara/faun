#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>

#include "FDetectorBC.h"
#include "FAPD.h"

#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

FDetectorBC::FDetectorBC() {
}

FDetectorBC::~FDetectorBC() {
}

void FDetectorBC::Dump() {
  std::cout << "fTraces[" << fTraces.size() << "]" << std::endl;
}

TList* FDetectorBC::Init() {
  TList *output = new TList();
  output->SetName( "BC" );
  fTraces.push_back(new FTrace());
  if(1) {
    fQA_TRACE   = new TH2D( "BC_TRACE", ";slice;value", 1024,-0.5,1023.5, 200,0,512);
    fQA_SIGNALS = new TH1D( "BC_SIGNALS", ";Signal",     300,-0.5,299.5);
    output->Add( fQA_TRACE );
    output->Add( fQA_SIGNALS );
  }
  fCommonNoise0.clear();
  for(unsigned int j=0; j!=1024; ++j)
    fCommonNoise0.push_back(0);
  return output;
}

void FDetectorBC::DoQA() {
  if(fQA_TRACE)
    for(unsigned int j=0; j!=1024; ++j)
      fQA_TRACE->Fill(j,fTraces[0]->GetDataSlice(j));
  if(fQA_SIGNALS)
    fQA_SIGNALS->Fill(fTraces[0]->Signal());
}

void FDetectorBC::Reset() {
  fCorrupt = false;
  fTraces[0]->Reset();
}

void FDetectorBC::Read() {
  fTraces[0]->SetRange(0,300,440,460);
  /*double chk = */fTraces[0]->ComputeSignal(fCommonNoise0);
  //std::cout << "BC chk " << chk << std::endl;
}
