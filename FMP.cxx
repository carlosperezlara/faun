#include "FMP.h"

FMP::FMP(){
  fHigh = new FADC();
  fLow = new FADC();
  Reset();
}

FMP::~FMP(){
}

double FMP::Energy() {
  short hen = fHigh->Get();
  short len = fLow->Get();
  if(hen<180) return fHigh->Energy();
  return fHigh->Energy();
  //if(len>50) return fLow->Energy();
  //return ( fHigh->Energy() + fLow->Energy() ) / 2.0;
}
