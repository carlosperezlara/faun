#include "FMP.h"

FMP::FMP(){
  Reset();
}

FMP::~FMP(){
}

double FMP::Energy() {
  short hen = fHigh.Get();
  short len = fLow.Get();
  if(hen<100) return fHigh.Energy();
  if(len>50) return fLow.Energy();
  return ( fHigh.Energy() + fLow.Energy() ) / 2.0;
}
