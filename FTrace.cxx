#include "FTrace.h"
#include <vector>
#include <cmath>

FTrace::FTrace() {
  fSignal = -1;
}

FTrace::~FTrace(){
}

void FTrace::ComputeSignal(std::vector<double> noise){
  fSignal = -1;
  if(fRange[0]>fRange[1]) return;
  if(fRange[2]>fRange[3]) return;
  if(fRange[1]>fRange[3]) return;
  if(fData.size()<fRange[3]-1) return;
  if(fData.size()<noise.size()) return;
  double a = 0;
  for( std::vector<double>::size_type i = fRange[0]; i != fRange[1]; ++i )
    a  += fData[i] - noise[i];
  a = a / (fRange[1]-fRange[0]);
  double b = 0;
  for( std::vector<double>::size_type i = fRange[2]; i != fRange[3]; ++i )
    b  += fData[i] - noise[i];
  b = b / (fRange[3]-fRange[2]);
  fSignal = b - a;
  return;
}
