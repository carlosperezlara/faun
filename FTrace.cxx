#include "FTrace.h"
#include <vector>
#include <cmath>

FTrace::FTrace() {
  fSignal = -1;
}

FTrace::~FTrace(){
}

int FTrace::ComputeSignal(std::vector<double> noise, double residual){
  fSignal = -1;
  if(fRange[0]>fRange[1]) return 9999;
  if(fRange[2]>fRange[3]) return 9999;
  if(fRange[1]>fRange[3]) return 9999;
  if(fData.size()<fRange[3]-1) return 9999;
  if(fData.size()<noise.size()) return 9999;
  double a;
  for( std::vector<double>::size_type i = fRange[0]; i != fRange[1]; ++i )
    a  += fData[i] - noise[i];
  double b = 0;
  for( std::vector<double>::size_type i = fRange[2]; i != fRange[3]; ++i )
    b  += fData[i] - noise[i];
  fSignal = b/(fRange[3]-fRange[2]) - a/(fRange[1]-fRange[0]) - residual;
  return a;
}
