#include "FTrace.h"
#include <vector>
#include <cmath>
#include <iostream>

FTrace::FTrace() {
  fSignal = 0;
  fBaseline = 0;
}

FTrace::~FTrace(){
}

float FTrace::ComputeBaseline(const std::vector<float>& noise) {
  // compute baseline based on specified range in trace
  // if noise is present use it
  fBaseline = 0;
  if(fRange[0]>fRange[1]) return 0;
  if(fData.size()<fRange[1]-1) return 0;
  for( std::vector<float>::size_type i = fRange[0]; i != fRange[1]; ++i ) {
    fBaseline += fData[i];
    if(noise.size()==fData.size())
      fBaseline -= noise[i];
  }
  fBaseline /= fRange[1]-fRange[0];
  return fBaseline;
}

float FTrace::ComputeSignal(const std::vector<float>& noise, const float& residual) {
  // compute signal based on specified range in trace
  // if noise is present use it
  fSignal = 0;
  if(fRange[2]>fRange[3]) return 0;
  if(fData.size()<fRange[3]-1) return 0;
  for( std::vector<float>::size_type i = fRange[2]; i != fRange[3]; ++i ) {
    fSignal += fData[i];
    if(noise.size()==fData.size())
      fSignal -= noise[i];
  }
  fSignal /= fRange[3]-fRange[2];
  fSignal -= fBaseline;
  //std::cout << "<>|" << fSignal;
  fSignal -= residual;
  //std::cout << "==>" << fSignal << std::endl;
  return fSignal;
}
