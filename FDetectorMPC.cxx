#include "FDetectorMPC.h"
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>

FDetectorMPC::FDetectorMPC() {
}

FDetectorMPC::~FDetectorMPC() {
}

void FDetectorMPC::EstimateCommonNoise() {
  fCommonNoise.clear();
  for(std::vector<FAPD>::size_type i = 0; i!=fCrystals.size(); ++i) {
    fCommonNoise.push_back(0);
  }
  return;
}

void FDetectorMPC::ReadEnergy() {
  EstimateCommonNoise();
  fEnergy = 0;
  for(std::vector<FAPD>::size_type i = 0; i!=fCrystals.size(); ++i) {
    fCrystals[i].SetRange(0,75,125,675);
    fCrystals[i].ComputeSignal(fCommonNoise);
    fEnergy += fCrystals[i].Energy();
  }
}

double FDetectorMPC::CentroidX() {
  std::vector<FAPD>::size_type size = fCrystals.size();
  if(size<1) return -999.;
  double x = 0;
  for(std::vector<FAPD>::size_type i = 0; i!=size; ++i) {
    x += fCrystals[i].Energy()*fCrystalX[i];
  }
  return x/size;
}

double FDetectorMPC::CentroidY() {
  std::vector<FAPD>::size_type size = fCrystals.size();
  if(size<1) return -999.;
  double x = 0;
  for(std::vector<FAPD>::size_type i = 0; i!=size; ++i) {
    x += fCrystals[i].Energy()*fCrystalY[i];
  }
  return x/size;
}

void FDetectorMPC::InputGains(char *file) {
  std::ifstream calib;
  calib.open(file);
  double tmp;
  for(std::vector<FAPD>::size_type i=0; i!=fCrystals.size(); ++i) {
    calib >> tmp;
    fCrystals[i].SetGain( tmp );
    printf("Gain read for DetectorMPC || crystal %lu <=> %f\n",i,tmp);
    if(!calib.good()) break;
  }
  calib.close();
  return;
}
