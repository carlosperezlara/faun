#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>

#include "FDetectorEX.h"
#include "FMP.h"

#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

FDetectorEX::FDetectorEX() {
}

FDetectorEX::~FDetectorEX() {
}

void FDetectorEX::Dump() {
  std::cout << "fMinipads[" << fMinipads.size() << "]" << std::endl;
}

void FDetectorEX::ReadEnergy(){
}

void FDetectorEX::Reset(){
}

void FDetectorEX::DoQA(){
}

void FDetectorEX::InputGains(char*){
}

double FDetectorEX::CentroidX(){
  return 0;
}

double FDetectorEX::CentroidY(){
  return 0;
}



