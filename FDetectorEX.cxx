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

TList* FDetectorEX::Init(){ //FDetector has a funtion named Init which returns TList
  TList *myQA = new TList();
  return myQA;
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



