#include <iostream>
#include "FMaze.h"
#include "FTrace.h"
#include "FDetectorMPC.h"

#include "TList.h"

FMaze::FMaze() {
  fMPC = new FDetectorMPC();
  fSCI = new FTrace();
}
FMaze::~FMaze() {
  delete fMPC;
  delete fSCI;
}

void FMaze::Dump() {
  std::cout << "MPC detector" << std::endl;
  std::cout << "============" << std::endl;
  fMPC->Dump();
  std::cout << std::endl;
}

void FMaze::Init() {
  fForest = new TList();
  fForest->SetOwner();
  fForest->Add( fMPC->Init() );
}

void FMaze::WriteOutput(const char *out) {
  fForest->SaveAs(out);
}

void FMaze::Reset() {
  fSCI->Reset();
  fMPC->Reset();
}

void FMaze::Exec() {
  fMPC->ReadEnergy();
  if(fMPC->Corrupt()) return;

  fMPC->DoQA();
}
