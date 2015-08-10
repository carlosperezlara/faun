#include <iostream>
#include "FMaze.h"
#include "FDetectorBC.h"
#include "FDetectorMPC.h"

#include "TList.h"

FMaze::FMaze() {
  fBC = new FDetectorBC();
  fMPC = new FDetectorMPC();
}
FMaze::~FMaze() {
  delete fMPC;
  delete fBC;
}

void FMaze::Dump() {
  std::cout << "MPC detector" << std::endl;
  std::cout << "============" << std::endl;
  fMPC->Dump();
  std::cout << std::endl;
  std::cout << "BC detector" << std::endl;
  std::cout << "===========" << std::endl;
  fBC->Dump();
  std::cout << std::endl;
}

void FMaze::Init() {
  fForest = new TList();
  fForest->SetOwner();
  fForest->Add( fBC->Init() );
  fForest->Add( fMPC->Init() );
  TList *myQA = new TList();
  fQA_BC_MPC = new TH2D("fQA_BC_MPC",";BC Signal; MPC Energy",250,0,250,120,0,60);
  myQA->Add( fQA_BC_MPC );

  fForest->Add( myQA );
}

void FMaze::WriteOutput(const char *out) {
  fForest->SaveAs(out);
}

void FMaze::Reset() {
  fBC->Reset();
  fMPC->Reset();
  fQA_BC_MPC->Fill( fBC->Signal(), fMPC->Energy() );
}

void FMaze::Exec() {
  fBC->Read();
  if(fBC->Corrupt()) return;
  fMPC->ReadEnergy();
  if(fMPC->Corrupt()) return;

  fBC->DoQA();
  fMPC->DoQA();

  fQA_BC_MPC->Fill( fBC->Signal(), fMPC->Energy() );
}
