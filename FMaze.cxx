#include <iostream>
#include "FMaze.h"
#include "FDetectorBC.h"
#include "FDetectorMPC.h"
#include "FDetectorEX.h"
#include "FAPD.h"

#include "TList.h"

FMaze::FMaze() {
  fBC = new FDetectorBC();
  fMPC = new FDetectorMPC();
  fEX = new FDetectorEX();
}
FMaze::~FMaze() {
  delete fMPC;
  delete fBC;
  delete fEX;
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
  std::cout << "EX detector" << std::endl;
  std::cout << "===========" << std::endl;
  fEX->Dump();
  std::cout << std::endl;
}

void FMaze::Init() {
  fForest = new TList();
  fForest->SetOwner();
  fForest->Add( fBC->Init() );
  fForest->Add( fMPC->Init() );
  fForest->Add( fEX->Init() );
  TList *myQA = new TList();
  fQA_Events = new TH1D("fQA_Events","fQA_Events",4,-0.5,3.5);
  myQA->Add( fQA_Events );
  fQA_BC_MPC = new TH2D("fQA_BC_MPC",";BC Signal; MPC Energy (GeV)",250,0,250,250,0,60);
  myQA->Add( fQA_BC_MPC );
  fQA_EX_MPC = new TH2D("fQA_EX_MPC",";EX Energy (GeV); MPC Energy (GeV)",250,0,0.3,250,0,60);
  myQA->Add( fQA_EX_MPC );
  fQA_TotalE = new TH1D("fQA_TotalE","MPC+EX Energy (GeV); Energy (GeV)",250,0,60);
  myQA->Add( fQA_TotalE );
  fForest->Add( myQA );
  fMPCtree = new TNtuple("MPCtree","MPCtree",
			 "c0:c1:c2:c3:c4:c5:c6:c7:c8:c9:c10:c11:c12:c13:c14:c15:c16:c17:c18:c19:c20:c21:c22:c23");
  fForest->Add( fMPCtree );
  Dump();
}

void FMaze::WriteOutput(const char *out) {
  fForest->SaveAs(out);
}

void FMaze::Reset() {
  fBC->Reset();
  fMPC->Reset();
  fEX->Reset();
}

void FMaze::Exec() {
  fQA_Events->Fill( 0 );
  fBC->Read();
  fMPC->ReadEnergy();
  fEX->ReadEnergy();

  if(fBC->Corrupt()) return;
  if(fMPC->Corrupt()) return;
  if(fEX->IsPresent() && fEX->Corrupt()) return;

  //Float_t *crys = new Float_t [23];
  //for(int i=0; i!=23; ++i)
  //  crys[i] = fMPC->GetCrystal(i)->Signal();
  //fMPCtree->Fill( crys );
  //delete[] crys;

  fQA_Events->Fill( 1 );
  fQA_BC_MPC->Fill( fBC->Signal(), fMPC->Energy() );
  fQA_EX_MPC->Fill( fEX->Energy(), fMPC->Energy() );
  fQA_TotalE->Fill( fEX->Energy() + fMPC->Energy() ) ;
  double mpc_e = fMPC->Energy();
  double bc_e = fBC->Signal();

  if(0) { // particle selector
    float ox = 45;
    float oy = 6;
    float ax = 15;
    float ay = 3;
    int ntrk = 2; // number of tracks
    float xsgm = (bc_e-ntrk*ox)*(bc_e-ntrk*ox)/ax/ax;
    float ysgm = (mpc_e-ntrk*oy)*(mpc_e-ntrk*oy)/ay/ay;
    if( (xsgm+ysgm) > 1.0 ) return;
    printf("PASSED\n");
  }

  fQA_Events->Fill( 2 );
  fBC->DoQA();
  fMPC->DoQA();
  fEX->DoQA();
}
