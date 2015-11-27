#ifndef __FMAZE_H__
#define __FMAZE_H__

#include "FDetectorBC.h"
#include "FDetectorMPC.h"
#include "FDetectorEX.h"

#include "TList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TNtuple.h"

class FMaze{
 public:
  FMaze();
  ~FMaze();
  FDetectorBC* BC() {return fBC;}
  FDetectorMPC* MPC() {return fMPC;}
  FDetectorEX* EX() {return fEX;}
  void Init();
  void Reset();
  void Exec();
  void Dump();
  void WriteOutput(const char *);

 protected:
  TList *fForest;
  FDetectorBC *fBC;
  FDetectorMPC *fMPC;
  FDetectorEX *fEX;
  TNtuple *fMPCtree;
  TH1D *fQA_Events;
  TH2D *fQA_BC_MPC;
  TH2D *fQA_EX_MPC;
  TH1D *fQA_TotalE;
};

#endif /*__FMAZE_H__*/
