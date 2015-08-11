#ifndef __FMAZE_H__
#define __FMAZE_H__

#include "FDetectorBC.h"
#include "FDetectorMPC.h"
#include "FDetectorEX.h"

#include "TList.h"
#include "TH2D.h"

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
  TH2D *fQA_BC_MPC;
};

#endif /*__FMAZE_H__*/
