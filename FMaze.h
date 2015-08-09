#ifndef __FMAZE_H__
#define __FMAZE_H__

#include "FTrace.h"
#include "FDetectorMPC.h"

#include "TList.h"

class FMaze{
 public:
  FMaze();
  ~FMaze();
  FTrace* SCI() {return fSCI;}
  FDetectorMPC* MPC() {return fMPC;}
  void Init();
  void Reset();
  void Exec();
  void Dump();
  void WriteOutput(const char *);

 protected:
  TList *fForest;
  FTrace *fSCI;
  FDetectorMPC *fMPC;
};

#endif /*__FMAZE_H__*/
