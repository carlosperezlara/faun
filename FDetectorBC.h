#ifndef __FDetectorBC_H__
#define __FDetectorBC_H__

#include <vector>

#include "FAPD.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TList.h"

class FDetectorBC{
 public:
  FDetectorBC();
  ~FDetectorBC();
  TList* Init();
  FTrace* GetTrace(unsigned int idx) {return idx<fTraces.size()?fTraces[idx]:NULL;}
  void Read();
  void Reset();
  void Dump();
  void DoQA();
  bool Corrupt() {return fCorrupt;}

 protected:
  std::vector<FTrace*> fTraces;
  std::vector<double> fCommonNoise0;
  bool fCorrupt;
  TH2D *fQA_TRACE;
  TH1D *fQA_SIGNALS;
};

#endif /*__FDetectorBC_H__*/
