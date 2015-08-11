#ifndef __FDetectorEX_H__
#define __FDetectorEX_H__

#include <vector>

#include "FMP.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TList.h"

class FDetectorEX{
 public:
  FDetectorEX();
  ~FDetectorEX();
  TList* Init();
  FMP* GetMinipad(unsigned int idx) {return idx<fMinipads.size()?fMinipads[idx]:NULL;}
  double Energy() {return fEnergy;}
  double CentroidX();
  double CentroidY();
  void ReadEnergy();
  void InputGains(char *);
  void Reset();
  void Dump();
  void DoQA();
  bool Corrupt() {return fCorrupt;}

 protected:
  std::vector<FMP*> fMinipads;
  std::vector<double> fMinipadX;
  std::vector<double> fMinipadY;
  std::vector<int> fMinipadZ; //Layer index
  bool fCorrupt;
  double fEnergy;
};

#endif /*__FDetectorEX_H__*/
