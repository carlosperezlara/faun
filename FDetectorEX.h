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
  void Fill(std::vector<short> v, short cid[2], short clk[2]);
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
  std::vector<FMP*>   fMinipads;
  std::vector<short>    fMinipadK;
  std::vector<double> fMinipadX;
  std::vector<double> fMinipadY;
  std::vector<double> fMinipadZ;
  bool fCorrupt;
  double fEnergy;
  short fClkFor;
  short fClkBck;
  short fCIDFor;
  short fCIDBck;

};

#endif /*__FDetectorEX_H__*/
