
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
  std::vector<double> fCommonNoise0;
  double fCommonNoise1;
  bool fCorrupt;
  double fEnergy;
  TH2D *fQA_APDS;
  TH2D *fQA_SIGNALS;
  TH2D *fQA_ENERGIES;
  TH2D *fQA_CENTROID;
  TH2D *fQA_ENERGY;
  TProfile *fQA_CMN0;
  TH2D *fQA_CMN0_N;
  TH2D *fQA_CMN0_S[23];
  TH1D *fQA_CMN1_N;
  TH1D *fQA_CMN1;

  void EstimateCommonNoise0();
  void EstimateCommonNoise1();
};

#endif /*__FDetectorEX_H__*/
