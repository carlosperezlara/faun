#ifndef __FDetectorMPC_H__
#define __FDetectorMPC_H__

#include <vector>

#include "FAPD.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TList.h"

class FDetectorMPC{
 public:
  FDetectorMPC();
  ~FDetectorMPC();
  TList* Init();
  FAPD* GetCrystal(unsigned int idx) {return idx<fCrystals.size()?fCrystals[idx]:NULL;}
  double Energy() {return fEnergy;}
  double CentroidX();
  double CentroidY();
  void ReadEnergy();
  void InputGains();
  void Reset();
  void Dump();
  void DoQA();
  bool Corrupt() {return fCorrupt;}

 protected:
  bool fDOQA;
  std::vector<FAPD*> fCrystals;
  std::vector<double> fCrystalX;
  std::vector<double> fCrystalY;
  std::vector<float> fCommonNoise0;
  float fCommonNoise1;
  bool fCorrupt;
  double fEnergy;
  float fBaseline[23];
  int fOrd[23];
  TH2D *fQA_BL;
  TH2D *fQA_SIGNALS;
  TH2D *fQA_ENERGIES;
  TH2D *fQA_CENTROID;
  TH2D *fQA_ENERGY;
  TH2D *fQA_CMN0;
  TH1D *fQA_CMN0_N;
  TH2D *fQA_F100[23];
  TH2D *fQA_S[23];
  TH2D *fQA_CMN0_S[23];
  TH1D *fQA_CMN1_N;
  TH1D *fQA_CMN1;

  void ComputeSignals();
  void Sort();
  void EstimateCommonNoise0();
  void EstimateCommonNoise1();
  int getneig(int b[23], int nb);
};

#endif /*__FDetectorMPC_H__*/
