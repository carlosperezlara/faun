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
  void InputGains(char *);
  void Reset();
  void Dump();
  void DoQA();

 protected:
  std::vector<FAPD*> fCrystals;
  std::vector<double> fCrystalX;
  std::vector<double> fCrystalY;
  std::vector<double> fCommonNoise;
  double fEnergy;
  TH2D *fQA_APDS;
  TH2D *fQA_SIGNALS;
  TH2D *fQA_ENERGIES;
  TH2D *fQA_CENTROID;
  TH1D *fQA_ENERGY;
  TProfile *fQA_CMN;

  void EstimateCommonNoise();
};

#endif /*__FDetectorMPC_H__*/
