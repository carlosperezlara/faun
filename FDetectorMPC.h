#ifndef __FDetectorMPC_H__
#define __FDetectorMPC_H__

#include "FAPD.h"
#include <vector>

class FDetectorMPC{
 public:
  FDetectorMPC();
  ~FDetectorMPC();
  void AddCrystal(double x, double y, FAPD c) {
    fCrystalX.push_back(x);
    fCrystalY.push_back(y);
    fCrystals.push_back(c); }
  double Energy() {return fEnergy;}
  double CentroidX();
  double CentroidY();
  void ReadEnergy();
  void InputGains(char *);
 protected:
  std::vector<FAPD> fCrystals;
  std::vector<double> fCrystalX;
  std::vector<double> fCrystalY;
  std::vector<double> fCommonNoise;
  double fEnergy;
  void EstimateCommonNoise();
};

#endif /*__FDetectorMPC_H__*/
