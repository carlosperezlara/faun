#ifndef __FMP_H__
#define __FMP_H__

#include "FADC.h"

class FMP {
 public:
  FMP();
  ~FMP();
  FADC* GetHigh() { return fHigh; }
  FADC* GetLow() { return fLow; }
  void Reset()
  { fHigh->Set(0); fLow->Set(0); fHigh->SetCMN(0); fLow->SetCMN(0); }
  void Fill(short h, short l)
  { fHigh->Set(h); fLow->Set(l); }
  void SetGains(double h, double l)
  { fHigh->SetGain(h); fLow->SetGain(l); }
  void SetPedestals(double hm, double hs, double lm, double ls)
  { fHigh->SetPedestal(hm,hs); fLow->SetPedestal(lm,ls); }
  double Energy();
  double EnergyHigh() { return fHigh->Energy(); }
  double EnergyLow() {  return fLow->Energy(); }

 protected:
  FADC *fHigh;
  FADC *fLow;
};

#endif /*__FMP_H__*/
