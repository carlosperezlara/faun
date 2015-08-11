#ifndef __FMP_H__
#define __FMP_H__

#include "FADC.h"

class FMP {
 public:
  FMP();
  ~FMP();
  void Reset() { fHigh.Set(0); fLow.Set(0); }
  void Fill(short h, short l) { fHigh.Set(h); fLow.Set(l); }
  void SetGains(double h, double l) { fHigh.SetGain(h); fLow.SetGain(l); }
  double Energy();

 protected:
  FADC fHigh;
  FADC fLow;
};

#endif /*__FMP_H__*/
