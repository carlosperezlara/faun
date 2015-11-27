#ifndef __FAPD_H__
#define __FAPD_H__
#include "FTrace.h"

class FAPD: public FTrace{
 public:
  FAPD();
  ~FAPD();
  double Energy() { return fSignal*fGain; }
  void SetGain(double g) { fGain=g; }
  double GetGain() { return fGain; }
 protected:
  double fGain;
};

#endif /*__FAPD_H__*/
