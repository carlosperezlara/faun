#ifndef __FAPD_H__
#define __FAPD_H__
#include "FTrace.h"

class FAPD: public FTrace{
 public:
  FAPD();
  ~FAPD();
  double Energy() {return fSignal*fGain;}

 protected:
  double fGain;
  virtual void ComputeSignal();
  virtual void ComputeTime();
};

#endif /*__FAPD_H__*/
