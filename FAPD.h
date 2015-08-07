#include "FTrace.h"

class FAPD::FTrace {
 public:
  FAPD();
  ~FAPD();
  double Energy() {return fSignal*fGain;}

 protected:
  double fGain;
  virtual double ComputeSignal();
  virtual double ComputeTime();
};
