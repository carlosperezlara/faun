#ifndef __FADC_H__
#define __FADC_H__

class FADC {
 public:
  FADC();
  ~FADC();
  short Get() {return fADC;}
  void Set(short v) { fADC=v; }
  void SetGain(double g) { fGain = g; }
  double Energy() { return fADC*fGain; }

 protected:
  short fADC;
  double fGain;
};

#endif /*__FADC_H__*/
