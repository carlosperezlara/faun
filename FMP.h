#ifndef __FMP_H__
#define __FMP_H__

class FMP {
 public:
  FMP();
  ~FMP();
  void Fill(short h, short l) { fDataHigh=h; fDataLow=l; }
  void SetGain(double g) { fGain = g; }
  double Energy() { return fEnergy; }

 protected:
  short fDataHigh;
  short fDataLow;
  double fGain;
  double fEnergy;
};

#endif /*__FMP_H__*/
