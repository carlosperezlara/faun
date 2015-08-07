#ifndef __FMINIPAD_H__
#define __FMINIPAD_H__
class FMinipad {
 public:
  FMinipad();
  ~FMinipad();

 protected:
  short fADChigh;
  short fADClow;
  float fSIGhigh;
  float fSIGlow;
  float fEnergy;
};

#endif /*__FMINIPAD_H__*/
