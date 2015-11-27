#ifndef __FADC_H__
#define __FADC_H__

class FADC {
 public:
  FADC();
  ~FADC();
  short Get() { return fADC; }
  double GetGain() { return fGain; }
  double GetPMean() { return fPMean; }
  double GetPSigma() { return fPSigma; }
  double GetCMN() { return fCMN; }
  double RawPed() { return fADC-fPMean; }

  void Set(short v) { fADC=v; }
  void SetGain(double g) { fGain = g; }
  void SetPedestal(double a, double b) { fPMean=a; fPSigma=b; }
  void SetCMN(double v) { fCMN=v; }

  double Energy() { return (fADC-fPMean-fCMN)*fGain; }
  double SEnergy() { return fPSigma*fGain; }
  

 protected:
  short fADC;
  double fGain;
  double fPMean;
  double fPSigma;
  double fCMN;
};

#endif /*__FADC_H__*/
