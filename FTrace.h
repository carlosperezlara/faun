#ifndef __FTRACE_H__
#define __FTRACE_H__
#include <vector>

class FTrace {
 public:
  FTrace();
  virtual ~FTrace();
  void Fill(float val) { fData.push_back(val); }
  void Reset() {fData.clear(); }
  unsigned int NumberOfSlices() {return fData.size(); }
  float Signal() { return fSignal; }
  float Baseline() { return fBaseline; }
  void SetRange(int a,int b,int c,int d) {
    fRange[0]=a; fRange[1]=b; fRange[2]=c; fRange[3]=d;
  }
  virtual float ComputeBaseline(const std::vector<float>& noise);
  virtual float ComputeSignal(const std::vector<float>& noise, const float& res=0.0);
  float GetDataSlice(unsigned int idx) {return idx<fData.size()?fData[idx]:0;}

 protected:
  std::vector<float> fData;
  unsigned int fRange[4];
  float fSignal;
  float fBaseline;
};

#endif /*__FTRACE_H__*/
