#ifndef __FTRACE_H__
#define __FTRACE_H__
#include <vector>

class FTrace {
 public:
  FTrace();
  virtual ~FTrace();
  void Fill(double val) { fData.push_back(val); }
  void Reset() {fData.clear(); }
  double Signal() { return fSignal; }
  void SetRange(int a,int b,int c,int d) {
    fRange[0]=a; fRange[1]=b; fRange[2]=c; fRange[3]=d;
  }
  virtual int ComputeSignal(std::vector<double>,double res=0);
  double GetDataSlice(unsigned int idx) {return idx<fData.size()?fData[idx]:0;}

 protected:
  std::vector<double> fData;
  unsigned int fRange[4];
  double fSignal;
};

#endif /*__FTRACE_H__*/
