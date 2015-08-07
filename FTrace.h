#ifndef __FTRACE_H__
#define __FTRACE_H__
#include <vector>

class CTrace {
  public;
  FTrace();
  ~FTrace();
  double Signal() {return fSignal;}
  double Time() {return fTime;}

  protected;
  std::vector<double> voltage;
  void ComputeSignal();
  void ComputeTime();
  double fSignal;
  double fTime;
};

#endif /*__FTRACE_H__*/
