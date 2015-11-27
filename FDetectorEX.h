#ifndef __FDetectorEX_H__
#define __FDetectorEX_H__

#include <vector>

#include "FMP.h"

#include "TProfile.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TList.h"

class FDetectorEX{
 public:
  FDetectorEX();
  ~FDetectorEX();
  TList* Init();
  void Fill(std::vector<short> v, short cid[2], short clk[2]);
  FMP* GetMinipad(unsigned int idx) {return idx<fMinipads.size()?fMinipads[idx]:NULL;}
  double Energy(int lyr=-1);
  double CentroidX(int);
  double CentroidY(int);
  void ReadEnergy();
  void InputGains();
  void EstimateCommonNoise();
  void Reset();
  void Dump();
  void DoQA();
  bool Corrupt() {return fCorrupt;}
  bool IsPresent() {return fIsPresent;}

 protected:
  bool fDOQA;
  bool fIsPresent;
  std::vector<FMP*>   fMinipads;
  std::vector<int>    fLowMap;
  std::vector<int>    fHighMap;
  std::vector<double> fMinipadX;
  std::vector<double> fMinipadY;
  std::vector<double> fMinipadZ;
  std::vector<double> fCommonNoise0;
  std::vector<double> fCommonNoise1;
  std::vector<int>    fStateMap;
  bool fCorrupt;
  double fEnergy;
  short fClkFor;
  short fClkBck;
  short fCIDFor;
  short fCIDBck;

  TH2D *fQA_CellID;
  TH2D *fQA_Clock;
  TH2D *fQA_RawMap;
  TH2D *fQA_PedMap;
  TH2D *fQA_RawPedMap;
  TH2D *fQA_MIPMap;
  TH2D *fQA_EneMap;
  TH3D *fQA_HL;
  TH2D *fQA_EneLyr;
  TH1D *fQA_Energy;
  TH3D *fQA_Centroid;
  TH3D *fQA_Layer0;
  TH3D *fQA_Layer1;
  TProfile *fQA_CommonNoise;
  TH2D *fQA_NCN;
};

#endif /*__FDetectorEX_H__*/
