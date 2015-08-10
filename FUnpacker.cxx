//#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include <fstream>

#include "pmonitor.h"
#include "Event.h"
#include "FTrace.h"
#include "FAPD.h"
#include "FMP.h"
#include "FDetectorMPC.h"
#include "FDetectorBC.h"
#include "FMaze.h"

using namespace std;

ofstream fileEvType("EvType.log");
ofstream fileHybrid("Hybrid.log");
ofstream fileMPCEX0("MPCEX0.log");
FMaze *FaunsMaze;

int b2G(int num) {
  return (num >> 1) ^ num;
}

int pinit() {
  cout << "pINIT" << endl;
  FaunsMaze = new FMaze();
  FaunsMaze->Init();
  //FaunsMaze->Dump();
}

int plistopen(const char *out) {
  FaunsMaze->WriteOutput( out );
}

int process_event (Event * e) {
  FaunsMaze->Reset();
  if(1) {
    fileEvType.close();
    fileHybrid.close();
    fileMPCEX0.close();
  }
  fileEvType << e->getEvtType() << endl;
  if (e->getEvtType() > 7) return 0;
  int thisRun = e->getRunNumber();
  static int skipentries=20;
  if(skipentries>0) {
    if(skipentries==20) cout << "SKIPPING FIRST: ";
    cout << "X";
    skipentries--;
    //if(e->getEvtSequence() < 20 && e->getEvtType()==1) return 0;
    if(skipentries==0) cout << "EVENTS." << endl;
    return 0;
  }
  // COUNTING
  static int entries=0;
  entries++;
  if(entries%1000==0) {
    cout << "Proccessed " << entries;
    time_t now = time(0);
    char* dt = ctime(&now);      
    cout << " at " << dt << endl;
  }
  // TRACKERS
  Packet *p1010 = e->getPacket(1010);
  if(p1010) {
    int nhb = p1010->iValue(0,"NHYBRIDS");
    fileHybrid << "****" << endl;
    fileHybrid << nhb << endl;
    if(nhb<8) {
      fileHybrid << "Wrong number of hybrids: " << nhb << " ==> abort." << endl;
      delete p1010;
      return 0;
    }
    for(int i=1; i!=nhb-1; ++i) {
      if(p1010->iValue(0,"NSAMPLES") != p1010->iValue(i,"NSAMPLES")) {
	fileHybrid << "Data length mismatch: " << p1010->iValue(i,"NSAMPLES") << " out of " << p1010->iValue(0,"NSAMPLES") << " on hybrid " << i << " ==> abort." << endl;
	delete p1010;
	return 0;
      }
    }
    for(int i=0; i!=nhb; ++i) {
      if(p1010->iValue(i,"NSAMPLES")<24) {
	fileHybrid << "ERROR in Sample Length: " << p1010->iValue(i,"NSAMPLES") << " ==> abort." << endl;
	delete p1010;
	return 0;
      }
    }
    static bool FirstEvent = true;
    if(FirstEvent) {
      fileHybrid << " Reports for " << nhb << " Hybrids:" << endl;
      for (int i=0; i!=nhb; ++i) {
	fileHybrid << i << ": " << p1010->iValue(i,"NSAMPLES") << endl;
      }
      FirstEvent = false;
    }
    /*
    // Read the tracker raw data into the raw array...
    int rindex=0;
    for (int JINX =0; JINX<8; JINX++) {
      for ( int j = 0; j< 128; j++) {
	for (int i=0; i<min(p->iValue(JINX,"NSAMPLES"),24); i++) { // limit to 24 samples...
	  //  For now we'll use a trivial routine to analyze waveform.
	  //  We'll just sum across the samples...
	  AStrip::Raw[rindex].push_back(p->iValue(j,i,JINX));
	}
	rindex++;
      }	      
    }
    */
  }

  // SCINTILLATOR
  FDetectorBC *theBC = FaunsMaze->BC();
  Packet *p1020 = e->getPacket(1020);
  if(p1020) {
    int samples = p1020->iValue(0, "SAMPLES"); 
    for( int i=0; i!=samples; ++i )
      theBC->GetTrace(0)->Fill(-p1020->rValue(i,0));
  }
  if(p1020) delete p1020;

  // MPC
  FDetectorMPC *theMPC = FaunsMaze->MPC();
  Packet *p2001 = e->getPacket(2001);
  if(p2001) {
    for(int c=0; c!=23; ++c) {
      for (int i = 0; i< 1024; i++) {
	theMPC->GetCrystal(c)->Fill( p2001->iValue(i,c) ); // prev -1
      }
    }
  }
  if(p2001) delete p2001;

  // EX
  Packet *p1030 = e->getPacket(1030);
  int cellIDF=-1;
  if(p1030) {
    int data[4000];
    int nw;
    p1030->fillIntArray(data,4000,&nw,"RAW");
    fileMPCEX0 << "******" << endl;
    fileMPCEX0 << "MPCEX_1030 NW " << nw << endl;
    int INDEXER=0;
    int bytes[16000];
    for(int i=0; i!=nw; ++i) {
      int n = data[i];
      bytes[INDEXER+3] = (n >> 24) & 0xFF;
      bytes[INDEXER+2] = (n >> 16) & 0xFF;
      bytes[INDEXER+1] = (n >> 8) & 0xFF;
      bytes[INDEXER+0] = n & 0xFF;
      INDEXER+=4;
    }
    fileMPCEX0 << "0: ";
    for(int i=0; i<nw*4; ++i) {
      fileMPCEX0 << hex << setfill('0') << setw(2) << bytes[i] << " ";
      if((i+1)%16==0)
	fileMPCEX0 << endl << dec << i+1 << ":" << hex;
    }
    fileMPCEX0 << dec << endl;
    bool unique = true;
    vector<int> starts;
    starts.clear();
    starts.push_back(32);
    int marker = bytes[31];
    int index = 31+128;
    while (index<4*nw && unique && starts.size()<8) {
      int count=0;
      int match=0;
      for (int i=index; i!=index+7; ++i) {
	if(bytes[i]==marker) {
	  match = i;
	  count++;
	}
      }
      if(count!=1) {
	unique = false;
      } else {
	starts.push_back(match+1);
	index = match+128;
      }
    }
    if(starts.size()!=8) unique=false;
    int clock = bytes[4*nw-8] & 0x7;
    //clock = (clock+3)%8;
    if(unique) {
    //  int ctr=0;
    //  for(int i=0; i!=starts.size(); ++i) {
    //	  for(int j=starts[i]; j!=starts[i]+128; ++j) {
    //	    AMinipad::Raw[ctr] = bytes[j];
    //	    AMinipad::Clock[ctr] = clock;
    //	    ctr++;
    //	  }
    //  }
      cellIDF = b2G( bytes[31] );
    } else {
      fileMPCEX0 << "Not Unique Result from USB0...abort event " << entries << endl;
      delete p1030;
      return 0;
    }
    delete p1030;
  }
  if(cellIDF==63) {
    fileMPCEX0 << "CELLID == 63 from USB0... ==> abort." << endl;
    return 0;
  }
  
  Packet *p1031 = e->getPacket(1031);
  int cellIDB=-1;
  if(p1031) {
    int data[4000];
    int nw;
    p1031->fillIntArray(data,4000,&nw,"RAW");
    fileMPCEX0 << "******" << endl;
    fileMPCEX0 << "MPCEX_1031 NW " << nw << endl;
    int INDEXER=0;
    int bytes[16000];
    for (int i=0; i<nw; i++) {
      int n = data[i];
      bytes[INDEXER+3] = (n >> 24) & 0xFF;
      bytes[INDEXER+2] = (n >> 16) & 0xFF;
      bytes[INDEXER+1] = (n >> 8) & 0xFF;
      bytes[INDEXER+0] = n & 0xFF;
      INDEXER+=4;
    }
    fileMPCEX0 << "0: ";
    for(int i=0; i<nw*4; ++i) {
      fileMPCEX0 << hex << setfill('0') << setw(2) << bytes[i] << " ";
      if((i+1)%16==0)
	fileMPCEX0 << endl << dec << i+1 << ":" << hex;
    }
    fileMPCEX0 << dec << endl;
    bool unique = true;
    vector<int> starts;
    starts.clear();
    starts.push_back(32);
    int marker = bytes[31];
    int index = 31+128;
    while (index<4*nw && unique && starts.size()<8) {
      int count =0;
      int match =0;
      for(int i=index; i<index+7; i++) {
	if(bytes[i]==marker) {
	  match = i;
	  count++;
	}
      }
      if(count!=1) {
	unique = false;
      } else {
	starts.push_back(match+1);
	index = match+128;
      }
    }
    if(starts.size()!=8) unique=false;
    int clock = bytes[4*nw-8] & 0x7;
    //clock = (clock+6)%8;
    if(unique) {
    //  int ctr=1024;
    //  for (int i=0; i<starts.size(); i++) {
    //	for (int j=starts[i]; j<starts[i]+128; j++) {
    //	  AMinipad::Raw[ctr] = bytes[j];
    //	  AMinipad::Clock[ctr] = clock;
    //	  ctr++;
    //	}
    //}
      cellIDB = b2G( bytes[31] );
    } else {
      fileMPCEX0 << "Not Unique Result from USB1...abort event " << entries << endl;
      delete p1031;
      return 0;
    }
    delete p1031;
  }
  if(cellIDF==63) {
    fileMPCEX0 << "CELLID == 63 from USB1... ==> abort." << endl;
    return 0;
  }

  // all unpacked
  FaunsMaze->Exec();

  return 0;
}
