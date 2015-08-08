#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include "pmonitor.h"
#include "FTMUnpacker.h"
#include "FTrace.h"
#include "FAPD.h"
#include "FMP.h"
#include "FDetectorMPC.h"

int first_event_nr = 10;

using namespace std;

int process_event (Event * e) {
  if (e->getEvtType() > 7) return 0;
  int thisRun = e->getRunNumber();

  //Use nr-th event as starting point of process_event
  // and skip every event before that
  static int skipentries=0;
  skipentries++;
  static bool SetSkipEntry = true;
  if(SetSkipEntry) {
    cout << " Skipping " << first_event_nr << " events." << endl; 
    SetSkipEntry = false;
  }
  if(e->getEvtSequence() < first_event_nr && e->getEvtType()==1)
    return 0;


  // COUNTING
  static int entries=0;
  entries++;
  if(entries%1000==0) {
    cout << "Processed " << entries;
    time_t now = time(0);
    char* dt = ctime(&now);      
    cout << " at " << dt << endl;
  }

  // Here is a data integrity check for the SRS system
  static bool PriorPassed = true;
  Packet *p1010 = e->getPacket(1010);
  if(p1010) {
    int nhb = p1010->iValue(0,"NHYBRIDS");
    if(nhb < 8) {
      if(PriorPassed) {
	cout << "Wrong number of hybrids: " << nhb << ", abort event." << endl;
      } else {
	cout << ".";
      }
      delete p1010;
      PriorPassed = false;
      return 0;
    }
    for(int i=1; i!=nhb-1; ++i) {
      if(p1010->iValue(0,"NSAMPLES") != p1010->iValue(i,"NSAMPLES")) {
	if(PriorPassed) {
	  cout << "Data length mismatch: " << p1010->iValue(i,"NSAMPLES")
	       << " out of " << p1010->iValue(0,"NSAMPLES")
	       << " on hybrid " << i << " abort event." << endl;
	} else {
	  cout << ".";
	}
	delete p1010;
	PriorPassed = false;
	return 0;
      }
    }
    if(!PriorPassed) cout << endl;
    PriorPassed = true;
    for (int i=0; i!=nhb; i++) {
      if (p1010->iValue(i,"NSAMPLES")<24) {
	cout << "ERROR in Sample Length, abort run" << endl;
	exit(0);
      }
    }
    static bool FirstEvent = true;
    if(FirstEvent) {
      cout << "Reports for " << nhb << " Hybrids:" << endl;
      for (int i=0; i!=nhb; ++i) {
	cout << i << ": " << p1010->iValue(i,"NSAMPLES") << endl;
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
  FTrace theSCI;
  Packet *p1020 = e->getPacket(1020);
  if(p1020) {
    int samples = p1020->iValue(0, "SAMPLES"); 
    for( int i=0; i!=samples; ++i )
      theSCI.Fill(-p1020->rValue(i,0));
  }
  if(p1020) delete p1020;

  // MPC
  FDetectorMPC theMPC;
  Packet *p2001 = e->getPacket(2001);
  if(p2001) {
    for(int c=0; c!=23; ++c) {
      FAPD crys;
      for (int i = 0; i< 1024; i++) {
	crys.Fill( -p2001->iValue(i,c) );
      }
      theMPC.AddCrystal(0,0,crys);
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
    printf("MPCEX_1030 NW %d\n",nw);
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
    //  disambiguation 
    bool unique = true;
    vector<int> starts;
    starts.clear();
    starts.push_back(32);
    int marker = bytes[31];
    int index = 31+128;
    while (index<4*nw && unique && starts.size()<8) {
      int count =0;
      int match =0;
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
    //if(unique) {
    //  int ctr=0;
    //  for(int i=0; i!=starts.size(); ++i) {
    //	for(int j=starts[i]; j!=starts[i]+128; ++j) {
    //	  AMinipad::Raw[ctr] = bytes[j];
    //	  AMinipad::Clock[ctr] = clock;
    //	  ctr++;
    //	}
    //}
    cellIDF = b2G( bytes[31] );
  } else {
    cout << "Not Unique Result from USB0...abort event " << entries << endl;
    delete p1030;
    return 0;
  }
  if(p1030) delete p1030;
  
  Packet *p1031 = e->getPacket(1031);
  int cellIDB=-1;
  if(p1031) {
    int data[4000];
    int nw;
    p1031->fillIntArray(data,4000,&nw,"RAW");
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
	cout << " " << i;
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
    //if(unique) {
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
    cout << "Not Unique Result from USB1...abort event "
	 << entries << endl;
    delete p1031;
    return 0;
  }
  if(p1031) delete p1031;
  if(cellIDF==63) return 0;
  if(cellIDB==63) return 0;
  return 0;  
}

int b2G(int num) {
  return (num >> 1) ^ num;
}

