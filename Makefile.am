AUTOMAKE_OPTIONS = foreign

INCLUDES = -I$(includedir) -I$(ONLINE_MAIN)/include -I$(OFFLINE_MAIN)/include  -I`root-config --incdir`

lib_LTLIBRARIES = \
  libFaun.la

noinst_HEADERS = \
  FaunLinkDef.h 

include_HEADERS = \
	FTMUnpacker.h \
	FTrace.h \
	FAPD.h \
	FMP.h \
	FDetectorMPC.h

libFaun_la_SOURCES = \
	FTMUnpacker.cxx \
	FTrace.cxx \
	FAPD.cxx \
	FMP.cxx \
	FDetectorMPC.cxx \
	Faun_Dict.C

libFaun_la_LDFLAGS =  -L$(ONLINE_MAIN)/lib -L$(OFFLINE_MAIN)/lib -L$(ROOTSYS)/lib -lpmonitor -fPIC

libFaun_la_LIBADD = \
  -lfun4allfuncs \
  -L$(libdir) \
  -L$(ONLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib \
  -L$(ROOTSYS)/lib \
 -lGpad -lHist -lGraf -lGraf3d -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lRIO -lNet -lThread -lCore -lCint -lm -ldl -lSpectrum

Faun_Dict.C: \
	FTMUnpacker.h \
	FTrace.h \
	FAPD.h \
	FMP.h \
	FDetectorMPC.h \
	FaunLinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*