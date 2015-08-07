libFaun.so:
	g++ -c -Wall -Werror -fpic FTrace.cxx 
	g++ -c -Wall -Werror -fpic FAPD.cxx 
	g++ -c -Wall -Werror -fpic FMinipad.cxx 
	g++ -shared -o libFaun.so FTrace.o FAPD.o FMinipad.o


