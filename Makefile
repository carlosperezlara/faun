libFaun.so:
	g++ -c -Wall -Werror -fpic FTrace.cxx 
	g++ -c -Wall -Werror -fpic FAPD.cxx 
	g++ -c -Wall -Werror -fpic FMP.cxx 
	g++ -c -Wall -Werror -fpic FDetectorMPC.cxx 
	g++ -shared -o libFaun.so FTrace.o FMP.o FAPD.o FDetectorMPC.o
	rm *.o

clean:
	rm *.o *.so
