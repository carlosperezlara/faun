libFaun.so:
	g++ -c -Wall -Werror -fpic FTrace.cxx 
	g++ -c -Wall -Werror -fpic FAPD.cxx 
	g++ -c -Wall -Werror -fpic FMP.cxx 
	g++ -shared -o libFaun.so FTrace.o FMP.o FAPD.o 
	rm *.o

clean:
	rm *.o *.so
