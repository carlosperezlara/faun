libFaun.so:
	g++ -c -Wall -Werror -fpic FTrace.cxx 
	g++ -c -Wall -Werror -fpic FAPD.cxx 
	g++ -c -Wall -Werror -fpic FADC.cxx 
	g++ -c -Wall -Werror -fpic FMP.cxx 
	g++ -c -Wall -Werror -fpic FDetectorMPC.cxx `root-config --libs` -I`root-config --incdir`
	g++ -c -Wall -Werror -fpic FDetectorBC.cxx `root-config --libs` -I`root-config --incdir`
	g++ -c -Wall -Werror -fpic FDetectorEX.cxx `root-config --libs` -I`root-config --incdir`
	g++ -c -Wall -Werror -fpic FMaze.cxx  `root-config --libs` -I`root-config --incdir`
	g++ -shared -o libFaun.so FTrace.o FMP.o FAPD.o FDetectorMPC.o FMaze.o
	rm *.o

clean:
	rm *.o *.so
