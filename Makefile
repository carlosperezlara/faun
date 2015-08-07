libFaun.so:
	g++ -c -Wall -Werror -fpic FTrace.cxx
	g++ -shared -o libFaun.so FTrace.o
