platform :=$(shell uname)
hardware :=$(shell uname -i)

ifeq ($(platform),Linux)
	# Linux
	CC=g++
	CCC=g++
	CCOPTS= -Wall -DLinux -DIntel 
	OS=linux
	IFLAGS = -I/usr/include/python2.7/ -I.

	ifeq ($(hardware),x86_64)
		IFLAGS = -I/usr/include/python2.7/Python.h -I. 
	endif
endif

#ifeq ($(platform),CYGWIN_NT-5.1)
#	# Windows
#
#endif

#ifeq ($(platform),Darwin)
#	# Mac
#endif

TCM:
	g++ -c  *.c  $(CCOPTS) $(IFLAGS)
	g++ -shared *.o -o _TCM.so
	cp _TCM.so ..

clean:
	rm *.o
	rm *.so
