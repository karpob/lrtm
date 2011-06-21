platform :=$(shell uname)
hardware :=$(shell uname )

ifeq ($(platform),Linux)
	# Linux
	CC=g++
	CCC=g++
	CCOPTS= -Wall -DLinux -DIntel 
	OS=linux
	IFLAGS = -I/usr/include/python2.7/ -I.

	ifeq ($(hardware),x86_64)
		IFLAGS = -I/usr/include/python2.7/ -I. 
	endif
endif

#ifeq ($(platform),CYGWIN_NT-5.1)
#	# Windows
#
#endif

ifeq ($(platform),Darwin)
	# Mac
	CCOPTS= -Wall -DDarwin -m64 -DPIC -arch x86_64
	IFLAGS = -I/Users/karpob2/py/include/python2.7/
	#LFLAGS = /sw/lib/python2.7/config/libpython2.7.a
endif

TCM:
	g++ -c  *.c  $(CCOPTS) $(IFLAGS)
	ld -arch x86_64 -bundle -flat_namespace -undefined suppress -o _TCM.so *.o  
	cp _TCM.so ..

clean:
	rm *.o
	rm *.so
