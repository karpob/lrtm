platform :=$(shell uname )
hardware :=$(shell uname -a )

ifeq ($(platform),Linux)
	# Linux
	CC=g++
	CCC=g++
	CCOPTS= -Wall -DLinux -DIntel 
	OS=linux
	IFLAGS = -I/usr/include/python2.7/ -I.
	LD=g++ *.o -shared -o _TCM.so
	ifeq ($(hardware),x86_64)
		CCOPTS= -Wall -DLinux -DIntel -m64 -fPIC 
	endif
endif

#ifeq ($(platform),CYGWIN_NT-5.1)
#	# Windows

#endif

ifeq ($(platform),Darwin)
	# Mac
	CCOPTS= -Wall -DDarwin -m64 -DPIC -arch x86_64
	IFLAGS = -I/sw/include/python2.7/
	LD=ld -arch x86_64 -bundle -flat_namespace -undefined suppress -o _TCM.so *.o
endif

TCM:
	g++ -c  *.c  $(CCOPTS) $(IFLAGS)
	$(LD)
	cp _TCM.so ..

clean:
	rm *.o
	rm *.so
