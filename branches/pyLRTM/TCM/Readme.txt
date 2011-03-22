

If you change the interface to the C code, you will need to run swig.
cd src
$swig -python TCM.i

An then make
$make
and then

$cp TCM.py ..

At the core of what the makefile does...
g++ -c -fpic *.c -I/tgrs05/scratch/local/include/python2.6/
g++ -shared *.o /tgrs05/scratch/local/lib/python2.6/config/libpython2.6.a -o _TCM.so




