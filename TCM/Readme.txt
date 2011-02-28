swig -python TCM.i
g++ -c -fpic *.c -I/tgrs05/scratch/local/include/python2.6/
g++ -shared *.o /tgrs05/scratch/local/lib/python2.6/config/libpython2.6.a -o _TCM.so

python test_tcm.py

I think we still need to worry about how the "inner python" is called to load compressibility.
A makefile would be nice too...

