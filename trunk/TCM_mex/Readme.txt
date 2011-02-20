Back in my day we passed things on in temporary text files: welcome to 2010, we don't do that anymore!

Before running program with TCM, you now must compile using:
/Applications/MATLAB_R2010a.app/bin/mex src/*.C -I/opt/local/include/python2.5/ /opt/local/lib/python2.5/config/libpython2.5.a -o TCM

You will need to vary this according to your flavor of *nix. Windows...eh sorry I don't have a machine that will do that...probably similar?

-I is "include" search path (the path to the Python.h file).

Depending on your distribution you'll need to change site-packages location currently in the C source in LAYERS.C
to_site_packages="/opt/local/lib/python2.5/site-packages";

In future versions pass this on as an arg from the matlab script?
In summary you need 3 things:
1. Path to your Python.h file
2. path+filename of your libpythonX.Xa file (where X.X is your version of python)
3. path of your python site-packages.

Dependencies:
1. numpy
2. matplotlib
3. scipy

At somepoint I should remove matplotlib deps (I'm not sure if I cleaned them all out to by pure numpy deps)



