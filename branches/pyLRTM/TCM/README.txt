This code generates Jovian profiles presented in Karpowicz and Steffes, 2012.

To run the code one must have numpy, scipy, and matlplotlib installed.

To run the code (under linux, or mac) in a terminal simply change to this 
directory and type:

$python jupiterCasesTCMonly.py

Under windows one can simply double click on the jupiterCasesTCMonly.py script.


The code will run through all cases presented in Karpowicz, and Steffes, 2012.
The following plots will be generated:
CloudsLindal.pdf
CloudsSeiff.pdf
CloudsSeiff2.pdf
ConstituentsLindal.pdf
ConstituentsSeiff.pdf
Residual.pdf
TP.pdf

Code which is used by jupiterCasesTCMonly.py includes functions in the pyTCM 
directory, python_compressibility direcotry, along with DeBoerTCM.py. 

The source code used to compute $c_p$, and $P$ can be obtained using 
python_compressibility/calc_Cp/getCp.py, and python\_compressibility/calc_Cp/getPreal.py.


All other directories are legacy code from the original C implementation
of DeBoer, 1995's C TCM model.
  
