from newQlib import NAprocess
from scipy.io import loadmat
import numpy
a=loadmat('1_vac1.mat')
B=a['S21']
freq=a['faxis']
tran=numpy.zeros(len(freq))
numswp=30
poo=NAprocess(B,freq,tran)
print poo
