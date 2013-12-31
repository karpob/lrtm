import cPickle
from newQlib import NAprocess,processVSWR
import numpy
file_handle=open('allGasSpectra.pkl')
Spectra=cPickle.load(file_handle)
file_handle.close()


file_handle=open('allGasSpectraAxis.pkl')
Axis=cPickle.load(file_handle)
file_handle.close()

file_handle=open('allGasS11.pkl')
S11=cPickle.load(file_handle)
file_handle.close()


#processed frequency objects
pfo={}
jj=0
jj=0
allExperiments={}
Experiment={}
press=[]
for i in Spectra.keys():
        s=Spectra[i]
        ax=Axis[i]
        ss=S11[i]
        for pressureLevel in s.keys():
                
                for resonance in s[pressureLevel]:
                        print i, pressureLevel,jj
                        pfo=processVSWR(ss[pressureLevel][jj],resonance)#NAprocess(resonance,ax[pressureLevel][jj],numpy.zeros(1601))
                        press.append(pfo)
                        jj+=1
                jj=0
                Experiment[pressureLevel]=press
                press=[]
        allExperiments[i]=Experiment
        Experiment={}        

file_handle=open('allGasVSWR.pkl','w')
cPickle.dump(allExperiments,file_handle)
file_handle.close()                                        
