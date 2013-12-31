import cPickle
from newQlib import NAprocess
import numpy
file_handle=open('allGasSpectra.pkl')
Spectra=cPickle.load(file_handle)
file_handle.close()

file_handle=open('allGasSpectra.pklS12.pkl')
Spectra2=cPickle.load(file_handle)
file_handle.close()

file_handle=open('allGasSpectraAxis.pkl')
Axis=cPickle.load(file_handle)
file_handle.close()

file_handle=open('allGasS11.pkl')
S11=cPickle.load(file_handle)
file_handle.close()

file_handle=open('allGasSpectra.pklS22.pkl')
S22=cPickle.load(file_handle)
file_handle.close()

file_handle=open('allTransSpectra.pkl')
Trans=cPickle.load(file_handle)
file_handle.close()

file_handle=open('allTransSpectra.pklS12.pkl')
Trans2=cPickle.load(file_handle)
file_handle.close()

file_handle=open('all_data.pkl')
import gc
gc.enable()
#processed frequency objects
pfo={}
jj=0
jj=0
allExperiments={}
Experiment={}
press=[]
plevels=[1,]


for i in Spectra.keys():
        if(i<20):
                s=Spectra[i]
                ax=Axis[i]
                ss=S11[i]
                tt=Trans[i]
        else:
                s=Spectra2[i]
                ax=Axis[i]
                ss=S22[i]
                tt=Trans2[i]         
    ex=cPickle.load(file_handle)
    InterpSpectra=ex['nh3InterpSpec']
    for kk in ex.keys():
           if('nh3' in kk): print kk        
        for pressureLevel in s.keys():
                        
                for resonance in s[pressureLevel]:
                        print i, pressureLevel,resonance.shape,InterpSpectra.shape,ex['nh3Q'].shape,ex['f'].shape
                    print jj,pressureLevel
                        gamma=10**(ss[pressureLevel][jj]/20)
                        pfo=NAprocess(resonance,gamma,ax[pressureLevel][jj],numpy.mean(tt[pressureLevel][jj],axis=1).reshape([1601,1]),jj,i,pressureLevel,InterpSpectra[jj,pressureLevel-1,:,:],ex['nh3PeakFa'][jj,pressureLevel-1],ex['nh3Q'][jj,pressureLevel-1],'Gas')
                        press.append(pfo)
                        jj+=1
                jj=0
                Experiment[pressureLevel]=press
                press=[]
        allExperiments[i]=Experiment
        Experiment={}
        Spectra[i]=[]
        S11[i]=[]
        Trans[i]=[]
        Axis[i]=[]
        ex=[]
        gc.collect()
                
file_handle.close()    
file_handle2=open('allGasQs.pkl','w')
cPickle.dump(allExperiments,file_handle2)

file_handle2.close()                                        
