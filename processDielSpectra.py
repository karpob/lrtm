import cPickle
from newQlib import NAprocess
import numpy,gc
gc.enable()
file_handle=open('allDielSpectra.pkl')
Spectra=cPickle.load(file_handle)
file_handle.close()


file_handle=open('allDielSpectra.pklS12.pkl')
Spectra2=cPickle.load(file_handle)
file_handle.close()

file_handle=open('allDielSpectra.pklS22.pkl')
S22=cPickle.load(file_handle)
file_handle.close()


file_handle=open('allDielSpectraAxis.pkl')
Axis=cPickle.load(file_handle)
file_handle.close()
file_handle=open('allDielS11.pkl')
S11d=cPickle.load(file_handle)
file_handle.close()

file_handle=open('allGasS11.pkl')
S11g=cPickle.load(file_handle)
file_handle.close()

file_handle=open('allTransSpectra.pkl')
Trans=cPickle.load(file_handle)
file_handle.close()

file_handle=open('allTransSpectra.pklS12.pkl')
Trans2=cPickle.load(file_handle)
file_handle.close()



file_handle=open('all_data.pkl')

#processed frequency objects
pfo={}
jj=0
jj=0
allExperiments={}
Experiment={}
press=[]
#print Spectra.keys()

for i in Spectra.keys():
        if(i<20):
                s=Spectra[i]
                ax=Axis[i]
                ssd=S11d[i]
                ssg=S11g[i]
                tTrans=Trans[i]
        else:
                s=Spectra2[i]
                ax=Axis[i]
                ssd=S22[i]
                ssg=S22[i]
                tTrans=Trans2[i]
                        
        ex=cPickle.load(file_handle)
        InterpSpectra=ex['ArInterpSpec']
        for pressureLevel in s.keys():
                
                for resonance in s[pressureLevel]:
                        
                        pfo=NAprocess(resonance,10**(ssd[pressureLevel][jj]/20),ax[pressureLevel][jj],numpy.mean(tTrans[pressureLevel][jj],axis=1).reshape([1601,1]),jj,i,pressureLevel,InterpSpectra[jj,pressureLevel-1,:,:],ex['ArPeakFa'][jj,pressureLevel-1],ex['ArQ'][jj,pressureLevel-1],'Diel')
                        press.append(pfo)
                        jj+=1
                jj=0
                Experiment[pressureLevel]=press
                press=[]
        allExperiments[i]=Experiment
        Experiment={}
        Spectra[i]=[]
file_handle.close()                
file_handle2=open('allDielQs.pkl','w')
cPickle.dump(allExperiments,file_handle2)
file_handle2.close()     
                                   
