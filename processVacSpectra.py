import cPickle
from newQlib import NAprocess
import numpy
file_handle=open('allVacSpectra.pkl')
Spectra=cPickle.load(file_handle)
file_handle.close()


file_handle=open('allVacAxis.pkl')
Axis=cPickle.load(file_handle)
file_handle.close()

file_handle=open('allVacS11.pkl')
S11=cPickle.load(file_handle)
file_handle.close()

file_handle=open('allTransSpectra.pkl')
Trans=cPickle.load(file_handle)
file_handle.close()

jj=0
allExperiments={}
Experiment={}
press=[]
EE=[15,]
PP=[3,]

for i in Spectra.keys():
        s=Spectra[i]
        ax=Axis[i]
        ss=S11[i]
	tt=Trans[i]
        for pressureLevel in s.keys():
                
                for resonance in s[pressureLevel]:
                        print i, pressureLevel,jj
			gamma=10**(ss[pressureLevel][jj]/20)
                        #pfo=NAprocess(resonance,0.0*ss[pressureLevel][jj],ax[pressureLevel][jj],numpy.zeros(1601),jj)
			pfo=NAprocess(resonance,0.0*gamma,ax[pressureLevel][jj],numpy.mean(tt[pressureLevel][jj],axis=1).reshape([1601,1]),jj)
                        press.append(pfo)
                        jj+=1
                jj=0
                Experiment[pressureLevel]=press
                press=[]
        allExperiments[i]=Experiment
        Experiment={}
file_handle=open('allVacQs.pkl','w')
cPickle.dump(allExperiments,file_handle)
file_handle.close()                                        
