import cPickle
from newQlib import calculateAbsorptionAndErrs
import numpy
import pylab
file_handle=open('allGasQs.pkl')
gasOb=cPickle.load(file_handle)
file_handle.close()

file_handle=open('allDielQs.pkl')
dielOb=cPickle.load(file_handle)
file_handle.close()

file_handle=open('allVacQs.pkl')
vacOb=cPickle.load(file_handle)
file_handle.close()

file_handle=open('allGasVSWR.pkl')
vswrGas=cPickle.load(file_handle)
file_handle.close()   

file_handle=open('allDielVSWR.pkl')
vswrDiel=cPickle.load(file_handle)
file_handle.close()

file_handle=open('allGasTimes.pkl')
timez=cPickle.load(file_handle)
file_handle.close()

allExperiments={}
press=[]
Experiment={}
aVals=[]
fVals=[]
for ExperimentNum in [9,]:#gasOb.keys():#experiment number
        for p in gasOb[ExperimentNum].keys():# pressure
                for resonance in range(0,12):
                        #print 'yes.'
                        #print ExperimentNum,p,resonance, gasOb[ExperimentNum][p][resonance]['Qe'][0],dielOb[ExperimentNum][p][resonance]['Qe'][0],timez[ExperimentNum][p][resonance]#,vacOb[ExperimentNum][resonance][:]
                        
                        abObs=calculateAbsorptionAndErrs(gasOb[ExperimentNum][p][resonance],dielOb[ExperimentNum][p][resonance],vacOb[ExperimentNum][1][resonance],timez[ExperimentNum][p][resonance],vswrGas[ExperimentNum][p][resonance],vswrDiel[ExperimentNum][p][resonance])
                        
                        if(ExperimentNum==9):
                                if(p==1):
                                        
                                        #print gasOb[ExperimentNum][p][resonance]['Qa'][0],dielOb[ExperimentNum][p][resonance]['Qa'][0]#gasOb[ExperimentNum][p][resonance]['Levela'],dielOb[ExperimentNum][p][resonance]['Levela']
                                        print abObs['AbsorptionK'],gasOb[ExperimentNum][p][resonance]['PeakFa']#,gasOb[ExperimentNum][p][resonance]['Qe'][0],gasOb[ExperimentNum][p][resonance]['PeakFe'],vacOb[ExperimentNum][1][resonance]['Qe'][0],vacOb[ExperimentNum][1][resonance]['PeakFe']
                     #vswrDiel[ExperimentNum][p][resonance],vswrGas[ExperimentNum][p][resonance]
                    aVals.append(abObs['AbsorptionK'])
                    fVals.append(gasOb[ExperimentNum][p][resonance]['PeakFa'])
                        press.append(abObs)
                Experiment[p]=press
                press=[]
        allExperiments[ExperimentNum]=Experiment
        
        Experiment={}
pylab.plot(numpy.asarray(fVals)/1.0e9,numpy.asarray(aVals),'x')
print aVals
pylab.semilogy()
pylab.ylim([1e-3,1e1])
pylab.show()
for a in aVals:
        print a
f=open('allProcessedAbsorption','w')
cPickle.dump(allExperiments,f)
f.close()