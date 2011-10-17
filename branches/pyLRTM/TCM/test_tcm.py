from pyTCM.TCM import runTCM
from scipy.io import loadmat
import numpy
import pylab
filez=['Depleted_Water_Clouds.mat']
for mat in filez:
        IN={}
	MatlabRuns=loadmat(mat)
	print MatlabRuns.keys()
	dz =1#.MatlabRuns['dz']
	case=MatlabRuns['case_select']
	IN['XHe']=float(MatlabRuns['XHe'][:,case-1])
	tcmeI=MatlabRuns['tcme'][:]
	idx=numpy.nonzero(tcmeI[:,0]==1.0)
	IN['T_targ_i']=165.#float(tcmeI[idx,1])
	
	IN['XH2S'] =float(MatlabRuns['XH2S'][:,case-1])
	IN['XNH3'] =float(MatlabRuns['XNH3'][:,case-1])
	IN['XH2O'] =float(MatlabRuns['XH2O'][:,case-1])
	IN['XCH4'] =float(MatlabRuns['XCH4'][:,case-1])
	IN['XPH3'] =float(MatlabRuns['XPH3'][:,case-1])
	
	IN['XCO'] =float(MatlabRuns['XCO'][:])
	IN['P_temp'] =float(MatlabRuns['P_temp'][:])	
	IN['T_temp'] =1653.9901123
	IN['g0'] =float(MatlabRuns['g0_i'][:])
	IN['R0'] =float(MatlabRuns['R0e_i'][:])
	IN['P0'] =float(MatlabRuns['P0_i'][:])
        IN['P_targ'] =float(MatlabRuns['P_targ_i'][:])
	IN['P_term'] =float(MatlabRuns['P_term_i'][:])
	IN['use_lindal']='Y'
	
	IN['SuperSatSelf1']=0.0
	IN['SuperSatSelf2']=0.0
	IN['SuperSatSelf3']=0.0
	IN['SuperSatSelf4']=0.0
	IN['supersatNH3']=0.0
	IN['supersatH2S']=0.0
	IN['AutoStep']=False
	IN['AutoStep_constant']=int(MatlabRuns['AutoStep_constant'][:])
	IN['fp']=666.#MatlabRuns['fp']
	IN['dP_init'] =float(MatlabRuns['dP_init'][:])
	IN['dP_fine'] =float(MatlabRuns['dP_fine'][:])
	IN['P_fine_start']=float(MatlabRuns['P_fine_start'][:])
	IN['P_fine_stop']=float(MatlabRuns['P_fine_stop'][:])
	IN['use_dz'] =False#MatlabRuns['use_dz'][:]
	IN['frain'] =float(MatlabRuns['frain'][:])
	IN['select_ackerman']=0#MatlabRuns['select_ackerman']
	Profiles=[]
	#data=MatlabRuns['Data'][:]

	

	layers,others=runTCM(IN)
	for item in layers['P']:
	        print item
	        
	#print 'll\n',layers
	#print 'oo\n',others
        """
	#data=numpy.zeros([nlevels,22])
	data=MatlabRuns['tcme']
	layerKeys=['P','T','z','XH2','XHe','XH2S','XNH3','XH2O','XCH4','XPH3','DNH4SH','DH2S','DNH3','DH2O','DCH4','DPH3','DSOL','g','mu','DSOL_NH3','P_real']
	layer={}
	for key in layerKeys:layer[key]=[]			

	for j in range(0, nlevels+1):
		for i in range(0,len(layerKeys)):
			layer[layerKeys[i]].append(TCM.getFloatValues(i,j))
					
	#Profiles.append(layer)

	offsets=[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0];

	offset=0
	
	#print 'data,Profiles',numpy.shape(data),numpy.shape(numpy.asarray(Profiles[0]))
	
	for k in range(0,len(layerKeys)):#len(layerKeys)):
		offset=offset+offsets[k]
		print numpy.shape(layer[layerKeys[k]])	
		#if(k==20): 
		#	print Profiles[i][layerKeys[k]][0:nlevels-1]
		#	print data[i,1][1:nlevels,k+offset]
		print 'k offset',type(k),type(offset)
		#print k,layerKeys[k],k+offset,data[:,k+offset],layer[layerKeys[k]]#numpy.max(numpy.max(data[1:nlevels,k+offset]-data1[0:nlevels-1,k])),numpy.min(numpy.min(data[1:nlevels,k+offset]-Profiles[0:nlevels-1,k]))
		if(layerKeys[k]=='T'):
			pylab.figure()
			#
			#print numpy.shape(layer[layerKeys[k]])
		        print layerKeys[k],data[0,k+offset],layer[layerKeys[k]][0]
			#pylab.plot(data[:,k+offset]-layer[layerKeys[k]],'b')
			
			#pylab.plot(layer[layerKeys[k]],'k')
			#print numpy.shape(data[:,k])
			#print numpy.shape(layer[layerKeys[k]])
			
			pylab.show()
		if(layerKeys[k]=='DSOL'):
			pylab.figure()
		#	#data[:,k]-
			pylab.plot(layer[layerKeys[k]],'r')
			pylab.plot(data[:,k+offset],'b')
		#	#pylab.plot(data[0:194,k]-layer[layerKeys[k]],'r')
		#	#pylab.plot(layer[layerKeys[k]],'g')
			pylab.savefig('poop.pdf')
		if(layerKeys[k]=='P'):
			print k,layerKeys[k],data[0,k+offset],layer[layerKeys[k]][0]
		if(layerKeys[k]=='P_real'):

			print k,layerKeys[k],data[0,k+offset],layer[layerKeys[k]][0]
		#	pylab.figure()
		#	pylab.plot(data[:,k]-layer[layerKeys[k]],'g')
		#	#pylab.plot(data[0:194,k]-layer[layerKeys[k]],'g')
		#	pylab.show()
	offset=0	
	i+=1
	"""
