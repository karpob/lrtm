import TCM
from scipy.io import loadmat
import numpy
import pylab
#
#filez=['testTCM/Depleted_Ammonia_Clouds.mat','testTCM/Depleted_Water_Clouds.mat','testTCM/Enhanced_Ammonia_Clouds.mat','testTCM/Enhanced_Water_Clouds.mat','testTCM/Mean_Lindal_Clouds.mat']
#filez=['deBoer_check/Depleted_Ammonia_Clouds.mat','deBoer_check/Depleted_Water_Clouds.mat','deBoer_check/Enhanced_Ammonia_Clouds.mat','deBoer_check/Enhanced_Water_Clouds.mat','deBoer_check/Mean_Lindal_Clouds.mat']
filez=['Mean_Lindal_Clouds.mat']
for mat in filez:
	MatlabRuns=loadmat(mat)
	dz =1#.MatlabRuns['dz']
	case=MatlabRuns['case_select']
	
	XHe_i =MatlabRuns['XHe'][:,case-1]
	tcmeI=MatlabRuns['tcme'][:]
	idx=numpy.nonzero(tcmeI[:,0]==1.0)
	T_targ_i=165.#float(tcmeI[idx,1])
	#print "HEEEEEEEEEEEEEEEEEEY!", float(tcmeI[idx,1])
	XH2S_i =MatlabRuns['XH2S'][:,case-1]
	XNH3_i =MatlabRuns['XNH3'][:,case-1]
	XH2O_i =MatlabRuns['XH2O'][:,case-1]
	XCH4_i =MatlabRuns['XCH4'][:,case-1]
	XPH3_i =MatlabRuns['XPH3'][:,case-1]
	print XNH3_i,XH2O_i
	print MatlabRuns['XH2O'][:]
	print MatlabRuns['XNH3'][:]
	XCO =MatlabRuns['XCO'][:]
	P_temp =float(MatlabRuns['P_temp'][:])	
	T_temp =1653.9901123#float(MatlabRuns['T_temp'][:])
	g0_i =float(MatlabRuns['g0_i'][:])
	R0_i =float(MatlabRuns['R0e_i'][:])
	P0_i =float(MatlabRuns['P0_i'][:])
	#T_targ_i =float(MatlabRuns['T_targ_i'][:])
	P_targ_i =float(MatlabRuns['P_targ_i'][:])
	P_term_i =float(MatlabRuns['P_term_i'][:])
	use_lindal_in=1
	n_lindal_pts_i =61
	SuperSatSelf1_i=0.0
	SuperSatSelf2_i=0.0
	SuperSatSelf3_i=0.0
	SuperSatSelf4_i=0.0
	supersatNH3_i=0.0
	supersatH2S_i=0.0
	AutoStep_constant=int(MatlabRuns['AutoStep_constant'][:])
	fp =666.;#MatlabRuns['fp']
	dP_init =float(MatlabRuns['dP_init'][:])
	dP_fine =float(MatlabRuns['dP_fine'][:])
	P_fine_start=float(MatlabRuns['P_fine_start'][:])
	P_fine_stop=float(MatlabRuns['P_fine_stop'][:])
	use_dz =0#MatlabRuns['use_dz'][:]
	frain =float(MatlabRuns['frain'][:])
	select_ackerman=0#MatlabRuns['select_ackerman']
	Profiles=[]
	#data=MatlabRuns['Data'][:]

	print XHe_i
	print P0_i
	print P_temp,T_temp
	#pylab.plot(MatlabRuns['tcme'][:,3])
	#pylab.show()
	print 'stuff I need',dz
	print float(XHe_i)
	print float(XH2S_i)
        print float(XNH3_i)
        print float(XH2O_i)
        print float(XCH4_i)
        print float(XPH3_i)
        print float(XCO)
        print P_temp
        print T_temp
        print g0_i
        print R0_i
        print P0_i
        print T_targ_i
        print P_targ_i
        print P_term_i
        print use_lindal_in
        print n_lindal_pts_i
        print SuperSatSelf1_i                                             
        print SuperSatSelf2_i
        print SuperSatSelf3_i
        print SuperSatSelf4_i
        print supersatNH3_i
        print supersatH2S_i
        print AutoStep_constant
        print float(fp)
        print dP_init
        print dP_fine
        print P_fine_start
        print P_fine_stop
        print use_dz
        print frain
        print 'ackerman',select_ackerman

	nlevels=TCM.intoTheVoid( dz,
		  float(XHe_i),
		  float(XH2S_i),
                  float(XNH3_i),
                  float(XH2O_i),
                  float(XCH4_i),
                  float(XPH3_i),
                  float(XCO),
                  P_temp,
                  T_temp,
                  g0_i,
                  R0_i,
                  P0_i,
                  T_targ_i,
                  P_targ_i,
                  P_term_i,
                  use_lindal_in,
                  n_lindal_pts_i,
                  SuperSatSelf1_i,                                             
                  SuperSatSelf2_i,
                  SuperSatSelf3_i,
                  SuperSatSelf4_i,
                  supersatNH3_i,
                  supersatH2S_i,
                  AutoStep_constant,
                  float(fp),
                  dP_init,
                  dP_fine,
                  P_fine_start,
                  P_fine_stop,
                  use_dz,
                  frain,
                  select_ackerman)

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
