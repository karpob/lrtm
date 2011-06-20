import TCM
from scipy.io import loadmat
from scipy.io import netcdf
import numpy
import pylab

filez=['Depleted_Ammonia_Clouds.mat',
       'Depleted_Water_Clouds.mat',
       'Enhanced_Ammonia_Clouds.mat',
       'Enhanced_Water_Clouds.mat',
       'Mean_Lindal_Clouds.mat']

for mat in filez:
	MatlabRuns=loadmat(mat)
	dz =0#MatlabRuns['dz']
	case=MatlabRuns['case_select']
	print 'mr. case?',str(case[0])
	XHe_i =MatlabRuns['XHe'][:,case-1]
	tcmeI=MatlabRuns['tcme'][:]
	idx=numpy.nonzero(tcmeI[:,0]==1.0)
	XH2S_i =MatlabRuns['XH2S'][:,case-1]
	XNH3_i =MatlabRuns['XNH3'][:,case-1]
	XH2O_i =MatlabRuns['XH2O'][:,case-1]
	XCH4_i =MatlabRuns['XCH4'][:,case-1]
	XPH3_i =MatlabRuns['XPH3'][:,case-1]
	XCO =MatlabRuns['XCO'][:]
	P_temp =float(MatlabRuns['P_temp'][:])	
	T_temp =float(MatlabRuns['T_temp'][:])
	g0_i =float(MatlabRuns['g0_i'][:])
	R0_i =float(MatlabRuns['R0e_i'][:])
	P0_i =float(MatlabRuns['P0_i'][:])
	T_targ_i =float(MatlabRuns['T_targ_i'][:])
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
	fp =MatlabRuns['fp']
	dP_init =0.25#float(MatlabRuns['dP_init'][:])
	dP_fine =0.25#float(MatlabRuns['dP_fine'][:])
	P_fine_start=float(MatlabRuns['P_fine_start'][:])
	P_fine_stop=float(MatlabRuns['P_fine_stop'][:])
	use_dz =0#MatlabRuns['use_dz'][:]
	frain =float(MatlabRuns['frain'][:])
	select_ackerman=0#MatlabRuns['select_ackerman']
	Profiles=[]
	print "XHe, XH2S, XNH3, XH2O, XCH4, XPH3, XCO"
	print float(XHe_i),float(XH2S_i),float(XNH3_i),float(XH2O_i),float(XCH4_i),float(XPH3_i),float(XCO)
        

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

	data=MatlabRuns['tcme']
	layerKeys=['P','T','z','XH2','XHe','XH2S','XNH3','XH2O','XCH4','XPH3','DNH4SH','DH2S','DNH3','DH2O','DCH4','DPH3','DSOL','g','mu','DSOL_NH3','P_real']
	layer={}
	for key in layerKeys:layer[key]=[]			

	for j in range(0, nlevels+1):
		for i in range(0,len(layerKeys)):
			layer[layerKeys[i]].append(TCM.getFloatValues(i,j))
					
	f=netcdf.netcdf_file(mat+'.nc','w')
	f.history='created for '+str(case[0])+'case.'
	for key in layerkeys:
	        f.createDimension(key,'f',len(layer[key]))
        f.close()
        	
