def maintamone(Raydirection,Rayorigin,tcme,tcmp,ao,bo,co,f,no_ph3,
        select_h2h2_model,select_ammonia_model,select_water_model,
        include_clouds,N_ring_one,Nphi,BWHM,refractivity_source,
        cassini_pattern,cassini_data_path):
        import numpy,os
        from findrefindex import findrefindex
        from findellipseradiusvector import findellipseradiusvector
        from findraypath import findraypath
        from findkappa import findkappa
        from ftau import ftau
        from ftam import ftam
        from fwght import fwght
        from beamsample import beamsample
        #from findrefindex import findrefindex
	#  maintamone.py -> main routine for LRTM function calculates Tb, and
	#                  weighting functions
	#
	#                           Input:
	#                                 -->Raydirection: Direction of Beam in
	#                                                  terms of unit vector [X Y Z]
	#                                 -->Rayorigin: Position of spacecraft in [X Y Z]
	#                                 -->tcme: Thermochemical model values for Equatorial region (values not squished in altitude (z))
	#                                 -->tcmp: Themrochemical model values for Polar region (values squished in altitude (z))
	#                                 -->ao: Ellipse (cm) value along X
	#                                 -->bo: Ellipse (cm) value along Y
	#                                 -->co: Ellipse (cm) value along Z
	#                                 -->f: Frequency in GHz
	#                                 -->no_ph3: Flag to include/exclude Hoffman Phosphine decay module 
	#                                           no_PH3=1-->active will decay phosphine
	#                                           no_PH3=0-->inactive, will use PH3 profile provided by TCM
	#                                 -->select_h2h2_model: Select model for collisionally induced absorption by H2
	#                                 -->select_ammonia_model: Select model for ammonia absorption
	#                                 -->select_water_model: Select model for water absorption
	#                                 -->include_clouds: Include or exclude absorption by clouds
	#                                 -->N_ring_one: Number of rays in first ring of beam battern (sampling param)
	#                                 -->Nphi: Number of phi rings (sampling param)
	#                                 -->BWHM: Beamwidth half maximum in degrees (3dB beamwidth of antenna)
	#                                 -->refractivity_source: Select between original refractivity, and all inclusive refractivity options                     
	#
	#                            
                                      
	#global CRITICALFLAG
	CRITICALFLAG=0				# If is '1' then critical refraction reached for that ray
	USEBEAM=1

	# Flip order from low atltitude to high over to top-down
	recordlength1=numpy.max(tcme.shape)
	recordlength2=numpy.max(tcmp.shape)

	if(recordlength1>recordlength2):
    		recordlength=recordlength2

	if(recordlength2>recordlength1):
    		recordlength=recordlength1

	if(recordlength1==recordlength2):
    		recordlength=recordlength1

	
	#
	# Extract Prameters, flip around
	# Add a zero for space (as in No atmosphere P=0)
	P=numpy.r_[0.,numpy.flipud(tcme[:,0])]
	T=numpy.r_[2.7,numpy.flipud(tcme[:,1])]
	# third is the dR vector -handled elsewhere since need two
	major=1.e5*numpy.flipud(tcme[:,2]) # convert from km to cm
	minor=1.e5*numpy.flipud(tcmp[:,2]) # convert from km to cm
	xH2=numpy.r_[0.,numpy.flipud(tcme[:,3])]
	xHe=numpy.r_[0.,numpy.flipud(tcme[:,4])]
	xH2S=numpy.r_[0.,numpy.flipud(tcme[:,5])]
	xNH3=numpy.r_[0.,numpy.flipud(tcme[:,6])]
	xH2O=numpy.r_[0.,numpy.flipud(tcme[:,7])]
	xCH4=numpy.r_[0.,numpy.flipud(tcme[:,8])]
	xPH3=numpy.r_[0.,numpy.flipud(tcme[:,9])]
	# Get cloud Densities from TCME
	DNH4SH=numpy.r_[0.,numpy.flipud(tcme[:,11])]
	DH2S=numpy.r_[0.,numpy.flipud(tcme[:,12])]
	DNH3=numpy.r_[0.,numpy.flipud(tcme[:,13])]
	DH2O=numpy.r_[0.,numpy.flipud(tcme[:,14])]
	DSOL=numpy.r_[0.,numpy.flipud(tcme[:,17])]

	# cause ph3 decay (photolyse)
	if (no_ph3>0.):
    		look_out_for_decaying_ph3=1
    		print 'Computing Phosphine Decay.'
    		xPH3[1]=xPH3[1]
    		xPH3(recordlength)
    		decay=ph3decay(P,xPH3)
    		xPH3=decay


	smallestdr=min(len(major),len(minor))
	minor=minor[0:smallestdr]		# unify lengths
	major=major[0:smallestdr]    # unify lengths

	# Convert mole fraction to partial pressures
	P_H2=P*xH2
	P_He=P*xHe
	P_H2S=P*xH2S
	P_NH3=P*xNH3
	P_H2O=P*xH2O
	P_CH4=P*xCH4
	P_PH3=P*xPH3

	sph3=P_PH3.shape[0]

	# Find the refractive index profile
	refindex=findrefindex(T,P_H2,P_He,P_CH4,P_H2O,DSOL,f,refractivity_source)

	# Find the elliptical shells used by 'findraypath' (major/minor) are vectors of how radius changes with each index
	ellipses=findellipseradiusvector(ao,bo,co,major,minor)

	# IMAGING (mini)


	# DOES BORESIGHT MISS-IF SO=> END
	[intercept,internormal,d,t,masterindexa,missflag]=findraypath(recordlength,refindex,P,ellipses,Rayorigin,Raydirection)

	if missflag==1:
   		disp('Boresight Misses Planet')
   		Tbeam=2.7
   		weightingfactor=0
   		return	[Tbeam,zenith,wfwa,refindex,Tatma,intercepts_boresight,intercepts_b]			# This will help with full imaging (if this file becomes a function
   # that gets passed various raydirection (boresights)


	###########################################################################
	#  Calculate Zenith Angle
	#
	###########################################################################
	# what Hoffman calls zenith is the angle between normal to the surface,
	# and observer. Not what you'd think when he talks about nadir, and
	# zenith. This is typically what one thinks of from a "ground-based"
	# observation.
	###########################################################################
	zenith=numpy.arccos(numpy.dot(-internormal[0,:],Raydirection))*(180./numpy.pi)
        print 'zenith angle',zenith
	#Calculate absorption coeff as a function of TP 
	kappa=findkappa(f,T,P,P_H2,P_He,P_NH3,P_H2O,P_CH4,P_PH3,P_H2S,
                xH2,xHe,xNH3,xH2O,DNH4SH,DH2S,DNH3,DH2O,DSOL,
                select_h2h2_model,select_ammonia_model,select_water_model,include_clouds)
	print 'Absorption coefficients computed.'

	# Calculate Ta (Ray brightness temperature at boresight)
	if CRITICALFLAG==1:
   		print 'Critical Refraction'
   		disa=d
   		bs=1
   		[tau_a,tau]=ftau(kappa,d,masterindexa)	# Comes out with deepest first
   		[Tatma,wlayersa]=ftam(T,tau,tau_a,masterindexa)
   		sura=intercept[0,:]
   		intercepts_boresight=intercept
   		CRITICALFLAG=0
	else:
   		[tau_a,tau]=ftau(kappa,d,masterindexa)	# Comes out with deepest first
   		disa=d
   		bs=1
   		[Tatma,wlayersa]=ftam(T,tau,tau_a,masterindexa)
   		sura=intercept[0,:]
   		intercepts_boresight=intercept
   

	# BORESIGHT DIDN'T MISS

	# d is the pathlength, t is theta, masterindex are the indices of the raypaths
	# Need to gets tau's
	# tau_a is the tau of that layer, tau is the cumulative summations of those layers
	# Now do the beamspread and rotate beampattern towards planet along look
	# vector
	if(cassini_pattern==1):
    		[beamz,beam_weightz,beam_sum]=load_cassini_beampattern(cassini_data_path)
	else:
    		[beamz,beam_weightz,beam_sum]=beamsample(Nphi,N_ring_one,BWHM)


	[Vr1,Zr]=rotbeam(Raydirection,beamz)

	# initialize wlayers
	wlayersb=numpy.zeros(len(P)-1,len(Vr1))

	# masterindex for wght fnctns
	windexb=numpy.zeros(len(P)-1,len(Vr1))
	missb=0

	# Calculate Optical Depth values, brightness temperature and weighting
	# function along each ray path.
	#VRONE=0
	print 'Raypaths calculated, calculating brightness temperatures along rays.'
	for p in range(0,len(Vr1)):
   		bs=bs+1						# bs-beamspread- keeps track of beamspread samples
   		Rd=numpy.transpose(Vr1[:,p])
   		[intercept,internormal,d,t,masterindexb,missflag]=findraypath(recordlength,refindex,P,ellipses,Rayorigin,Rd)
   		windexb[0:masterindexb.shape[0],p]=masterindexb
   		
   		disb=d
   		if missflag==1:
      			print 'This beam-sample  (b) Misses Planet.'
      			Tatmb[p]=2.7
      			wlayersb[:,p]=0.
      			missb=missb+1.
      			surb[p,:]=numpy.array([[0., 0., 0.]])
   		elif CRITICALFLAG==1:
      			print 'Critical Refraction.'
      			[tau_a,tau]=ftau(kappa,d,masterindexb)	# Comes out with deepest first
      			[Tatmb[p],wtemp]=ftam(T,tau,tau_a,masterindexb)
      			wlayersb[0:wtemp.shape[0],p]=wtemp
      			surb[p,:]=intercept[0,:]
      			CRITICALFLAG=0
   		else:
      			[tau_a,tau]=ftau(kappa,d,masterindexb)	# Comes out with deepest first
      			[Tatmb[p],wlayersb[:,p]]=ftam(T,tau,tau_a,masterindexb)
      			surb[p,:]=intercept[0,:]
      			intercepts_b[p,:,:]=intercept
      
   


	# Apply Beam weights (Beam coupling) from beamsample.m
	print 'Applying beamweights.'
	Na=len(Tatma)
	Nb=len(Tatmb)

	Twa=numpy.sum(Tatma)*1./Na
	Twb=numpy.sum(Tatmb*beam_weightz)
	Tbeam=(Twa+Twb)/(1.+beam_sum)

	# Find weighting function
	wfa=1.*fwght(wlayersa,masterindexa)
	#wfb=beam_weightz.*fwght(wlayersb,windexb)

	# To make Matlab happy-make all the weight vectors the same size by padding with zeros
	# to the size of the Pressure vector
	sP=P.shape[0]
	sa=numpy.zeros([sP-wfa.shape[0],1])
	#sb=zeros(sP-size(wfb,1),1)


	# wfwa is already only 1 column
	wfwa=[wfasa]
	# turn into columns
	#wfwb=[wfbsb]
	return Tbeam,zenith,wfwa,refindex,Tatma,intercepts_boresight,intercepts_b
