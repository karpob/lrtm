def HanleySteffes(f,T,P,H2mr,Hemr,NH3mr):
	"""
	Final Hanley-Steffes 

	Input:
	        -->f: Frequency in GHz (vector)
	        --> T: Temperature in Kelvin, 
	        --> P: Pressure (bars)
	        -->H2mr: Hydrogen mole fraction
	        -->Hemr: Helium mole fraction
	        -->NH3mr: Ammonia mole fraction
	
	Output:
	        <-alphanh3: Absorption coefficient in dB/km
	
	        Opacity (alphanh3) is returned as a vector matching f
	        in dB/km.
	"""
	import numpy
	from scipy.io import loadmat
	#import pupynere
	#at some point convert database to netCDF..open format would be better.
	LineParameters=loadmat('rtm/nh3/nh3LineParameters/nh3lincat190Latest.mat')
	Eo=LineParameters['Eo']
	H2HeBroad=LineParameters['H2HeBroad']
	Io=LineParameters['Io']
	J=LineParameters['J']
	K=LineParameters['K']
	Kscale=LineParameters['Kscale']
	fo=LineParameters['fo']
	gammaNH3o=LineParameters['gammaNH3o']

	#load nh3lincat190Latest2	
	#Loads the Poynter Pickett (JPL) nh3 line catalog.  The workspace contains the arrays for Freq (fo) in GHz, intensity (Io) in inverse cm,
	#lower state energy (Eo) in inverse cm and self broadening due to NH3 (gammaNH3o) in GHz/bar.
	#Constants
	GHztoinv_cm=1./29.9792458	            #for converting GHz to inverse cm
	OpticaldepthstodB=434294.5				#converts from cm^-1 to dB/km
	torrperatm=760.
	bartoatm=0.987
	GHztoMHz=1000.
	hc=19.858252418e-24			#planks (J.s) light (cm/s)
	k=1.38e-23					#boltzmann's in J/K or N.m/K
	No=6.02297e23					#Avogadros Number [mole^-1]
	R=8.31432e7					#Rydberg's [erg/mole-K]
	To=300.			                #Ref temp for P/P Catalogue
	dynesperbar=1.e6				#dyne=bar/1e6
	coef=dynesperbar*No/R          #See Appendix D: Using the Poyter-Pickett Catalogs

	PH2=P*H2mr #partial pressures
	PHe=P*Hemr
	PNH3=P*NH3mr

	#calculate vector linewidth
	xi1=0.7653 
	xi2=2./3.
	xi3=1.
	xi12=0.7714 
	xi22=2./3.
	xi32=1.5591 
	Tdiv=To/T
	gnu1=1.6302 
	gnu2=0.75
	gnu3=0.8420 

	gH2=gnu1*PH2
	gHe=gnu2*PHe
	gNH3=gnu3*PNH3*gammaNH3o
	gamma=(gH2)*((Tdiv)**(xi1))+(gHe)*((Tdiv)**(xi2))+gNH3*(300/T)**(xi3)

	delt=-0.0363*gamma
	znu1=1.2747 
	znu2=0.3
	znu3=0.5413
	zH2=znu1*PH2
	zHe=znu2*PHe
	zNH3=znu3*PNH3*gammaNH3o
	zeta=(zH2)*((Tdiv)**(xi12))+(zHe)*((Tdiv)**(xi22))+zNH3*(300/T)**(xi32)

	zetasize=fo.shape[0]
	pst=delt							# answer in GHz
	#Coupling element, pressure shift and dnu or gamma are in GHz, need to convert brlineshape to inverse cm which is done below

	n=f.shape[1]  #returns the number of columns in f
	m=fo.shape[0] #returns the number of rows in fo
	
	# f1 f2 f3 f4 ....fn  n times where n is the number of frequency steps
	# f1 f2 f3 f4 ....fn				in the observation range                            
	# ...
	# f1 f2 f3 f4 ....fn 
	# m times where m is the number of spectral lines

	nones=numpy.ones([1,n])
	mones=numpy.ones([m,1])
	 
	f_matrix=mones*f
	fo_matrix=fo*nones

	# The 10^6 allows use of P(bar) for P(dynes/cm^2)

	eta=3./2.			# for symmetric top molecule
	expo=-(1./T-1./To)*Eo*hc/k
	ST=Io*numpy.exp(expo)	# S(T) =S(To)converted for temperature
	Con=0.9316
	alpha_noshape=Con*coef*(PNH3/To)*((To/T)**(eta+2))*ST  
	#Alpha Max Found

	#Ben Reuven lineshape calculated by the brlineshape function gives the answer in GHz	
	#Here we change from GHz to inverse cm.
	lineshape=0 # 0 for Ben Reuven lineshape, 1 for Gross lineshape
	dnu_matrix=gamma*nones    
	ce_matrix=zeta*nones
	pst_matrix=pst*nones
	if (lineshape<0.5):
		Aa=(2./numpy.pi)*(numpy.power((f_matrix/fo_matrix),2))			

		Bb=(dnu_matrix-ce_matrix)*(numpy.power(f_matrix,2))
		Cc=dnu_matrix+ce_matrix
		Dd=(numpy.power((fo_matrix+pst_matrix),2)) + (numpy.power(dnu_matrix,2))-(numpy.power(ce_matrix,2))
		Ee=numpy.power(f_matrix,2)
		Jj=numpy.power((fo_matrix+pst_matrix),2)
		Gg=numpy.power(dnu_matrix,2)
		Hh=numpy.power(ce_matrix,2)
		Ii=4.*(numpy.power(f_matrix,2))*(numpy.power(dnu_matrix,2))
		Ff=Aa*(Bb+Cc*Dd)/((numpy.power((Ee-Jj-Gg+Hh),2))+Ii)	
	else:
		Aa=(2./numpy.pi)*2*(numpy.power(f_matrix,2))*dnu_matrix			

		Bb=numpy.power((numpy.power(fo_matrix,2)-numpy.power(f_matrix,2)),2)
		Cc=4.*numpy.power(dnu_matrix,2)*numpy.power(f_matrix,2)# watch this one.
		Ff=Aa/(Bb+Cc)  	
	

	Fbr=(1./GHztoinv_cm)*Ff

	alpha_noshape_matrix=alpha_noshape*nones
	br_alpha_matrix=alpha_noshape_matrix*Fbr

	# Uncomment to add effect of rotational lines

	# foRot=[140141.8974 572498.0678 1168451.6877 1214858.6009 1215245.1833 1763525.3868 1763601.8638 1763821.3719 1808935.5500 1810377.7915 2357210.3572 2357726.7497 2358563.2307 2400017.6324 2400578.3942 2402264.8766 2405121.2992 2948410.6478 2948669.3987 2949480.4272 2950814.5336]'/1000
	# IoRot=[2.9612788E-24 1.0606676E-20 2.9273813E-20 8.2904728E-20 3.1664986E-20 2.0686132E-19 9.3579193E-20 6.1713375E-20 9.8487088E-20 6.5024841E-20 1.6891931E-19 1.4258679E-19 1.8191952E-19 3.6718024E-19 1.7525858E-19 1.4807417E-19 1.8922534E-19 4.5287599E-19 2.2129892E-19 2.0430514E-19 3.4039196E-19]'
	# EoRot=[984.1708 0.3967 16.5667 19.4932 15.7763 60.0163 56.3125 45.1906 55.542 44.3993 115.8816 104.7871 86.2612 118.8411 115.1399 104.0254 85.465 198.8972 195.2146 184.1564 165.6913]'
	# expoRot=-(1/T-1/To)*EoRot*hc/k
	# STRot=IoRot.*exp(expoRot)	
	# gammaNH3oRot=[25 0 16.3 0 16.3 0 12.5 20.45 12.5 20.45 11.2 17.3 22.45 0 11.2 17.3 22.45 0 9.8 14.75 19.35]'
	# #gammaNH3oRot=10*ones(21,1)
	# 
	# gNH3rot=gnu3*PNH3*gammaNH3oRot
	# 
	# gammaRot=(gH2)*((Tdiv)^(xi1))+(gHe)*((Tdiv)^(xi2))+gNH3rot*(300/T)^(xi3)
	# 
	# n=size(f,2)  #returns the number of columns in f
	# m=size(foRot,1)
	# nones=ones(1,n)
	# mones=ones(m,1) 
	# f_matrix=mones*f
	# fo_matrix=foRot*nones
	# dnu_matrix=gammaRot*nones
	# Aa=4*(f_matrix.^2).*dnu_matrix			
	# Bb=(fo_matrix.^2-f_matrix.^2).^2	
	# Cc=4*f_matrix.^2.*dnu_matrix.^2
	# FRot=Aa./(Bb+Cc)
	# FbrRot=(1/GHztoinv_cm).*FRot
	# 
	# alpha_rot=coef*(PNH3/To)*((To/T)^(eta+2))*STRot*nones.*FbrRot  
	alpha_opdep=numpy.sum(br_alpha_matrix,axis=0)#+0.55*sum(alpha_rot,1)

	#sums up the element in the matrix to calculate the alpha in optical depths or inverse cm
	#alpha_opdep=sum(br_alpha_matrix,1)

	#C=1
	alphanh3=alpha_opdep*434294.5#*0.907
	#answer in inverse cm converted to dB/km
	return alphanh3
