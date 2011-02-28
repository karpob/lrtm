def OrtonCIA(f,T,P_H2,P_He,P_CH4):
	"""
	This function uses data supplied by Glenn Orton from his quantum scattering
	code for H2-H2, H2-He, and H2-CH4. The data is interpolated using interp2
 	with spline option (along frequency, and temperature grid).
	Input
	     f     --> Frequency in GHz
	     T     --> Temperature in K
	     P_H2  --> Patial Pressure from H2 (bars)
	     P_He  --> Partial Pressure from He (bars)
	     P_CH4 --> Partial Pressure from CH4 (bars)
	     normal_or_equilibrium --> select normal (1) or eq(0) hydrogen

	     Output
	     alpha <-- total CIA absorption from H2
	     alpha_H2_prime  <-- CIA absorption from H2-H2
	     alpha_He_prime  <-- CIA absorption from H2-He
	     alpha_CH4_prime <-- CIA absorption from H2-CH4
	"""
	from numpy import exp,linspace,log
	import scipy.interpolate 
	from ReadasciiCIA import ReadasciiCIA
	import scipy,numpy,pylab,math
	################################################################
	# Convert stuff to wavenumbers and amagats...let the fun begin.#
	################################################################
	GHztoHz=1.0e9
	f=f*GHztoHz
	c=2.99792458e10 #cm/sec
	nu=f/c

	Lo=2.68719e19 #Loschmidt number molecules/cm^3 at stp

	Pascal_per_bar=100000.
	cm3_per_m3=1.e6
	grams_per_kg=1000.
	molecules_per_mole=6.0221415e23

	RH2=4124.18 #Nm/Kg/K
	RHe=2077.
	RCH4=518.3

	rho_he=(Pascal_per_bar*P_He)/(RHe*T)
	rho_he=rho_he/cm3_per_m3
	rho_he=rho_he*grams_per_kg

	rho_h2=(Pascal_per_bar*P_H2)/(RH2*T)
	rho_h2=rho_h2/cm3_per_m3
	rho_h2=rho_h2*grams_per_kg

	rho_ch4=(Pascal_per_bar*P_CH4)/(RCH4*T)
	rho_ch4=rho_ch4/cm3_per_m3
	rho_ch4=rho_ch4*grams_per_kg

	grams_per_mole_h2=2.*1.007822
	grams_per_mole_He=4.003
	grams_per_mole_CH4=16.04

	rho_h2=rho_h2/grams_per_mole_h2
	rho_h2=rho_h2*molecules_per_mole
	amagat_h2=rho_h2/Lo
	amagat_sq_h2=amagat_h2**2

	rho_he=rho_he/grams_per_mole_He
	rho_he=rho_he*molecules_per_mole
	amagat_he=rho_he/Lo
	amagat_sq_he=amagat_he**2

	rho_ch4=rho_ch4/grams_per_mole_CH4
	rho_ch4=rho_ch4*molecules_per_mole
	amagat_ch4=rho_ch4/Lo
	amagat_sq_ch4=amagat_ch4**2
	######################################################################    

	########################################################################
	# Enough having fun with amagats, time to read in the data             #
	# Read in data supplied by Glenn Orton from His quantum scattering code#    
	########################################################################
	
	nTemps,minTemp,maxTemp,wave_numbers,alpha_eH2,alpha_eH2_He,alpha_eH2_CH4=ReadasciiCIA('h2Parameters/cia.longrev08.data')

	#####################################################################
	#  Start interpolation scheme                                       #
	#####################################################################
	# "evenly spaced T in log space used for calculation 
	#   (ie. 10 temperatures in Orton's file)

	#Ti=[40,51.667,66.724,86.177,111.302,143.753,185.664,239.794,309.705,400]
	logTi=linspace(log(minTemp),log(maxTemp),num=nTemps)
	
	# interpolate values for given wavenumber (nu), and while we're at it multiply
	# by amagat^2 of each constituent
	f=scipy.interpolate.interp2d(logTi,wave_numbers,alpha_eH2,kind='cubic')
	
	#f=scipy.interpolate.RectBivariateSpline(logTi,wave_numbers,alpha_eH2)
    	alpha_H2_prime_1=amagat_sq_h2*exp(f(log(T),nu))
	
	#f=scipy.interpolate.RectBivariateSpline(logTi,wave_numbers,alpha_eH2_He)
    	#alpha_He_prime_1=amagat_he*amagat_h2*exp(f(log(T),nu))

	#f=scipy.interpolate.RectBivariateSpline(logTi,wave_numbers,alpha_eH2_CH4)
	#alpha_CH4_prime_1=amagat_h2*amagat_ch4*exp(f(log(T),nu))

	# Add up alpha's of H2-H2, H2-He, H2-CH4
	#alpha_1=alpha_H2_prime_1+alpha_He_prime_1+alpha_CH4_prime_1
	
    	#alpha=alpha_1
	
	return alpha_H2_prime_1#alpha,alpha_H2_prime_1,alpha_He_prime_1,alpha_CH4_prime_1
