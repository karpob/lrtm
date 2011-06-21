def refractivity_water_vapor_rueger(T,P):
	# This function gives the value for refractivity (N) for a given Temperature and Pressure
	#     P--> is in Bars
	#     T--> is in Kelvins
	#     N <-- N-units (unitless)
	if(T>0.):
		K2=71.97
		K3=375406. # Value given by Rueger 2002
		Pmb=P*1000 # convert to millibars
		N=K2*Pmb/T +K3*Pmb/(T*T)
	else:
     		N=0.
	return N

