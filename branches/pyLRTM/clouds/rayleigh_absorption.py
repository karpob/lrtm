def rayleigh_absorption(f,Cloud_Bulk_Density,Material_Density,Complex_Dielectric_Constant):
	import numpy
	#absorption for small spheres in the rayleigh limit. 
	# Cloud Bulk density in grams/cm^3
	# Material Density in grams/cm^3
	#Complex dielectric constant:
	# epp=imaginary component, ep is the real component	
	# 
	# Returns alpha (absorption coefficient in cm^-1 sometimes incorrectly
	# referred to as optical depths).

	c=2.99792458e8 #m/s
	
	epp=numpy.imag(Complex_Dielectric_Constant)
	ep=numpy.real(Complex_Dielectric_Constant)
	f_in_Hz=f*1.e9
	lambdaz_m=c/f_in_Hz #Need to convert this to centimeters to get absorption coefficient in cm^-1
	lambdaz=lambdaz_m*100 # Convert to centimeters
	alpha=((18.*numpy.pi)/lambdaz)*(Cloud_Bulk_Density/Material_Density)*(epp/(((ep+2.)**2)+epp**2))

	return alpha
