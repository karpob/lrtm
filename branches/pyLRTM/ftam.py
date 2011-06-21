def ftam(T,tau,tau_a,masterindex)
	import numpy
	# function ftam
	#
	#
	#            Calculate the brightness temperature, and weighting function for a given ray.
	#            
        #    VARIABLE DEFINITIONS:
	#                        
	#                         ->INPUT:	
	#                              ->  T: physical temperature of the layer
	#                              ->tau: total optical depth from bottom to current layer
	#                           -> tau_a: optical depth of each layer across ds
	#                     ->  masterindex: layer index
	#
	#                        <-OUTPUT:
	#                          <-Tatm: Brightiness temperature for each ray
	#                          <-wlayers: weighting function for each ray    

	#JPH function Tatm=ftam(T,tau,tau_a)
	#JPH T is the vector from ?top to bottom of thermal temperatures
	#JPH tau is the vector of integrated absorptions from top to bottom ?
	#JPH tau_a is the vector of absorptions for each layer top to bottom
	#JPH Units are kelvin and opdepths
	#JPH The first layer is the ray from the craft to the surface
	#JPH its tau, tau_a, T, are zero
	#JPH So T(1) should always be zero

	loss=numpy.exp(-tau)
	emit=T[masterindex]*(1.-numpy.exp(-tau_a))
	wlayers=loss*(1.-numpy.exp(-tau_a))
	layers=loss*emit

	Tatm=layers.sum()

	return Tatm,wlayers


