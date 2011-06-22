def ftau(kappas,s,masterindex):
        import numpy
	# function ftau
	#
	#
	#             Calculate the optical depth in each layer, and total optical depth up to each layer.
	#
	#             VARIABLE DEFINITIONS:
	#                          
	#                          ->  INPUT:
	#                          
	#                          ->  kappas: The values for absorption coefficient for each layer in cm^-1
	#                          ->       s: The values for ray pathlength for each layer
	#                          -> materindex: indices for each ray/ellipse intersection
	#                        
	#                          <- OUTPUT:
	#
	#                          <- tau_a: optical depth for each layer
	#                          <- tau: total optical depth from bottom to current layer 
	#

	#JPH function tau=ftau(kappa,s)
	#JPH kappa should be total absorption for each layer (a priori)
	#JPH and 's' is the profile of pathlengths in each layer
	#JPH kappa (opdepth/cm), s should be in cm
	#JPH the tau is t(b,c)=sum[(a=b to c) tau(a)
	#JPH where tau(a(0:i))= integral[(0 to i)kappa(s)*ds]
	#JPH So find tau(a)'s which is the tau for each layer
	#JPH then find the tau(b,c)'s which  is the integrated from b to c
	#JPH Must be column vector

	# tau(a)'s
	# 

	tau_a=kappas[masterindex]*s		# now units of opdepth

	# tau(b,c)'s
	# I think this works cumsum?  For column vector use cumsum(X,1)

	tau1=numpy.cumsum(tau_a,axis=0)

	st=tau1.shape[0]
	tau=numpy.zeros([st])
	tau[1:len(tau)]=tau1[0:st-1]
	return tau_a, tau
