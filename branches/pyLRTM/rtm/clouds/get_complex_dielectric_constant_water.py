def get_complex_dielectric_constant_water(f_in_GHz,T_in_K):
	import numpy
	f=f_in_GHz*1e9
	T=T_in_K-273.15

	ep_wo=88.045-0.4147*T+((6.295e-4)*T*T)+((1.075e-5)*T*T*T) #Static dielectric constant of pure H2O
	ep_winf=4.9 #  high frequency limit of ep_w
	tau_w=((1.1109e-10)-(3.824e-12)*T+(6.938e-14)*T*T-(5.096e-16)*T*T*T)/(2.*numpy.pi) # relaxation time for pure water

	ep_prime=ep_winf+(ep_wo-ep_winf)/(1.+numpy.power((2.*numpy.pi*f*tau_w),2.))

	ep_prime_prime=(2.*numpy.pi*f*tau_w*(ep_wo-ep_winf))/(1.+numpy.power((2.*numpy.pi*f*tau_w),2))

	complex_dielectric_constant=ep_prime+(1j)*ep_prime_prime
	return complex_dielectric_constant
