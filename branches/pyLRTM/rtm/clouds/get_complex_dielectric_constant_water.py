def get_complex_dielectric_constant_water(f_in_GHz,T_in_K):
	"""
	Function based on Ulaby, F. T., R. K. Moore, and A. K. Fung (1986), Microwave Remote Sensing Active
    and Passive: From Theory to Applications, vol. III, Artech House, Inc.
	"""
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
def get_complex_dielectric_constant_aq(f_in_GHz,T_in_K,C):
	"""
	Function based on Duong Master's Thesis.
	Frequency units in GHz,
	T units in K, converted in routine to C.
	C is the concentration (volume fraction).
	"""
	#convert to T_in_C
	T_in_C=T_in_K-273.15
	import numpy as np
	nu=f_in_GHz
	x=[]
	x.append(float(5.723))
	x.append(float(2.2379e-2))
	x.append(float(-7.1237e-4))
	x.append(float(5.0478e0))
	x.append(float(-7.0315e-2))
	x.append(float(6.0059e-4))
	x.append(float(3.6143e0))
	x.append(float(2.8841e-2))
	x.append(float(1.3652e-1))
	x.append(float(1.4825e-3))
	x.append(float(2.4166e-4))
	#here comes the stuff that duong did.
	x.append(float(-7.8e1))
	x.append(float(1.090e-2))
	x.append(float(5.860e-2))
	x.append(float(2.264e2))
	x.append(float(2.310e-2))
	x.append(float(1.290e1))
	x.append(float(2.477e1))
	
	ep_s=(3.70886e4 - 8.2168e1*T_in_C)/(4.21854e2+T_in_C)
	ep_1=x[0]+x[1]*T_in_C +x[2]*T_in_C*T_in_C
	nu_1=(45.+T_in_C)/(x[3]+x[4]*T_in_C+x[5]*T_in_C*T_in_C)
	ep_inf=x[6]+x[7]*T_in_C
	nu_2=(45.+T_in_C)/(x[8]+x[9]*T_in_C+x[10]*T_in_C*T_in_C)
	try:
		deltaNH3=(x[11]*C*nu**x[12])/(T_in_C**x[13])-(1j)*((x[14]*C*nu**x[15])/(T_in_C**x[16])+x[17]*C)
	except:
		deltaNH3=0.0
	
	term1=(ep_s-ep_1)/(1.+(1j)*(nu/nu_1))
	term2=(ep_1-ep_inf)/(1.+(1j)*(nu/nu_2))
	term3=ep_inf
	term4=0.#it;s really conductivity (asssuming it's zero) over 2*pi*ep_o*nu
	term5=deltaNH3
	
	ep=term1+term2+term3-term4+term5
	return ep
