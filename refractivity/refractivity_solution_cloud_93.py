def refractivity_solution_cloud_93(T,D_sol,rho,f):
        import numpy
	# Calculates the real (Np) and imaginary (Npp) parts of refractvitiy
	#
	#     Input Variables:
	#                  -->T: Temperature in Kelvins
	#                  -->D_Sol: Cloud Density in g/cm^3
	#                  -->rho: Density of Cloud Material in g/cm^3
	#                  -->f: frequency in GHz
	#
	#     Output Variables:
	#                  <--Np: Nprime or real part of refractivity
	#                  <--Npp: N double prime or imaginary part of refractivity

	D_sol_g_m3=(D_sol*1e6)/rho # convert to g/m^3 This is necessary, since N', N'' are small parts of total refractive index
	# Use Refractivity from Liebe, 1989 (most closely matches his notation in
	# paper), but 1993 will work too

	V=300./T
	fD=20.20-146.4*(V-1)+316.*numpy.power((V-1),2)
	fS=39.8*fD
	Eps=103.3*(V-1)+77.66
	Epinf=0.0671*Eps
	Eopt=3.52
	f_fD=f+(1j)*fD
	f_fS=f+(1j)*fS
	ZEp=Eps-f*((Eps-Epinf)/f_fD+(Epinf-Eopt)/f_fS)
	ZNw=1.5*(D_sol_g_m3)*(((ZEp-1.)/(ZEp+2.))-1.+3./(Eps+2.))
	Np=numpy.real(ZNw)
	Npp=numpy.imag(ZNw)
	print Np, Npp
	return Np,Npp
