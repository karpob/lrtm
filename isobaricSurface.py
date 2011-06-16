def isobaricSurface(theta,a,Rref,Rp,J,tau):
	import numpy
	GM=126712765.
	cosSq=numpy.power(numpy.cos(theta),2)
	cosQuad=numpy.power(numpy.cos(theta),4)
	oneSixteenth=1./16.
	fiveSixteenths=5./16.
	
	tau=tauRotation(J,a,Rp)
	
	omega=2.*numpy.pi/(tau)
	q=(omega*omega*a*a*a)/GM
	J=[J[0]*numpy.power(Rref/a,2)*numpy.power(q,-1),J[1]*numpy.power(Rref/a,4)*numpy.power(q,-2),J[2]*numpy.power(Rref/a,6)*numpy.power(q,-3)]
	r1 =-0.5*(1. + 3.*J[0])*cosSq
	r2=-0.25*(2. + 3.*J[0]-9.*J[0]*J[0]-15.*J[1])*cosSq+0.125*(6.+6.*(1-6.*J[0])*J[0]-35.*J[1])*cosQuad
	
	
	r3=-1.0*oneSixteenth*(8.-27.*J[0]*(2.*(1.-J[0])*J[0] -5.*J[1]) -45.*J[1] + 105.*J[2])*cosSq \
	   +fiveSixteenths*(6.-8.*J[1] + 3.*J[0]*(2. + 6.*J[0]*(3.*J[0] -1.) + 43.*J[1]) + 63.*J[2])*cosQuad \
	   -1./16.*(24. + 35.*J[1] + 18.*J[0]*(2. + J[0] + 21.*J[0]*J[0] + 35.*J[1]) + 231.*J[2])*cosQuad*cosSq
	
	r=a*(1.+r1*q+r2*q*q+r3*q*q*q)

	return r
