def alpha_ortho_h2(T,density):
	from helmholtz_energy_residual import helmholtz_energy_residual
	from ideal_helmholtz_energy import ideal_helmholtz_energy

	N_i=[ -6.83148,0.01,2.11505,4.38353,0.211292,-1.00939,0.142086,-0.87696,0.804927,-0.710775,0.0639688,0.0710858,-0.087654,0.647088]
	t_i=[0.7333,1,1.1372,0.5136,0.5638,1.6248,1.829,2.404,2.105,4.1,7.658,1.259,7.589,3.946]
	d_i=[1,4,1,1,2,2,3,1,3,2,1,3,1,1]
	p_i=[0,0,0,0,0,0,0,1,1,0/1,0/1,0/1,0/1,0/1];
	phi_i=[0,0,0,0,0,0,0,0,0,-1.169,-0.894,-0.04,-2.072,-1.306]
	Beta_i=[0,0,0,0,0,0,0,0,0,-0.4555, -0.4046,-0.0869,-0.4415,-0.5743] 
	gamma_i=[0,0,0,0,0,0,0,0,0,1.5444,0.6627,0.763,0.6587,1.4327] 
	D_i=[0,0,0,0,0,0,0,0,0,0.6366,0.3876,0.9437,0.3976,0.9626]
	Tc=32.22
        rho_c=15.445*2.01594 # mol/l -> kg/m^3 factors of 1000 liter converstion, and gram conversion cancel
        delta=density/rho_c
	tau=Tc/T
	R=8314.472
        R_hydrogen=R/2.01594
	
	ni=2.5   
	ti=0.0
	vi=[2.54151,-2.3661,1.00365,1.22447] 
        ui=[856,1444,2194,6968]
	ho = 6883.3960754547
	so = 141.315444516

	alpha_o=ideal_helmholtz_energy(ni,ti,vi,ui,R,ho,so,rho_c/2.01594,Tc,tau,delta)        
        alpha_r=helmholtz_energy_residual(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i)
        
        absolute=R_hydrogen*T*(alpha_o+alpha_r)
        nomalized=alpha_o+alpha_r
	ideal=alpha_o
	residual=alpha_r
        
	return absolute,nomalized,ideal,residual
