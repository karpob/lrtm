def alpha_normal_h2(T,density):
	from helmholtz_energy_residual import helmholtz_energy_residual
	from ideal_helmholtz_energy import ideal_helmholtz_energy
	from pylab import zeros
	N_i=[-6.93643,0.01,2.1101,4.52059,0.732564,-1.34086,0.130985,-0.777414,0.351944,-0.0211716,0.0226312,0.032187,-0.0231752,0.0557346];
	t_i=[0.6844,1,0.989,0.489,0.803,1.1444,1.409,1.754,1.311,4.187,5.646,0.791,7.249,2.986]
	d_i=[1,4,1,1,2,2,3,1,3,2,1,3,1,1]
	p_i=[0,0,0,0,0,0,0,1,1,0/1,0/1,0/1,0/1,0/1]
	phi_i=[0,0,0,0,0,0,0,0,0,1.685,0.489,0.103,2.506,1.607]
	Beta_i=[0,0,0,0,0,0,0,0,0, 0.171,0.2245,0.1304,0.2785,0.3967] 
	gamma_i=[0,0,0,0,0,0,0,0,0,0.7164,1.3444,1.4517,0.7204,1.5445] 
	D_i=[0,0,0,0,0,0,0,0,0,1.506,0.156,1.736,0.67,1.662]
	betas=zeros(len(gamma_i))
	Tc=33.145
        rho_c=15.508*2.01594 # mol/l -> kg/m^3 factors of 1000 liter converstion, and gram conversion cancel
        delta=density/rho_c
	tau=Tc/T
	R=8314.472
        R_hydrogen=R/2.01594
	
	ni=2.5   
	ti=0.0
        vi=[1.616,-0.4117,-0.792,0.758,1.217] 
        ui=[531,751,1989,2484,6859]
        
	R=8.314472
	ho = 7206.9069892047
	so = 143.4846187346

	alpha_o=ideal_helmholtz_energy(ni,ti,vi,ui,R,ho,so,rho_c/2.01594,Tc,tau,delta)
	
        
        alpha_r=helmholtz_energy_residual(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i)
        
        absolute=R_hydrogen*T*(alpha_o+alpha_r)
        nomalized=alpha_o+alpha_r
	ideal=alpha_o
	residual=alpha_r
        
	return absolute,nomalized,ideal,residual
