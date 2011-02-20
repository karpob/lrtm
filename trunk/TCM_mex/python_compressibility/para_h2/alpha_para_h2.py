def alpha_para_h2(T,density):
	from helmholtz_energy_residual import helmholtz_energy_residual
	from ideal_helmholtz_energy import ideal_helmholtz_energy

	N_i =[-7.33375,0.01,2.60375,4.66279,0.682390,-1.47078,0.135801,-1.05327,0.328239,-0.0577833,0.0449743,0.0703464,-0.0401766,0.119510]
	t_i=[0.6855,1,1,0.489,0.774,1.133,1.386,1.619,1.162,3.96,5.276,0.99,6.791,3.19]
	d_i=[1,4,1,1,2,2,3,1,3,2,1,3,1,1]
	p_i=[0,0,0,0,0,0,0,1,1,0.0/1.0,0.0/1.0,0.0/1.0,0.0/1.0,0.0/1.0]
	phi_i=[0,0,0,0,0,0,0,0,0,1.7437,-0.5516 ,-0.0634 ,-2.1341 ,-1.777]
	Beta_i=[0,0,0,0,0,0,0,0,0,0.194,-0.2019,-0.0301,-0.2383, -0.3253]
	gamma_i=[0,0,0,0,0,0,0,0,0,0.8048,1.5248,0.6648,0.6832,1.493]
	D_i=[0,0,0,0,0,0,0,0,0,1.5487,0.1785,1.28,0.6319,1.7104]
	Tc=32.938
        rho_c=15.538*2.01594 # mol/l -> kg/m^3 factors of 1000 liter converstion, and gram conversion cancel
        delta=density/rho_c
	tau=Tc/T
	R=8314.472
        R_hydrogen=R/2.01594
	
	ni=2.5   
	ti=0.0
        vi=[4.30256,13.0289,-47.7365,50.0013,-18.6261,0.993973,0.536078]
        ui=[499,826.5,970.8,1166.2,1341.4,5395,10185]

	
	ho = 8172.6404795208
	so = 150.0625663134

	alpha_o=ideal_helmholtz_energy(ni,ti,vi,ui,R,ho,so,rho_c/2.01594,Tc,tau,delta)
	
        
        alpha_r=helmholtz_energy_residual(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i)
        
        absolute=R_hydrogen*T*(alpha_o+alpha_r)
        nomalized=alpha_o+alpha_r
	ideal=alpha_o
	residual=alpha_r
        
	return absolute,nomalized,ideal,residual
