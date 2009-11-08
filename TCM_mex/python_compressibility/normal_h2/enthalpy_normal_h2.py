def enthalpy_normal_h2(T,density):
	from helmholtz_functions.ideal_helmholtz_energy_from_Cp_o_over_R import ideal_helmholtz_energy_from_Cp_o_over_R
	from helmholtz_functions.ideal_helmholtz_energy_dtau_from_Cp_o_over_R import ideal_helmholtz_energy_dtau_from_Cp_o_over_R
	from helmholtz_functions.ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R import ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R
	from helmholtz_functions.d_alpha_d_tau import d_alpha_d_tau 
	from helmholtz_functions.d_alpha_d_delta import d_alpha_d_delta
        from helmholtz_functions.helmholtz_energy_residual import helmholtz_energy_residual
        from helmholtz_functions.d_alpha_d_tau_d_tau import d_alpha_d_tau_d_tau
        from helmholtz_functions.d_alpha_d_delta_d_delta import d_alpha_d_delta_d_delta
        from helmholtz_functions.d_alpha_d_delta_d_tau import d_alpha_d_delta_d_tau
        from helmholtz_functions.ideal_helmholtz_energy_from_Cp_o_over_R_vector import ideal_helmholtz_energy_from_Cp_o_over_R_vector
        from helmholtz_functions.ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R_vector import ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R_vector
        from thermodynamic_functions.Virial_B_vector import Virial_B_vector
        from math import sqrt
        from pylab import zeros
        from numpy import asarray
	N_i=[-6.93643,0.01,2.1101,4.52059,0.732564,-1.34086,0.130985,-0.777414,0.351944,-0.0211716,0.0226312,0.032187,-0.0231752,0.0557346];
	t_i=[0.6844,1,0.989,0.489,0.803,1.1444,1.409,1.754,1.311,4.187,5.646,0.791,7.249,2.986]
	d_i=[1,4,1,1,2,2,3,1,3,2,1,3,1,1]
	p_i=[0,0,0,0,0,0,0,1,1,0/1,0/1,0/1,0/1,0/1]
	phi_i=[0,0,0,0,0,0,0,0,0,1.685,0.489,0.103,2.506,1.607]# Not, we need to use -sign here
	Beta_i=[0,0,0,0,0,0,0,0,0, 0.171,0.2245,0.1304,0.2785,0.3967] #Note, we need to use - sign here
	gamma_i=[0,0,0,0,0,0,0,0,0,0.7164,1.3444,1.4517,0.7204,1.5445] 
	D_i=[0,0,0,0,0,0,0,0,0,1.506,0.156,1.736,0.67,1.6620]
	Tc=33.145
        rho_c=15.508*2.01594 # mol/l -> kg/m^3 factors of 1000 liter converstion, and gram conversion cancel
        delta=density/rho_c
	tau=Tc/T
	R=8.314472
        R_hydrogen=R/2.01594
	
	ni=2.5   
	ti=0.0
        vi=[1.616,-0.4117,-0.792,0.758,1.217] 
        ui=[531,751,1989,2484,6859]
        
	R=8.314472           #J/mol/K
	ho = 7206.9069892047 #J/mol
	so = 143.4846187346  #J/mol/K
	n_power_terms=9
	n_gaussian_terms=5
	n_critical_terms=0
	
	n_ideal_gas_terms_pow=0
	n_ideal_gas_terms_exp=5
	

	RES_a     = zeros(len(ui))
        RES_b     = zeros(len(ui))
        RES_B     = zeros(len(ui))
        RES_C     = zeros(len(ui))
        RES_D     = zeros(len(ui))
        RES_A     = zeros(len(ui))

	To=273.15        
        rho_o=((0.001/(R*To))*1000)
        h2_component_params={'N_i':asarray([-6.93643,0.01,2.1101,4.52059,0.732564,-1.34086,0.130985,-0.777414,0.351944,-0.0211716,0.0226312,0.032187,-0.0231752,0.0557346]),
                             't_i':asarray([0.6844,1,0.989,0.489,0.803,1.1444,1.409,1.754,1.311,4.187,5.646,0.791,7.249,2.986]),
                             'd_i':asarray([1,4,1,1,2,2,3,1,3,2,1,3,1,1]),
                             'p_i':asarray([0,0,0,0,0,0,0,1,1,0/1,0/1,0/1,0/1,0/1]),
                             'phi_i':asarray([0,0,0,0,0,0,0,0,0,1.685,0.489,0.103,2.506,1.607]),
                             'Beta_i':asarray([0,0,0,0,0,0,0,0,0, 0.171,0.2245,0.1304,0.2785,0.3967]), 
                             'gamma_i':asarray([0,0,0,0,0,0,0,0,0,0.7164,1.3444,1.4517,0.7204,1.5445]), 
                             'D_i':asarray([0,0,0,0,0,0,0,0,0,1.506,0.156,1.736,0.67,1.6620]),
                             'Tc':float(33.145),
                             'rho_c':float(15.508*2.01594), # mol/l -> kg/m^3 factors of 1000 liter converstion, and gram conversion cancel
                             'R':float(8.314472),
                             'ni':float(2.5),   
                             'ti':float(0.0),
                             'vi':asarray([1.616,-0.4117,-0.792,0.758,1.217]), 
                             'ui':asarray([531,751,1989,2484,6859]),
                             #'R':float(8.314472)           #J/mol/K
                             'ho':float(7206.9069892047), #J/mol
                             'so':float(143.4846187346),  #J/mol/K
                             'n_power_terms':int(9),
                             'n_power_terms_wo_exp':int(7),
                             'n_power_terms_w_exp':int(2),
                             'n_gaussian_terms':int(5),
                             'n_critical_terms':int(0),
                             'n_ideal_gas_terms_pow':int(0),
                             'n_ideal_gas_terms_exp':int(5),
                             'RES_a':zeros(5),
                             'RES_b':zeros(5),
                             'RES_B':zeros(5),
                             'RES_C':zeros(5),
                             'RES_D':zeros(5),
                             'RES_A':zeros(5),
                             'To':float(273.15),        
                             'rho_o':float((0.001/(8.314472*273.15))*1000),
                             'M_amu':float(2.01594),
                             'ideal_eqn_type':'Cp'}
        tau_v=asarray([tau,tau])
        delta_v=asarray([delta,delta])
        
        Virial_b=Virial_B_vector(tau_v,h2_component_params,h2_component_params,h2_component_params,1.0,1.0,rho_c/2.01594,Tc,'pure')
        
	#Calculate Ideal Terms
	ideal=ideal_helmholtz_energy_from_Cp_o_over_R(ni,ti,vi,ui,R,ho,so,rho_c/2.01594,Tc,tau,delta,n_ideal_gas_terms_pow,n_ideal_gas_terms_exp,To,rho_o)
        ideal_v=ideal_helmholtz_energy_from_Cp_o_over_R_vector(ni,ti,vi,ui,R,ho,so,rho_c/2.01594,Tc,tau_v,delta_v,n_ideal_gas_terms_pow,n_ideal_gas_terms_exp,To,rho_o)
        #print ideal-ideal_v
        dalpha_o=ideal_helmholtz_energy_dtau_from_Cp_o_over_R(ni,ti,vi,ui,R,ho,so,rho_c/2.01594,Tc,tau,delta,n_ideal_gas_terms_pow,n_ideal_gas_terms_exp,To,rho_o)
	
	dalpha_o_tau_tau=ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R(ni,ti,vi,ui,R,ho,so,rho_c/2.01594,Tc,tau,delta,n_ideal_gas_terms_pow,n_ideal_gas_terms_exp,To,rho_o)
	dalpha_o_tau_tau_vector=ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R_vector(ni,ti,vi,ui,R,ho,so,rho_c/2.01594,Tc,tau_v,delta_v,n_ideal_gas_terms_pow,n_ideal_gas_terms_exp,To,rho_o)
	#print dalpha_o_tau_tau-dalpha_o_tau_tau_vector
	#Calculate Residual Terms
	residual=helmholtz_energy_residual(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	
	dalpha_tau=d_alpha_d_tau(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	
	dalpha_delta=d_alpha_d_delta(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	print T, Virial_b,dalpha_delta
	dalpha_tau_tau=d_alpha_d_tau_d_tau(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
        
        dalpha_delta_delta=d_alpha_d_delta_d_delta(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
        
        dalpha_delta_tau=d_alpha_d_delta_d_tau(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	
	
        #compressibility Density<-->Pressure
	Z=delta*dalpha_delta+float(1.0)
        #Enthalpy
	h1=R_hydrogen*T*1.0
	h2=R_hydrogen*T*tau*(dalpha_o+dalpha_tau)
	h3=R_hydrogen*T*delta*dalpha_delta
	h=h1+h2+h3
	#Entropy
	s=R_hydrogen*(tau*dalpha_o+tau*dalpha_tau-ideal-residual)
	#Absolute Helmholtz Energy
	absolute_helmholtz_energy=R_hydrogen*T*(ideal+residual)
	#Isochoric heat capacity
	cv=-1.0*R_hydrogen*(tau*tau*(dalpha_o_tau_tau+dalpha_tau_tau))
	#Isobaric heat capacity
	cp=cv+R_hydrogen*pow((1.0+delta*dalpha_delta-delta*tau*dalpha_delta_tau),2)/(1.0+2.0*delta*dalpha_delta+delta*delta*dalpha_delta_delta)
	
	# Speed of Sound?
	
	a1 = 1 + delta*dalpha_delta - delta*tau*dalpha_delta_tau;
        b1 = tau*tau*(dalpha_o_tau_tau + dalpha_tau_tau);
        w = 1 + 2*delta*dalpha_delta + delta*delta*dalpha_delta_delta - a1*a1/b1;
        speed_o_sound=sqrt(R_hydrogen*T*w*1000)
        
        return Z,h,s,absolute_helmholtz_energy,cv,cp,speed_o_sound,R_hydrogen*T*ideal,residual*R_hydrogen*T
