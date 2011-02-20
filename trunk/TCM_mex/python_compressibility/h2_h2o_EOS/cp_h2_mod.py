def cp_h2_mod(T,Pin):
	from helmholtz_functions.ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R import ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R
	from helmholtz_functions.d_alpha_d_tau import d_alpha_d_tau 
	from helmholtz_functions.d_alpha_d_delta import d_alpha_d_delta
        from helmholtz_functions.d_alpha_d_tau_d_tau import d_alpha_d_tau_d_tau
        from helmholtz_functions.d_alpha_d_delta_d_delta import d_alpha_d_delta_d_delta
        from helmholtz_functions.d_alpha_d_delta_d_tau import d_alpha_d_delta_d_tau
        from numpy import asarray,zeros
        from scipy.optimize import fsolve
        from residual_pressure import residual_pressure
        bars_to_kPa=100
        
        P=Pin*bars_to_kPa
        R=8.314472
        R_hydrogen=R/2.01594
        M_amu=2.01594
        density_guess=P/(R_hydrogen*T)
        rho_c=15.508*2.01594 # mol/l -> kg/m^3 factors of 1000 liter converstion, and gram conversion cancel
        delta_guess=density_guess/rho_c
        
        
        ########################################
        #Hydrogen parameters from Leachman,2007
        ########################################
	Tc=33.145
        
	tau=Tc/T
	R=8.314472
        R_hydrogen=R/2.01594
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
        #delta=density/rho_c
	tau=Tc/T
	
	
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
	n_power_terms_w_exp=0
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
	#############################################
	## End parameters -> should be data structure
	#############################################
	
	delta=fsolve(residual_pressure,delta_guess,args=(P,tau,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A,R,Tc,rho_c,M_amu))
        
	
	dalpha_o_tau_tau=ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R(ni,ti,vi,ui,R,ho,so,rho_c/2.01594,Tc,tau,delta,n_ideal_gas_terms_pow,n_ideal_gas_terms_exp,To,rho_o)
	dalpha_tau=d_alpha_d_tau(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	
	dalpha_delta=d_alpha_d_delta(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	dalpha_tau_tau=d_alpha_d_tau_d_tau(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
        
        dalpha_delta_delta=d_alpha_d_delta_d_delta(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
        
        dalpha_delta_tau=d_alpha_d_delta_d_tau(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	
	#Isochoric heat capacity
	cv=-1.0*R_hydrogen*(tau*tau*(dalpha_o_tau_tau+dalpha_tau_tau))
	#Isobaric heat capacity
	cp=cv+R_hydrogen*pow((1.0+delta*dalpha_delta-delta*tau*dalpha_delta_tau),2)/(1.0+2.0*delta*dalpha_delta+delta*delta*dalpha_delta_delta)
	#default is in kJ/kg-K convert to erg/K/mol for DeBoer's TCM
	#  kJ    1000J    1e7 ergs    1 kg    M_amu g      erg
	#----- * ----- * --------- * ------ * --------  = ------- 
	# kg*K    kJ      1J         1000 g    mol         K mol
	        
        return cp*M_amu*1e7
        
