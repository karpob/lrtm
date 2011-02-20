def enthalpy_methane(T,density):
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
        from helmholtz_functions.ideal_helmholtz_energy_dtau_from_Cp_o_over_R_vector import ideal_helmholtz_energy_dtau_from_Cp_o_over_R_vector
        from helmholtz_functions.ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R_vector import ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R_vector
        from thermodynamic_functions.Virial_B_vector import Virial_B_vector
        from math import sqrt
        from pylab import zeros
        from numpy import asarray
        
        N_i=[0.43679010280e-01,
             0.67092361990, 
             -0.17655778590e01, 
             0.85823302410, 
             -0.12065130520e01,
             0.51204672200, 
             -0.40000107910e-03 ,
             -0.12478424230e-01,
             0.31002697010e-01,
             0.17547485220e-02,
             -0.31719216050e-05,
             -0.22403468400e-05,
             0.29470561560e-06,
             0.18304879090, 
             0.15118836790,
             -0.42893638770,
             0.68940024460e-01,
             -0.14083139960e-01,
             -0.30630548300e-01,
             -0.29699067080e-01,
             -0.19320408310e-01,
             -0.11057399590,
             0.99525489950e-01,
             0.85484378250e-02, 
             -0.61505556620e-01,
             -0.42917924230e-01,
             -0.18132072900e-01,
             0.34459047600e-01,
             -0.23859194500e-02,
             -0.11590949390e-01,
             0.66416936020e-01,
             -0.23715495900e-01, 
             -0.39616249050e-01, 
             -0.13872920440e-01,
             0.33894895990e-01,
             -0.29273787530e-02,
             0.93247999460e-04,
             -0.62871715180e01,
             0.12710694670e02,
             -0.64239534660e01]  
	t_i= [-0.5, 
	       0.5,
	       1.0,
	       0.5,
	       1.0,
	       1.5,
	       4.5,
	       0.0,
	       1.0,
	       3.0,
	       1.0,
	       3.0,
	       3.0,
	       0.0,
	       1.0,
	       2.0,
	       0.0,
	       0.0,
	       2.0,
	       2.0,
	       5.0,
	       5.0,
	       5.0,
	       2.0,
	       4.0,
	       12.0,
	       8.0,
	       10.0,
	       10.0,
	       10.0,
	       14.0,
	       12.0,
	       18.0,
	       22.0,
	       18.0,
	       14.0,
	       2.0,
	       0.0,
	       1.0,
	       2.0]
	d_i=[1.0,
	     1.0,
	     1.0,
	     2.0,
	     2.0,
	     2.0,
	     2.0,
	     3.0,
	     4.0,
	     4.0,
	     8.0,
	     9.0,
	     10.0,
	     1.0,
	     1.0,
	     1.0,
	     2.0,
	     4.0,
	     5.0,
	     6.0,
	     1.0,
	     2.0,
	     3.0,
	     4.0,
	     4.0,
	     3.0,
	     5.0,
	     5.0,
	     8.0,
	     2.0,
	     3.0,
	     4.0,
	     4.0,
	     4.0,
	     5.0,
	     6.0,
	     2.0,
	     0.0,
	     0.0,
	     0.0]    
	p_i=[0,
	     0,
	     0,
	     0,
	     0,
	     0,
	     0,
	     0,
	     0,
	     0,
	     0,
	     0,
	     0,
	     1,
	     1,
	     1,
	     1,
	     1,
	     1,
	     1,
	     2,
	     2,
	     2,
	     2,
	     2,
	     3,
	     3,
	     3,
	     3,
	     4,
	     4,
	     4,
	     4,
	     4,
	     4,
	     4,
	     0/1,
	     0/1,
	     0/1,
	     0/1]  
		
	phi_i=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20.0,40.0,40.0,40.0]# Not, we need to use -sign here
	Beta_i=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,200,250,250,250] #Note, we need to use - sign here
	gamma_i=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.07,1.11,1.11,1.11] 
	D_i=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.0,1.0,1.0,1.0]
	Tc=190.564
	
        rho_c=10.139128*16.0428 # mol/l -> kg/m^3 factors of 1000 liter converstion, and gram conversion cancel
        delta=density/rho_c
	tau=Tc/T
	R=8.31451
        R_methane=R/16.0428
	
	ni=4.0016  
	ti=0  
        vi=[0.84490000e-02,4.6942000,3.4865000,1.6572000,1.4115000] 
        ui=[648.0,1957.0,3895.0,5705.0,15080.0]

	R=8.31451          #J/mol/K
	ho = 8295.6883966242294	 #J/mol -       -
	so = 28.384819963016852  #J/mol/K
	rho_o=26.326811491312679 #mol/L
	To=111.66720547358069
	n_power_terms=36
	n_gaussian_terms=4
	n_critical_terms=0
	
	n_ideal_gas_terms_pow=0
	n_ideal_gas_terms_exp=5
	

	RES_a     = zeros(len(ui))
        RES_b     = zeros(len(ui))
        RES_B     = zeros(len(ui))
        RES_C     = zeros(len(ui))
        RES_D     = zeros(len(ui))
        RES_A     = zeros(len(ui))

	      
        
        ch4_component_params={'N_i':asarray([0.43679010280e-01,0.67092361990, -0.17655778590e01, 0.85823302410, -0.12065130520e01,0.51204672200, -0.40000107910e-03 ,-0.12478424230e-01,0.31002697010e-01,0.17547485220e-02,-0.31719216050e-05,-0.22403468400e-05,0.29470561560e-06,0.18304879090, 0.15118836790,-0.42893638770,0.68940024460e-01,-0.14083139960e-01,-0.30630548300e-01,-0.29699067080e-01,-0.19320408310e-01,-0.11057399590,0.99525489950e-01,0.85484378250e-02, -0.61505556620e-01,-0.42917924230e-01,-0.18132072900e-01,0.34459047600e-01,-0.23859194500e-02,-0.11590949390e-01,0.66416936020e-01,-0.23715495900e-01, -0.39616249050e-01, -0.13872920440e-01,0.33894895990e-01,-0.29273787530e-02,0.93247999460e-04,-0.62871715180e+01,0.12710694670e+02,-0.64239534660e+01]),
                             't_i':asarray([-0.5, 0.5,1.0,0.5,1.0,1.5,4.5,0.0,1.0,3.0,1.0,3.0,3.0,0.0,1.0,2.0,0.0,0.0,2.0,2.0,5.0,5.0,5.0,2.0,4.0,12.0,8.0,10.0,10.0,10.0,14.0,12.0,18.0,22.0,18.0,14.0,2.0,0.0,1.0,2.0]),
                             'd_i':asarray([1.0,1.0,1.0,2.0,2.0,2.0,2.0,3.0,4.0,4.0,8.0,9.0,10.0,1.0,1.0,1.0,2.0,4.0,5.0,6.0,1.0,2.0,3.0,4.0,4.0,3.0,5.0,5.0,8.0,2.0,3.0,4.0,4.0,4.0,5.0,6.0,2.0,0.0,0.0,0.0]),
                             'p_i':asarray([0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,4,4,4,4,4,2,2,2,2]),
                             'phi_i':asarray([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,20.0,40.0,40.0,40.0]),
                             'Beta_i':asarray([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,200,250,250,250]), 
                             'gamma_i':asarray([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.07,1.11,1.11,1.11]), 
                             'D_i':asarray([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1]),
                             'Tc':float(190.564),
                             'rho_c':float(10.139128*16.0428), # mol/l -> kg/m^3 factors of 1000 liter converstion, and gram conversion cancel
                             'R':float(8.31451),
                             'ni':float(4.0016),   
                             'ti':float(0.0),
                             'vi':asarray([0.84490000e-02,4.6942000,3.4865000,1.6572000,1.4115000]), 
                             'ui':asarray([648.0,1957.0,3895.0,5705.0,15080.0]),
                             #'R':float(8.314472)           #J/mol/K
                             'ho':float(8295.6883966242294), #J/mol
                             'so':float(28.384819963016852),  #J/mol/K
                             'n_power_terms':int(36),
                             'n_power_terms_wo_exp':int(13),
                             'n_power_terms_w_exp':int(36-13),
                             'n_gaussian_terms':int(4),
                             'n_critical_terms':int(0),
                             'n_ideal_gas_terms_pow':int(0),
                             'n_ideal_gas_terms_exp':int(5),
                             'RES_a':zeros(4),
                             'RES_b':zeros(4),
                             'RES_B':zeros(4),
                             'RES_C':zeros(4),
                             'RES_D':zeros(4),
                             'RES_A':zeros(4),
                             'To':float(111.66720547358069),        
                             'rho_o':float(26.326811491312679),
                             'M_amu':float(16.0428),
                             'ideal_eqn_type':'Cp'}
        tau_v=asarray([tau,tau])
        delta_v=asarray([delta,delta])
        
        Virial_b=Virial_B_vector(tau_v,ch4_component_params,ch4_component_params,ch4_component_params,1.0,1.0,rho_c/16.0428,Tc,'pure')
        
	#Calculate Ideal Terms
	ideal=ideal_helmholtz_energy_from_Cp_o_over_R(ni,ti,vi,ui,R,ho,so,rho_c/16.0428,Tc,tau,delta,n_ideal_gas_terms_pow,n_ideal_gas_terms_exp,To,rho_o)
        ideal_v=ideal_helmholtz_energy_from_Cp_o_over_R_vector(ni,ti,vi,ui,R,ho,so,rho_c/16.0428,Tc,tau_v,delta_v,n_ideal_gas_terms_pow,n_ideal_gas_terms_exp,To,rho_o)
        #print ideal-ideal_v
        dalpha_o=ideal_helmholtz_energy_dtau_from_Cp_o_over_R(ni,ti,vi,ui,R,ho,so,rho_c/16.0428,Tc,tau,delta,n_ideal_gas_terms_pow,n_ideal_gas_terms_exp,To,rho_o)
	dalpha_o_v=ideal_helmholtz_energy_dtau_from_Cp_o_over_R_vector(ni,ti,vi,ui,R,ho,so,rho_c/16.0428,Tc,tau_v,delta_v,n_ideal_gas_terms_pow,n_ideal_gas_terms_exp,To,rho_o)
	print dalpha_o-dalpha_o_v
	dalpha_o_tau_tau=ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R(ni,ti,vi,ui,R,ho,so,rho_c/16.0428,Tc,tau,delta,n_ideal_gas_terms_pow,n_ideal_gas_terms_exp,To,rho_o)
	dalpha_o_tau_tau_vector=ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R_vector(ni,ti,vi,ui,R,ho,so,rho_c/16.0428,Tc,tau_v,delta_v,n_ideal_gas_terms_pow,n_ideal_gas_terms_exp,To,rho_o)
	#print dalpha_o_tau_tau-dalpha_o_tau_tau_vector
	#Calculate Residual Terms
	residual=helmholtz_energy_residual(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	
	dalpha_tau=d_alpha_d_tau(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	
	dalpha_delta=d_alpha_d_delta(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	
	dalpha_tau_tau=d_alpha_d_tau_d_tau(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
        
        dalpha_delta_delta=d_alpha_d_delta_d_delta(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
        
        dalpha_delta_tau=d_alpha_d_delta_d_tau(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	
	
        #compressibility Density<-->Pressure
	Z=delta*dalpha_delta+float(1.0)
        #Enthalpy
	h1=R_methane*T*1.0
	h2=R_methane*T*tau*(dalpha_o+dalpha_tau)
	h3=R_methane*T*delta*dalpha_delta
	h=h1+h2+h3
	#Entropy
	s=R_methane*(tau*dalpha_o+tau*dalpha_tau-ideal-residual)
	#Absolute Helmholtz Energy
	absolute_helmholtz_energy=R_methane*T*(ideal+residual)
	#Isochoric heat capacity
	cv=-1.0*R_methane*(tau*tau*(dalpha_o_tau_tau+dalpha_tau_tau))
	#Isobaric heat capacity
	cp=cv+R_methane*pow((1.0+delta*dalpha_delta-delta*tau*dalpha_delta_tau),2)/(1.0+2.0*delta*dalpha_delta+delta*delta*dalpha_delta_delta)
	
	# Speed of Sound?
	
	a1 = 1 + delta*dalpha_delta - delta*tau*dalpha_delta_tau;
        b1 = tau*tau*(dalpha_o_tau_tau + dalpha_tau_tau);
        w = 1 + 2*delta*dalpha_delta + delta*delta*dalpha_delta_delta - a1*a1/b1;
        speed_o_sound=sqrt(R_methane*T*w*1000)
        
        return Z,h,s,absolute_helmholtz_energy,cv,cp,speed_o_sound,R_methane*T*ideal,residual*R_methane*T
