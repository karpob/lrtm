def enthalpy_h2o(T,density):
	from helmholtz_functions.ideal_helmholtz_from_coef import ideal_helmholtz_from_coef
	from helmholtz_functions.ideal_helmholtz_from_coef_dtau import ideal_helmholtz_from_coef_dtau
	from helmholtz_functions.ideal_helmholtz_from_coef_dtau_dtau import ideal_helmholtz_from_coef_dtau_dtau
	from helmholtz_functions.d_alpha_d_tau import d_alpha_d_tau 
	from helmholtz_functions.d_alpha_d_delta import d_alpha_d_delta
        from helmholtz_functions.helmholtz_energy_residual import helmholtz_energy_residual
        from helmholtz_functions.d_alpha_d_tau_d_tau import d_alpha_d_tau_d_tau
        from helmholtz_functions.d_alpha_d_delta_d_delta import d_alpha_d_delta_d_delta
        from helmholtz_functions.d_alpha_d_delta_d_tau import d_alpha_d_delta_d_tau
        
        from helmholtz_functions.d_alpha_d_tau_vector import d_alpha_d_tau_vector 
	from helmholtz_functions.d_alpha_d_delta_vector import d_alpha_d_delta_vector
        from helmholtz_functions.ideal_helmholtz_from_coef_dtau_vector import ideal_helmholtz_from_coef_dtau_vector
        from thermodynamic_functions.enthalpy_vector import enthalpy_vector
        from helmholtz_functions.d_alpha_d_delta_d_delta_vector import d_alpha_d_delta_d_delta_vector
        from helmholtz_functions.d_alpha_d_delta_vector import d_alpha_d_delta_vector
        from helmholtz_functions.d_alpha_d_tau_d_tau_vector import d_alpha_d_tau_d_tau_vector
        from helmholtz_functions.d_alpha_d_delta_d_tau_vector import d_alpha_d_delta_d_tau_vector
        from helmholtz_functions.ideal_helmholtz_from_coef_vector import ideal_helmholtz_from_coef_vector
        from helmholtz_functions.ideal_helmholtz_from_coef_dtau_dtau_vector import ideal_helmholtz_from_coef_dtau_dtau_vector
        from helmholtz_functions.helmholtz_energy_residual_vector import helmholtz_energy_residual_vector
        from pylab import zeros
        from math import sqrt
        from numpy import asarray
	N_i=[0.12533547935523E-1,0.78957634722828E1, -0.87803203303561E1,
		0.31802509345418, 
		-0.26145533859358, 
		-0.78199751687981E-2,
		0.88089493102134E-2, 
		-0.66856572307965,
		0.20433810950965,
		-0.66212605039687E-4,
		-0.19232721156002,
		-0.25709043003438,
		0.16074868486251,
		-0.40092828925807E-1,
		0.39343422603254E-6, 
		-0.75941377088144E-5,
		0.56250979351888E-3, 
		-0.15608652257135E-4,
		0.11537996422951E-8, 
		0.36582165144204E-6, 
		-0.13251180074668E-11,
		-0.62639586912454E-9, 
		-0.10793600908932, 
		0.17611491008752E-1, 
		0.22132295167546, 
		-0.40247669763528,
		0.58083399985759, 
		0.49969146990806E-2, 
		-0.31358700712549E-1,
		-0.74315929710341,
		0.47807329915480, 
		0.20527940895948E-1,
		-0.13636435110343, 
		0.14180634400617E-1, 
		0.83326504880713E-2, 
		-0.29052336009585E-1,
		0.38615085574206E-1, 
		-0.20393486513704E-1,
		-0.16554050063734E-2,
		0.19955571979541E-2, 
		0.15870308324157E-3, 
		-0.16388568342530E-4, 
		0.43613615723811E-1, 
		0.34994005463765E-1, 
		-0.76788197844621E-1,
		0.22446277332006E-1, 
		-0.62689710414685E-4,
		-0.55711118565645E-9,
		-0.19905718354408, 
		0.31777497330738,
		-0.11841182425981,
		-0.31306260323435e2,
		0.31546140237781e2,
		-0.25213154341695e4,
		-0.14874640856724,
		0.31806110878444]
	t_i=[-0.5, 0.875,  1.0,  0.5,  0.75, 0.375,  1.0,  4.0,  6.0, 12.0,  1.0,
           5.0, 4.0  ,  2.0, 13.0,  9.0 , 3.0  ,  4.0, 11.0,  4.0, 13.0,  1.0,
           7.0, 1.0  ,  9.0, 10.0, 10.0 , 3.0  ,  7.0, 10.0, 10.0,  6.0, 10.0,
          10.0, 1.0  ,  2.0,  3.0,  4.0 , 8.0  ,  6.0,  9.0,  8.0, 16.0, 22.0,
          23.0,23.0  , 10.0, 50.0, 44.0, 46.0  , 50.0,  0.0,  1.0,  4.0]
	d_i=[1.0, 1.0, 1.0, 2.0,  2.0,  3.0,  4.0,  1.0,  1.0, 1.0, 2.0,  2.0,  3.0,  4.0,
          4.0, 5.0, 7.0, 9.0, 10.0, 11.0, 13.0, 15.0,  1.0, 2.0, 2.0,  2.0,  3.0,  4.0,
          4.0, 4.0, 5.0, 6.0,  6.0,  7.0,  9.0,  9.0,  9.0, 9.0, 9.0, 10.0, 10.0, 12.0,
          3.0, 4.0, 4.0, 5.0, 14.0,  3.0,  6.0,  6.0,  6.0, 3.0, 3.0,  3.0]
	p_i=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0,
          2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
          2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 4.0, 6.0, 6.0, 6.0, 6.0]
        phi_i=zeros(51+3)
	phi_i[51]=20.0 # Note, we need to use -sign here
	phi_i[52]=20.0 # Note, we need to use -sign here
	phi_i[53]=20.0 # Note, we need to use -sign here
	Beta_i=zeros(51+5)
	Beta_i[51]=150.0 #Note, we need to use - sign here
	Beta_i[52]=150.0 #Note, we need to use - sign here
	Beta_i[53]=250.0 #Note, we need to use - sign here
	Beta_i[54]=0.3  #Note, we need to use - sign here
	Beta_i[55]=0.3 #Note, we need to use - sign here
	gamma_i=zeros(len(p_i)+3)
	gamma_i[51]=1.21
	gamma_i[52]=1.21
	gamma_i[53]=1.25 
	D_i=zeros(51+3)
	D_i[51]=1.0
	D_i[52]=1.0
	D_i[53]=1.0 
	
        R=8.314472
	Tc  = 647.096     #K
        rho_c      =322.0#18.015268*17.8737279956 #322.000214485       #kg m^-3 #mol/liter 17.8737279956
        R_h2o = 0.46151805#R/18.015268#  #kJ kg^-1 K^-1
	
	delta=density/rho_c
	tau=Tc/T
	
	ideal_n = [-8.32044648201,6.6832105268,3.00632,0.012436,0.97315,1.27950,0.96956,0.24873]
        ideal_gamma = [0.0,0.0,0.0,1.28728967,3.53734222,7.74073708,9.24437796,27.5075105]
        
	
	n_power_terms=51#len(p_i)
	n_gaussian_terms=3
	n_critical_terms=2
	#print delta,tau
	#critcal terms
	RES_a=zeros(len(p_i)+5)
	RES_b=zeros(len(p_i)+5)
	RES_B=zeros(len(p_i)+5)
	RES_C=zeros(len(p_i)+5)
	RES_D=zeros(len(p_i)+5)
	RES_A=zeros(len(p_i)+5)        
        RES_a[54]=3.5
        RES_a[55]=3.5
        RES_b[54]=0.85
        RES_b[55]=0.95
        RES_B[54]=0.2
        RES_B[55]=0.2
        RES_C[54]=28.0
        RES_C[55]=32.0
        RES_D[54]=700.0
        RES_D[55]=800.0
        RES_A[54]=0.32
        RES_A[55]=0.32
        h2o_component_params={'N_i':asarray([0.12533547935523E-1,0.78957634722828E1, -0.87803203303561E1,0.31802509345418, -0.26145533859358, -0.78199751687981E-2,0.88089493102134E-2, -0.66856572307965,0.20433810950965,-0.66212605039687E-4,-0.19232721156002,-0.25709043003438,0.16074868486251,-0.40092828925807E-1,0.39343422603254E-6, -0.75941377088144E-5,0.56250979351888E-3, -0.15608652257135E-4,0.11537996422951E-8, 0.36582165144204E-6, -0.13251180074668E-11,-0.62639586912454E-9, -0.10793600908932, 0.17611491008752E-1, 0.22132295167546, -0.40247669763528,0.58083399985759, 0.49969146990806E-2, -0.31358700712549E-1,-0.74315929710341,0.47807329915480, 0.20527940895948E-1,-0.13636435110343, 0.14180634400617E-1, 0.83326504880713E-2, -0.29052336009585E-1,0.38615085574206E-1, -0.20393486513704E-1,-0.16554050063734E-2,0.19955571979541E-2, 0.15870308324157E-3, -0.16388568342530E-4, 0.43613615723811E-1, 0.34994005463765E-1, -0.76788197844621E-1,0.22446277332006E-1, -0.62689710414685E-4,-0.55711118565645E-9,-0.19905718354408, 0.31777497330738,-0.11841182425981,-0.31306260323435e2,0.31546140237781e2,-0.25213154341695e4,-0.14874640856724,0.31806110878444]),
                    't_i':asarray([-0.5, 0.875,  1.0,  0.5,  0.75, 0.375,  1.0,  4.0,  6.0, 12.0,  1.0,
                                  5.0, 4.0  ,  2.0, 13.0,  9.0 , 3.0  ,  4.0, 11.0,  4.0, 13.0,  1.0,
                                  7.0, 1.0  ,  9.0, 10.0, 10.0 , 3.0  ,  7.0, 10.0, 10.0,  6.0, 10.0,
                                  10.0, 1.0  ,  2.0,  3.0,  4.0 , 8.0  ,  6.0,  9.0,  8.0, 16.0, 22.0,
                                  23.0,23.0  , 10.0, 50.0, 44.0, 46.0  , 50.0,  0.0,  1.0,  4.0]),
                     'd_i':asarray([1.0, 1.0, 1.0, 2.0,  2.0,  3.0,  4.0,  1.0,  1.0, 1.0, 2.0,  2.0,  3.0,  4.0,
                                  4.0, 5.0, 7.0, 9.0, 10.0, 11.0, 13.0, 15.0,  1.0, 2.0, 2.0,  2.0,  3.0,  4.0,
                                  4.0, 4.0, 5.0, 6.0,  6.0,  7.0,  9.0,  9.0,  9.0, 9.0, 9.0, 10.0, 10.0, 12.0,
                                  3.0, 4.0, 4.0, 5.0, 14.0,  3.0,  6.0,  6.0,  6.0, 3.0, 3.0,  3.0]),
                     'p_i':asarray([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                                  1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0,
                                  2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
                                  2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 4.0, 6.0, 6.0, 6.0, 6.0]),
                     'phi_i':phi_i,
                     'Beta_i':Beta_i,
                     'gamma_i':gamma_i,
                     'D_i':D_i,
                     'R':float(8.314472),
                     'Tc':float(647.096),     #K
                     'rho_c':float(322.0),#18.015268*17.8737279956 #322.000214485       #kg m^-3 #mol/liter 17.8737279956
                     'M_amu':float(8.314472/0.46151805),
	             'ideal_n': asarray([-8.32044648201,6.6832105268,3.00632,0.012436,0.97315,1.27950,0.96956,0.24873]),
                     'ideal_gamma':asarray([0.0,0.0,0.0,1.28728967,3.53734222,7.74073708,9.24437796,27.5075105]),
                     'n_power_terms':int(51),#len(p_i)
                     'n_power_terms_wo_exp':int(7),
                     'n_power_terms_w_exp':int(51-7),
                     'n_gaussian_terms':int(3),
                     'n_critical_terms':int(2),
                     'RES_a':RES_a,
                     'RES_b':RES_b,
                     'RES_B':RES_B,
                     'RES_C':RES_C,
                     'RES_D':RES_D,
                     'RES_A':RES_A,
                     'ni':int(0),
                     'ti':int(0),
                     'vi':int(0),
                     'ui':int(0),
                     'ho':int(0),
                     'so':int(0),
                     'n_ideal_gas_terms_pow':int(0),
                     'n_ideal_gas_terms_exp':int(0),
                     'To':int(0),
                     'Po':int(0),
                     'rho_o':int(0),
                     'ideal_eqn_type':'Coef'}
	#Calculate Ideal Terms
	delta_v=asarray([delta,delta,delta])
        tau_v=asarray([tau,tau,tau])
        
	ideal=ideal_helmholtz_from_coef(delta,tau,ideal_n,ideal_gamma)
        ideal_v=ideal_helmholtz_from_coef_vector(delta_v,tau_v,h2o_component_params['ideal_n'],h2o_component_params['ideal_gamma'])
        #print ideal-ideal_v
        
        dalpha_o=ideal_helmholtz_from_coef_dtau(delta,tau,ideal_n,ideal_gamma)
	
	
	dalpha_o_vector=ideal_helmholtz_from_coef_dtau_vector(delta_v,tau_v,h2o_component_params['ideal_n'],h2o_component_params['ideal_gamma'])
	#print dalpha_o-dalpha_o_vector
	dalpha_o_tau_tau=ideal_helmholtz_from_coef_dtau_dtau(delta,tau,ideal_n,ideal_gamma)
	dalpha_o_tau_tau_v=ideal_helmholtz_from_coef_dtau_dtau_vector(delta_v,tau_v,h2o_component_params['ideal_n'],h2o_component_params['ideal_gamma'])
	#print dalpha_o_tau_tau-dalpha_o_tau_tau_v
	
	#Calculate Residual Terms
	residual=helmholtz_energy_residual(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	residual_v=helmholtz_energy_residual_vector(tau_v,delta_v,h2o_component_params)
	print residual-residual_v
	dalpha_tau=d_alpha_d_tau(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	
	dalpha_tau_vector=d_alpha_d_tau_vector(tau_v,delta_v,h2o_component_params)
	#print dalpha_tau-dalpha_tau_vector
	dalpha_delta=d_alpha_d_delta(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	dalpha_delta_vector=d_alpha_d_delta_vector(tau_v,delta_v,h2o_component_params)
	#print dalpha_delta_vector-dalpha_delta
	
	dalpha_tau_tau=d_alpha_d_tau_d_tau(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
        dalpha_tau_tau_vector=d_alpha_d_tau_d_tau_vector(tau_v,delta_v,h2o_component_params)
        #print dalpha_tau_tau-dalpha_tau_tau_vector
        dalpha_delta_delta=d_alpha_d_delta_d_delta(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
        dalpha_delta_delta_vector=d_alpha_d_delta_d_delta_vector(tau_v,delta_v,h2o_component_params)
        #print T, dalpha_delta_delta-dalpha_delta_delta_vector
      #  print dalpha_delta_delta
        dalpha_delta_tau=d_alpha_d_delta_d_tau(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	dalpha_delta_tau_vector=d_alpha_d_delta_d_tau_vector(tau_v,delta_v,h2o_component_params)
	#print dalpha_delta_tau_vector-dalpha_delta_tau
        #compressibility Density<-->Pressure
	Z=delta*dalpha_delta+float(1.0)
        #Enthalpy
	h1=R_h2o*T*1.0
	h2=R_h2o*T*tau*(dalpha_o+dalpha_tau)
	h3=R_h2o*T*delta*dalpha_delta
	h=h1+h2+h3
	h_v=enthalpy_vector(tau_v,delta_v,Tc,rho_c,'pure_substance',h2o_component_params,h2o_component_params,h2o_component_params,0,0)
	#print dalpha_o,dalpha_o_vector
	#print h,R_h2o*T*h_v
	#Entropy
	
	s=R_h2o*(tau*dalpha_o+tau*dalpha_tau-ideal-residual)
	#Absolute Helmholtz Energy
	absolute_helmholtz_energy=R_h2o*T*(ideal+residual)
	#Isochoric heat capacity
	cv=-1.0*R_h2o*(tau*tau*(dalpha_o_tau_tau+dalpha_tau_tau))
	#Isobaric heat capacity
	cp=cv+R_h2o*pow((1.0+delta*dalpha_delta-delta*tau*dalpha_delta_tau),2)/(1.0+2.0*delta*dalpha_delta+delta*delta*dalpha_delta_delta)
	
	#print dalpha_delta,dalpha_delta_delta,dalpha_delta_tau,dalpha_tau_tau
	# Speed of Sound?
	#print dalpha_delta_delta,dalpha_delta_tau
	#speed_o_sound_p1 = 1.0 + 2.0*delta*dalpha_delta + delta*delta*dalpha_delta_delta
        #speed_o_sound_p2 = pow((1.0 + delta*dalpha_delta - delta*tau*dalpha_delta_tau),2)
        #speed_o_sound_p3 = tau*tau*(dalpha_o_tau_tau + dalpha_tau_tau)
        #speed_o_sound= sqrt((R_h2o*T)*(speed_o_sound_p1 - speed_o_sound_p2/speed_o_sound_p3)*1000.0)
        #what it's probably not 
        #d_tau_tau
        #d_delta
        #d_o_tau_tau
        a1 = 1 + delta*dalpha_delta - delta*tau*dalpha_delta_tau;
        b1 = tau*tau*(dalpha_o_tau_tau + dalpha_tau_tau);
        w = 1 + 2*delta*dalpha_delta + delta*delta*dalpha_delta_delta - a1*a1/b1;
        speed_o_sound=sqrt(R_h2o*T*w*1000)
        return Z,h,s,absolute_helmholtz_energy,cv,cp,speed_o_sound
