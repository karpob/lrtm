def pressure_h2_h2o_mixture(params,T,P,density,x_h2,x_h2o):
        from pylab import zeros,shape
        from scipy.optimize import leastsq
        from residual_pressure import residual_pressure
        F_ij,beta_ij,phi_ij,sigma_ij,xi_ij=params
        
        # equations written to get values for a T (K), and density (kg/m^3)


        h2_N_i=[-6.93643,0.01,2.1101,4.52059,0.732564,-1.34086,0.130985,-0.777414,0.351944,-0.0211716,0.0226312,0.032187,-0.0231752,0.0557346];
        h2_t_i=[0.6844,1,0.989,0.489,0.803,1.1444,1.409,1.754,1.311,4.187,5.646,0.791,7.249,2.986]
        h2_d_i=[1,4,1,1,2,2,3,1,3,2,1,3,1,1]
        h2_p_i=[0,0,0,0,0,0,0,1,1,0/1,0/1,0/1,0/1,0/1]
        h2_phi_i=[0,0,0,0,0,0,0,0,0,1.685,0.489,0.103,2.506,1.607]# Note, we need to use -sign here
        h2_Beta_i=[0,0,0,0,0,0,0,0,0, 0.171,0.2245,0.1304,0.2785,0.3967] #Note, we need to use - sign here
        h2_gamma_i=[0,0,0,0,0,0,0,0,0,0.7164,1.3444,1.4517,0.7204,1.5445] 
        h2_D_i=[0,0,0,0,0,0,0,0,0,1.506,0.156,1.736,0.67,1.6620]
        h2_Tc=33.145
        h2_rho_c=15.508*2.01594 # mol/l -> kg/m^3 factors of 1000 liter converstion, and gram conversion cancel
        h2_R=8.314472
        h2_ni=2.5   
        h2_ti=0.0
        h2_vi=[1.616,-0.4117,-0.792,0.758,1.217] 
        h2_ui=[531,751,1989,2484,6859]
        h2_R=8.314472           #J/mol/K
        h2_ho = 7206.9069892047 #J/mol
        h2_so = 143.4846187346  #J/mol/K
        h2_n_power_terms=9
        h2_n_gaussian_terms=5
        h2_n_critical_terms=0
        h2_n_ideal_gas_terms_pow=0
        h2_n_ideal_gas_terms_exp=5
        h2_RES_a     = zeros(len(h2_ui))
        h2_RES_b     = zeros(len(h2_ui))
        h2_RES_B     = zeros(len(h2_ui))
        h2_RES_C     = zeros(len(h2_ui))
        h2_RES_D     = zeros(len(h2_ui))
        h2_RES_A     = zeros(len(h2_ui))
        h2_To=273.15        
        h2_rho_o=((0.001/(h2_R*h2_To))*1000)
        h2_M_amu=2.01594

        h2o_N_i=[0.12533547935523E-1,0.78957634722828E1, -0.87803203303561E1,
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
        h2o_t_i=[-0.5, 0.875,  1.0,  0.5,  0.75, 0.375,  1.0,  4.0,  6.0, 12.0,  1.0,
           5.0, 4.0  ,  2.0, 13.0,  9.0 , 3.0  ,  4.0, 11.0,  4.0, 13.0,  1.0,
           7.0, 1.0  ,  9.0, 10.0, 10.0 , 3.0  ,  7.0, 10.0, 10.0,  6.0, 10.0,
          10.0, 1.0  ,  2.0,  3.0,  4.0 , 8.0  ,  6.0,  9.0,  8.0, 16.0, 22.0,
          23.0,23.0  , 10.0, 50.0, 44.0, 46.0  , 50.0,  0.0,  1.0,  4.0]
        h2o_d_i=[1.0, 1.0, 1.0, 2.0,  2.0,  3.0,  4.0,  1.0,  1.0, 1.0, 2.0,  2.0,  3.0,  4.0,
          4.0, 5.0, 7.0, 9.0, 10.0, 11.0, 13.0, 15.0,  1.0, 2.0, 2.0,  2.0,  3.0,  4.0,
          4.0, 4.0, 5.0, 6.0,  6.0,  7.0,  9.0,  9.0,  9.0, 9.0, 9.0, 10.0, 10.0, 12.0,
          3.0, 4.0, 4.0, 5.0, 14.0,  3.0,  6.0,  6.0,  6.0, 3.0, 3.0,  3.0]
        h2o_p_i=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0,
          2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
          2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 4.0, 6.0, 6.0, 6.0, 6.0]
        h2o_phi_i=zeros(51+3)
        h2o_phi_i[51]=20.0 # Note, we need to use -sign here
        h2o_phi_i[52]=20.0 # Note, we need to use -sign here
        h2o_phi_i[53]=20.0 # Note, we need to use -sign here
        h2o_Beta_i=zeros(51+5)
        h2o_Beta_i[51]=150.0 #Note, we need to use - sign here
        h2o_Beta_i[52]=150.0 #Note, we need to use - sign here
        h2o_Beta_i[53]=250.0 #Note, we need to use - sign here
        h2o_Beta_i[54]=0.3  #Note, we need to use - sign here
        h2o_Beta_i[55]=0.3 #Note, we need to use - sign here
        h2o_gamma_i=zeros(len(h2o_p_i)+3)
        h2o_gamma_i[51]=1.21
        h2o_gamma_i[52]=1.21
        h2o_gamma_i[53]=1.25 
        h2o_D_i=zeros(51+3)
        h2o_D_i[51]=1.0
        h2o_D_i[52]=1.0
        h2o_D_i[53]=1.0 
        h2o_R=8.314472
        h2o_Tc  = 647.096     #K
        h2o_rho_c      =322.0#18.015268*17.8737279956 #322.000214485       #kg m^-3 #mol/liter 17.8737279956
        h2o_R_h2o = 0.46151805#R/18.015268#  #kJ kg^-1 K^-1
        h2o_M_amu=h2o_R/h2o_R_h2o
	
	
        h2o_ideal_n = [-8.32044648201,6.6832105268,3.00632,0.012436,0.97315,1.27950,0.96956,0.24873]
        h2o_ideal_gamma = [0.0,0.0,0.0,1.28728967,3.53734222,7.74073708,9.24437796,27.5075105]
        
	
        h2o_n_power_terms=51#len(p_i)
        h2o_n_gaussian_terms=3
        h2o_n_critical_terms=2
        #print delta,tau
        #critcal terms
        h2o_RES_a=zeros(len(h2o_p_i)+5)
        h2o_RES_b=zeros(len(h2o_p_i)+5)
        h2o_RES_B=zeros(len(h2o_p_i)+5)
        h2o_RES_C=zeros(len(h2o_p_i)+5)
        h2o_RES_D=zeros(len(h2o_p_i)+5)
        h2o_RES_A=zeros(len(h2o_p_i)+5)        
        h2o_RES_a[54]=3.5
        h2o_RES_a[55]=3.5
        h2o_RES_b[54]=0.85
        h2o_RES_b[55]=0.95
        h2o_RES_B[54]=0.2
        h2o_RES_B[55]=0.2
        h2o_RES_C[54]=28.0
        h2o_RES_C[55]=32.0
        h2o_RES_D[54]=700.0
        h2o_RES_D[55]=800.0
        h2o_RES_A[54]=0.32
        h2o_RES_A[55]=0.32

        h2o_ni=0
        h2o_ti=0
        h2o_vi=0
        h2o_ui=0
        h2o_ho=0
        h2o_so=0
        h2o_n_ideal_gas_terms_pow=0
        h2o_n_ideal_gas_terms_exp=0
        h2o_To=0
        h2o_Po=0
        h2o_rho_o=0
        
        
        h2o_h2_mix_N_i=[-0.245476271425e-1,
         -0.241206117483,
         -0.513801950309e-2,
         -0.239824834123e-1,
          0.259772344009,
         -0.172014123104,
          0.429490028551e-1,
          -0.202108593862e-3,
          -0.382984234857e-2,
          0.262992331354e-5]
        #KLUUUUUDGE
        for i in range(0,len(h2o_h2_mix_N_i)):
                h2o_h2_mix_N_i[i]=F_ij*h2o_h2_mix_N_i[i]
        h2o_h2_mix_t_i=[2,4,-2,1,4,4,4,0,4,-2]
        h2o_h2_mix_d_i=[1,1,1,2,3,4,5,6,6,8]
        h2o_h2_mix_p_i=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0.0]
        h2o_h2_mix_phi_i=0
        h2o_h2_mix_Beta_i=0
        h2o_h2_mix_gamma_i=0
        h2o_h2_mix_D_i=0
        h2o_h2_mix_R=8.314472
#h2o_h2_mix_Tc  = 647.096     #K
#h2o_h2_mix_rho_c      =322.0#18.015268*17.8737279956 #322.000214485       #kg m^-3 #mol/liter 17.8737279956
#h2o_h2_mix_R_h2o = 0.46151805#R/18.015268#  #kJ kg^-1 K^-1

# Get the mixture's equivalent Mass/amu
        h2o_h2_mix_M_amu=x_h2o*(h2o_R/h2o_R_h2o)+x_h2*h2_M_amu
	
	
        h2o_h2_mix_ideal_n = 0
        h2o_h2_mix_ideal_gamma = 0
        
	
        h2o_h2_mix_n_power_terms=10
        h2o_h2_mix_n_gaussian_terms=0
        h2o_h2_mix_n_critical_terms=0

        h2o_h2_mix_RES_a=0
        h2o_h2_mix_RES_b=0
        h2o_h2_mix_RES_B=0
        h2o_h2_mix_RES_C=0
        h2o_h2_mix_RES_D=0
        h2o_h2_mix_RES_A=0        


        h2o_h2_mix_ni=0
        h2o_h2_mix_ti=0
        h2o_h2_mix_vi=0
        h2o_h2_mix_ui=0
        h2o_h2_mix_ho=0
        h2o_h2_mix_so=0
        h2o_h2_mix_n_ideal_gas_terms_pow=0
        h2o_h2_mix_n_ideal_gas_terms_exp=0
        h2o_h2_mix_To=0
        h2o_h2_mix_Po=0
        h2o_h2_mix_rho_o=0
        
        
                
      
        Rmix=h2o_h2_mix_R/((x_h2o*h2o_M_amu)+(x_h2*h2_M_amu))
        
        
        #Here we need to iteratively solve (using total Pressure to calculate error) for the right density for a given set of mixture parameters, mole fractions, and total pressure
        #leastsq usage plsq=leastsq(residuals,p0,args=(y_meas,x))
        # Kludge! too lazy to change variable list, so I'm setting h2o_h2_mix_Tc and rho to 0
        
        
        h2_delta=(h2_M_amu/((x_h2o*h2o_M_amu)+(x_h2*h2_M_amu)))*(x_h2*density)/h2_rho_c # Need to scale by molecular weight
        h2_tau=h2_Tc/T
        
        h2o_delta=(h2o_M_amu/((x_h2o*h2o_M_amu)+(x_h2*h2_M_amu)))*(x_h2o*density)/h2o_rho_c #Need to scale by molecular weight
        h2o_tau=h2o_Tc/T
       
        error_in_pressure=residual_pressure(density,Pin,h2_N_i,h2_t_i,h2_d_i,h2_p_i,h2_phi_i,h2_Beta_i,h2_gamma_i,h2_D_i,h2_n_power_terms,h2_n_gaussian_terms,h2_n_critical_terms,h2_RES_a,h2_RES_b,h2_RES_B,h2_RES_C,h2_RES_D,h2_RES_A,h2_R,h2_Tc,h2_rho_c,h2_M_amu,h2o_N_i,h2o_t_i,h2o_d_i,h2o_p_i,h2o_phi_i,h2o_Beta_i,h2o_gamma_i,h2o_D_i,h2o_n_power_terms,h2o_n_gaussian_terms,h2o_n_critical_terms,h2o_RES_a,h2o_RES_b,h2o_RES_B,h2o_RES_C,h2o_RES_D,h2o_RES_A,h2o_R,h2o_Tc,h2o_rho_c,h2o_M_amu,h2o_h2_mix_N_i,h2o_h2_mix_t_i,h2o_h2_mix_d_i,h2o_h2_mix_p_i,h2o_h2_mix_phi_i,h2o_h2_mix_Beta_i,h2o_h2_mix_gamma_i,h2o_h2_mix_D_i,h2o_h2_mix_n_power_terms,h2o_h2_mix_n_gaussian_terms,h2o_h2_mix_n_critical_terms,h2o_h2_mix_RES_a,h2o_h2_mix_RES_b,h2o_h2_mix_RES_B,h2o_h2_mix_RES_C,h2o_h2_mix_RES_D,h2o_h2_mix_RES_A,h2o_h2_mix_R,h2o_h2_mix_Tc,h2o_h2_mix_rho_c,h2o_h2_mix_M_amu,beta_ij,phi_ij,sigma_ij,xi_ij,x_h2,x_h2o,T
              
        return error_in_pressure
        #
        #enthalpy_h2=enthalpy(h2_tau,h2_delta,h2_ni,h2_ti,h2_vi,h2_ui,h2_R,h2_ho,h2_so,h2_rho_c,h2_Tc,h2_n_ideal_gas_terms_pow,h2_n_ideal_gas_terms_exp,h2_To,h2_rho_o,h2_N_i,h2_t_i,h2_d_i,h2_p_i,h2_phi_i,h2_Beta_i,h2_gamma_i,h2_D_i,h2_n_power_terms,h2_n_gaussian_terms,h2_n_critical_terms,h2_RES_a,h2_RES_b,h2_RES_B,h2_RES_C,h2_RES_D,h2_RES_A,h2_M_amu,0,0,'Cp')
        #enthalpy_h2o=enthalpy(h2o_tau,h2o_delta,h2o_ni,h2o_ti,h2o_vi,h2o_ui,h2o_R,h2o_ho,h2o_so,h2o_rho_c,h2o_Tc,h2o_n_ideal_gas_terms_pow,h2o_n_ideal_gas_terms_exp,h2o_To,h2o_rho_o,h2o_N_i,h2o_t_i,h2o_d_i,h2o_p_i,h2o_phi_i,h2o_Beta_i,h2o_gamma_i,h2o_D_i,h2o_n_power_terms,h2o_n_gaussian_terms,h2o_n_critical_terms,h2o_RES_a,h2o_RES_b,h2o_RES_B,h2o_RES_C,h2o_RES_D,h2o_RES_A,h2o_M_amu,h2o_ideal_n,h2o_ideal_gamma,'Coef')
        #h_idmix=x_h2o*enthalpy_h2o+x_h2*enthalpy_h2
