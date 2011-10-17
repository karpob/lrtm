def get_density(P,T):
        # this will take the equation of stat for water vapor, and solve for a water density at a given pressure and
        # temperature
        from scipy.optimize import fsolve
        from thermodynamic_functions.residual_pressure_pure_substance import residual_pressure_pure_substance
        from numpy import zeros,asarray
        bars_to_kPa=100.0
        kilograms_to_grams=1000.0
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
        h2o_gamma_i=zeros(51+3)
        h2o_gamma_i[51]=1.21
        h2o_gamma_i[52]=1.21
        h2o_gamma_i[53]=1.25 
        h2o_D_i=zeros(51+3)
        h2o_D_i[51]=1.0
        h2o_D_i[52]=1.0
        h2o_D_i[53]=1.0
        h2o_RES_a=zeros(51+5)
        h2o_RES_b=zeros(51+5)
        h2o_RES_B=zeros(51+5)
        h2o_RES_C=zeros(51+5)
        h2o_RES_D=zeros(51+5)
        h2o_RES_A=zeros(51+5)        
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
        h2o_params={'N_i':asarray([0.12533547935523E-1,0.78957634722828E1, -0.87803203303561E1,0.31802509345418, -0.26145533859358, -0.78199751687981E-2,0.88089493102134E-2, -0.66856572307965,0.20433810950965,-0.66212605039687E-4,-0.19232721156002,-0.25709043003438,0.16074868486251,-0.40092828925807E-1,0.39343422603254E-6, -0.75941377088144E-5,0.56250979351888E-3, -0.15608652257135E-4,0.11537996422951E-8, 0.36582165144204E-6, -0.13251180074668E-11,-0.62639586912454E-9, -0.10793600908932, 0.17611491008752E-1, 0.22132295167546, -0.40247669763528,0.58083399985759, 0.49969146990806E-2, -0.31358700712549E-1,-0.74315929710341,0.47807329915480, 0.20527940895948E-1,-0.13636435110343, 0.14180634400617E-1, 0.83326504880713E-2, -0.29052336009585E-1,0.38615085574206E-1, -0.20393486513704E-1,-0.16554050063734E-2,0.19955571979541E-2, 0.15870308324157E-3, -0.16388568342530E-4, 0.43613615723811E-1, 0.34994005463765E-1, -0.76788197844621E-1,0.22446277332006E-1, -0.62689710414685E-4,-0.55711118565645E-9,-0.19905718354408, 0.31777497330738,-0.11841182425981,-0.31306260323435e2,0.31546140237781e2,-0.25213154341695e4,-0.14874640856724,0.31806110878444]),
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
                     'phi_i':h2o_phi_i,
                     'Beta_i':h2o_Beta_i,
                     'gamma_i':h2o_gamma_i,
                     'D_i':h2o_D_i,
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
                     'RES_a':h2o_RES_a,
                     'RES_b':h2o_RES_b,
                     'RES_B':h2o_RES_B,
                     'RES_C':h2o_RES_C,
                     'RES_D':h2o_RES_D,
                     'RES_A':h2o_RES_A,
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
                   
        guess=(1e5/1e3)*P/((h2o_params['R']/h2o_params['M_amu'])*T)
       
        density=fsolve(residual_pressure_pure_substance,guess,args=(asarray([P*bars_to_kPa]),asarray([T]),h2o_params))
        #density is in terms of kg/m**3 return g/m**3
        return density*kilograms_to_grams
