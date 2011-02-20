def residual_enthalpy_data(p,actual_enthalpy,P,T,x_h2,x_h2o):
        #from numpy import shape
        from numpy import log10,abs,asarray,shape,zeros
        from thermodynamic_functions.Excess_Enthalpy_vector import Excess_Enthalpy_vector        
        from scipy.optimize import leastsq,fsolve
        from residual_pressure import residual_pressure
        from scale_tau_and_delta_kw import scale_tau_and_delta_kw

        h2_params={'N_i':asarray([-6.93643,0.01,2.1101,4.52059,0.732564,-1.34086,0.130985,-0.777414,0.351944,-0.0211716,0.0226312,0.032187,-0.0231752,0.0557346]),
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
               
        BetaT=p[0]#p[0]#0.35027854              
        BetaV=p[1] # p[1]#-0.61593331
        GammaT=p[2]#p[2]#1.4401281         
        GammaV=p[3]#p[3]#-0.12750122
              
           
   
        N_i=asarray([p[4],p[5],p[6],p[7]])#asarray([7.92247566, -7.37889621,48.18464436,-9.3777652])
        t_i=asarray([p[8],p[9],p[10],p[11]])#asarray([2,-1,1.75,2])
        d_i=asarray([p[12],p[13],p[14],p[15]])#asarray([1,3,3,4])#[1,1,3,2])              
        #N_i=asarray([p[4],p[7],p[10],p[13]])#asarray([p[4],p[7],p[10],p[13]])
        #t_i=asarray([p[5],p[8],p[11],p[14]])#asarray([p[5],p[8],p[11],p[14]])
        #d_i=asarray([p[6],p[9],p[12],p[15]])#asarray([p[6],p[9],p[12],p[15]])
        p_i=asarray([p[16],p[17],p[18],p[19]])#asarray([0.0,0.0,0.0,0.0])                                                 
        #print p,shape(N_i),shape(t_i),shape(d_i),shape(p_i)
        
        Rmix=x_h2o*h2o_params['R']+x_h2*h2_params['R']
        M_amu=x_h2o*h2o_params['M_amu']+x_h2*h2_params['M_amu']                              
        
        n_power_terms_wo_exp=0
        n_power_terms_w_exp=4
        
        mix_params={'N_i':N_i,
                   'BetaT':BetaT,
                   'BetaV':BetaV,
                   'GammaT':GammaT,
                   'GammaV':GammaV,
                   't_i':t_i,
                   'd_i':d_i,
                   'p_i':p_i,
                   'n_power_terms_wo_exp':n_power_terms_wo_exp,
                   'n_power_terms_w_exp':n_power_terms_w_exp,
                   'phi_i':0,
                   'Beta_i':0,
                   'gamma_i':0,
                   'D_i':0,
                   'R':Rmix,
                   'M_amu':M_amu,
                   'ideal_n': 0,
                   'ideal_gamma': 0,
                   'n_gaussian_terms':0,
                   'n_critical_terms':0,
                   'RES_a':0,
                   'RES_b':0,
                   'RES_B':0,
                   'RES_C':0,
                   'RES_D':0,
                   'RES_A':0,     
                   'ni':0,
                   'ti':0,
                   'vi':0,
                   'ui':0,
                   'ho':0,
                   'so':0,
                   'n_ideal_gas_terms_pow':0,
                   'n_ideal_gas_terms_exp':0,
                   'To':0,
                   'Po':0,
                   'rho_o':0}
        #density_guess=P/((Rmix/M_amu)*T)
        density_guess=P/((Rmix/M_amu)*T)
        density=fsolve(residual_pressure,density_guess,args=(P,T,h2_params,h2o_params,mix_params,x_h2,x_h2o))           
        #[density,message]=leastsq(residual_pressure,density_guess,args=(P,T,h2_params,h2o_params,mix_params,x_h2,x_h2o))
        [tau,delta,Tc,rho_c]=scale_tau_and_delta_kw(T,density/mix_params['M_amu'],x_h2,x_h2o,h2_params['Tc'],h2o_params['Tc'],h2_params['rho_c']/h2_params['M_amu'],h2o_params['rho_c']/h2o_params['M_amu'],mix_params['BetaT'],mix_params['BetaV'],mix_params['GammaT'],mix_params['GammaV'])                                 
        calculated_enthalpy=Excess_Enthalpy_vector(density,T,tau,delta,Tc,rho_c,h2_params,h2o_params,mix_params,x_h2,x_h2o)
        
        err=actual_enthalpy-calculated_enthalpy
        #print err
        return err
