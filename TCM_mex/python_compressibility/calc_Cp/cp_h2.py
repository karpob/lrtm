from thermodynamic_functions.C_p_vector import C_p_vector
from thermodynamic_functions.C_p_vector_mix import C_p_vector_mix
from residual_pressure_pure_substance import residual_pressure_pure_substance
from scipy.optimize import fsolve
from helmholtz_functions.mBWR_to_Helmholtz import mBWR_to_Helmholtz
from numpy import asarray,zeros
from scale_tau_and_delta_kw import scale_tau_and_delta_kw
#f=open('output_T_P.txt','r')
#data=f.readline()
#f.close()
#T,P_h2,P_he,P_ch4,P_h2o=data.split("\t")
#bars to kPa
Bars_to_kPa=100
T=asarray([500,500])
P_h2=asarray([1,1])*Bars_to_kPa
P_he=asarray([1,1])*Bars_to_kPa
P_ch4=asarray([1,1])*Bars_to_kPa
P_h2o=asarray([1,1])*Bars_to_kPa

parameters_H2={'N_i':asarray([-6.93643,0.01,2.1101,4.52059,0.732564,-1.34086,0.130985,-0.777414,0.351944,-0.0211716,0.0226312,0.032187,-0.0231752,0.0557346]),
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
                             
mBWR_Coef=[0.4558980227431e-4,#1 
                   0.1260692007853e-2,#2
                   -0.7139657549318e-2,#3
                   0.9728903861441e-2,#4
                   -0.1589302471562e-1,#5
                   0.1454229259623e-5,#6
                   -0.4708238429298e-4,#7
                   0.1132915232587e-2,#8
                   0.2410763742104e-2,#9
                   -0.5093547838381e-8,#10
                   0.2699726927900e-5,#11
                   -0.3954146691114e-4,#12`
                   0.1551961438127e-8,#13
                   0.1050712335785e-7,#14
                   -0.5501158366750e-7,#15
                   -0.1037673478521e-9,#16
                   0.6446881346448e-12,#17
                   0.3298960057071e-10,#18
                   -0.3555585738784e-12,#19
                   -0.6885401367690e-2,#20
                   0.9166109232806e-2,#21
                   -0.6544314242937e-5,#22
                   -0.3315398880031e-4,#23
                   -0.2067693644676e-7,#24
                   0.3850153114958e-7,#25
                   -0.1399040626999e-10,#26
                   -0.1888462892389e-11,#27
                   -0.4595138561035e-14,#28
                   0.6872567403738e-14,#29
                   -0.6097223119177e-18,#30
                   -0.7636186157005e-17,#31
                   0.3848665703556e-17]#32
#note T, and density are extra inputs, totally unused..., make them 0's for now.. 
#helium uses R=8.314310                                               
[N_i_he,d_i_he,t_i_he,p_i_he]=mBWR_to_Helmholtz(mBWR_Coef,0,0,17.3990,5.19530,8.314310)
                           
parameters_he={'N_i':N_i_he,
               'd_i':d_i_he,
               't_i':t_i_he,
               'p_i':p_i_he,
               'M_amu':float(4.0026),
               'R':float(8.314310), 
               'To':4.230359714841141,  #Kelvin
               'rho_o':31.163394763964778, # mol/L  
               'ho':108.78863197310453,    # J/mol
               'so':3.6929233790579463 ,   # J/mol/K
               'Tc':5.19530,# Kelvin
               'rho_c':17.3990*4.0026, #(mol/L)-> kg/m^3
                'ni':float(2.5),
                'ti':float(0.0),
                'vi':asarray([0.0,0.0]),
                'ui':asarray([0.0,0.0]),
                'n_ideal_gas_terms_pow':1,
                'n_ideal_gas_terms_exp':0,
                'n_power_terms':int(80),
                'n_power_terms_wo_exp':int(32),
                'n_power_terms_w_exp':int(48),
                'n_gaussian_terms':int(0),
                'n_critical_terms':int(0),
                'phi_i':float(0),
                'Beta_i':float(0), 
                'gamma_i':float(0), 
                'D_i':float(0),
                'RES_a':float(0),
                'RES_b':float(0),
                'RES_B':float(0),
                'RES_C':float(0),
                'RES_D':float(0),
                'RES_A':float(0),
                'ideal_eqn_type':'Cp'}
                
parameters_ch4={'N_i':asarray([0.43679010280e-01,0.67092361990, -0.17655778590e01, 0.85823302410, -0.12065130520e01,0.51204672200, -0.40000107910e-03 ,-0.12478424230e-01,0.31002697010e-01,0.17547485220e-02,-0.31719216050e-05,-0.22403468400e-05,0.29470561560e-06,0.18304879090, 0.15118836790,-0.42893638770,0.68940024460e-01,-0.14083139960e-01,-0.30630548300e-01,-0.29699067080e-01,-0.19320408310e-01,-0.11057399590,0.99525489950e-01,0.85484378250e-02, -0.61505556620e-01,-0.42917924230e-01,-0.18132072900e-01,0.34459047600e-01,-0.23859194500e-02,-0.11590949390e-01,0.66416936020e-01,-0.23715495900e-01, -0.39616249050e-01, -0.13872920440e-01,0.33894895990e-01,-0.29273787530e-02,0.93247999460e-04,-0.62871715180e+01,0.12710694670e+02,-0.64239534660e+01]),
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

parameters_h2o={'N_i':asarray([0.12533547935523E-1,0.78957634722828E1, -0.87803203303561E1,0.31802509345418, -0.26145533859358, -0.78199751687981E-2,0.88089493102134E-2, -0.66856572307965,0.20433810950965,-0.66212605039687E-4,-0.19232721156002,-0.25709043003438,0.16074868486251,-0.40092828925807E-1,0.39343422603254E-6, -0.75941377088144E-5,0.56250979351888E-3, -0.15608652257135E-4,0.11537996422951E-8, 0.36582165144204E-6, -0.13251180074668E-11,-0.62639586912454E-9, -0.10793600908932, 0.17611491008752E-1, 0.22132295167546, -0.40247669763528,0.58083399985759, 0.49969146990806E-2, -0.31358700712549E-1,-0.74315929710341,0.47807329915480, 0.20527940895948E-1,-0.13636435110343, 0.14180634400617E-1, 0.83326504880713E-2, -0.29052336009585E-1,0.38615085574206E-1, -0.20393486513704E-1,-0.16554050063734E-2,0.19955571979541E-2, 0.15870308324157E-3, -0.16388568342530E-4, 0.43613615723811E-1, 0.34994005463765E-1, -0.76788197844621E-1,0.22446277332006E-1, -0.62689710414685E-4,-0.55711118565645E-9,-0.19905718354408, 0.31777497330738,-0.11841182425981,-0.31306260323435e2,0.31546140237781e2,-0.25213154341695e4,-0.14874640856724,0.31806110878444]),
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

#here we are using the ideal gas law as a proxy for density. The final pressure will be defined by the equation of state.
#In other words this is a kludge to keep density fixed, then calculate the real pressure.	                         

density_h2=P_h2/((parameters_H2['R']/parameters_H2['M_amu'])*T)
#density_h2=fsolve(residual_pressure_pure_substance,density_guess_h2,args=(P_h2,T,parameters_H2))
delta_h2=density_h2/parameters_H2['rho_c']

density_h2o=P_h2o/((parameters_h2o['R']/parameters_h2o['M_amu'])*T)
#density_h2o=fsolve(residual_pressure_pure_substance,density_guess_h2o,args=(P_h2o,T,parameters_h2o))
delta_h2o=density_h2o/parameters_h2o['rho_c']

density_he=P_he/((parameters_he['R']/parameters_he['M_amu'])*T)
#density_he=fsolve(residual_pressure_pure_substance,density_guess_he,args=(P_he,T,parameters_he))
delta_he=density_he/parameters_he['rho_c']

density_ch4=P_he/((parameters_ch4['R']/parameters_ch4['M_amu'])*T)
#density_ch4=fsolve(residual_pressure_pure_substance,density_guess_ch4,args=(P_ch4,T,parameters_ch4))
delta_ch4=density_ch4/parameters_ch4['rho_c']


#note we're calculating x_h2 and x_ch4 according to the mixture of just H2+CH4, not the whole thing He+h2o+etc...
molar_density_h2=density_h2/parameters_H2['M_amu']
molar_density_ch4=density_ch4/parameters_ch4['M_amu']
x_h2=molar_density_h2/(molar_density_h2+molar_density_ch4)
x_ch4=molar_density_ch4/(molar_density_ch4+molar_density_h2)


x_ch4=asarray([0.1,0.1])
x_h2=asarray([0.9,0.9])
M_amu_ch4=x_ch4*parameters_ch4['M_amu']+x_h2*parameters_H2['M_amu']
density_mol=asarray([0.2925151700533289,0.2925151700533289])
density=M_amu_ch4*density_mol
T=asarray([300,300])

N_i_ch4=asarray([-0.25157134971934,-0.62203841111983e-2,0.88850315184396e-1,-0.35592212573239e-1])
BetaT_ch4=1.0
BetaV_ch4=1.0
GammaT_ch4=1.352643115
GammaV_ch4=1.018702573
t_i_ch4=asarray([2.0,-1.0,1.75,1.4])
d_i_ch4=asarray([1.0,3.0,3.0,4.0])
p_i_ch4=asarray([0.0,0.0,0.0,0.0])                              
Rmix_ch4=x_ch4*parameters_ch4['R']+x_h2*parameters_H2['R']
M_amu_ch4=x_ch4*parameters_ch4['M_amu']+x_h2*parameters_H2['M_amu']                              
n_power_terms_wo_exp_ch4=4
n_power_terms_w_exp_ch4=0

h2_ch4_mix_params={'N_i':N_i_ch4,
                   'BetaT':float(1.0),
                   'BetaV':float(1.0),
                   'GammaT':float(1.352643115),
                   'GammaV':float(1.018702573),
                   't_i':t_i_ch4,
                   'd_i':d_i_ch4,
                   'p_i':p_i_ch4,
                   'n_power_terms_wo_exp':n_power_terms_wo_exp_ch4,
                   'n_power_terms_w_exp':n_power_terms_w_exp_ch4,
                   'phi_i':0,
                   'Beta_i':0,
                   'gamma_i':0,
                   'D_i':0,
                   'R':Rmix_ch4,
                   'M_amu':M_amu_ch4,
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

Cp_h2=C_p_vector(delta_h2,parameters_H2['Tc']/T,parameters_H2)
Cp_h2o=C_p_vector(delta_h2o,parameters_h2o['Tc']/T,parameters_h2o)
Cp_ch4=C_p_vector(delta_ch4,parameters_ch4['Tc']/T,parameters_ch4)
Cp_He=C_p_vector(delta_he,parameters_he['Tc']/T,parameters_he)
#density=molar_density_h2[0]*parameters_H2['M_amu']+molar_density_ch4[0]*parameters_ch4['M_amu']
[tau,delta,Tc,rho_c]=scale_tau_and_delta_kw(T,density/h2_ch4_mix_params['M_amu'],x_h2,x_ch4,parameters_H2['Tc'],parameters_ch4['Tc'],parameters_H2['rho_c']/parameters_H2['M_amu'],parameters_ch4['rho_c']/parameters_ch4['M_amu'],h2_ch4_mix_params['BetaT'],h2_ch4_mix_params['BetaV'],h2_ch4_mix_params['GammaT'],h2_ch4_mix_params['GammaV']) 
print density,T
CP_mix=C_p_vector_mix(density,T,delta,tau,parameters_H2,parameters_ch4,h2_ch4_mix_params,x_h2,x_ch4)
f=open('input_Cp.txt','w')
print density,x_h2[0],x_ch4[0],CP_mix[0]
f.write(str(Cp_h2[0])+'\t'+str(Cp_h2o[0])+'\t'+str(Cp_ch4[0])+'\t'+str(Cp_He[0]))
f.close()


