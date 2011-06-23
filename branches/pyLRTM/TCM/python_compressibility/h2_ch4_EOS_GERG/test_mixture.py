from pylab import *
from numpy import *
from thermodynamic_functions.pressure_vector import pressure_vector
from thermodynamic_functions.Excess_Enthalpy_vector import Excess_Enthalpy_vector
from thermodynamic_functions.enthalpy_vector import enthalpy_vector
from thermodynamic_functions.Virial_B_vector import Virial_B_vector
from thermodynamic_functions.Virial_C_vector import Virial_C_vector
from scale_tau_and_delta_kw import scale_tau_and_delta_kw


#test_data=asarray([[0.100000000000,0.900000000000,300.000000000,0.732404627240,1.00000000000,8741.44012790,28.8280156048],
#                  [0.100000000000,0.900000000000,300.000000000,0.732404627240,1.00000000000,2557.00393313,28.8280156048]])
x_h2=test_data[:,1]
x_ch4=test_data[:,0]
T=test_data[:,2]
density=test_data[:,4]
P_refprop=test_data[:,3]
enthalpy_refprop=test_data[:,5]
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





ch4_params={'N_i':asarray([0.43679010280e-01,0.67092361990, -0.17655778590e01, 0.85823302410, -0.12065130520e01,0.51204672200, -0.40000107910e-03 ,-0.12478424230e-01,0.31002697010e-01,0.17547485220e-02,-0.31719216050e-05,-0.22403468400e-05,0.29470561560e-06,0.18304879090, 0.15118836790,-0.42893638770,0.68940024460e-01,-0.14083139960e-01,-0.30630548300e-01,-0.29699067080e-01,-0.19320408310e-01,-0.11057399590,0.99525489950e-01,0.85484378250e-02, -0.61505556620e-01,-0.42917924230e-01,-0.18132072900e-01,0.34459047600e-01,-0.23859194500e-02,-0.11590949390e-01,0.66416936020e-01,-0.23715495900e-01, -0.39616249050e-01, -0.13872920440e-01,0.33894895990e-01,-0.29273787530e-02,0.93247999460e-04,-0.62871715180e+01,0.12710694670e+02,-0.64239534660e+01]),
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
                             

N_i=asarray([-0.25157134971934,-0.62203841111983e-2,0.88850315184396e-1,-0.35592212573239e-1])
BetaT=1.0
BetaV=1.0
GammaT=1.352643115
GammaV=1.018702573
t_i=asarray([2.0,-1.0,1.75,1.4])
d_i=asarray([1.0,3.0,3.0,4.0])
p_i=asarray([0.0,0.0,0.0,0.0])                              
Rmix=x_ch4*ch4_params['R']+x_h2*h2_params['R']
M_amu=x_ch4*ch4_params['M_amu']+x_h2*h2_params['M_amu']                              
n_power_terms_wo_exp=4
n_power_terms_w_exp=0

h2_ch4_mix_params={'N_i':N_i,
                   'BetaT':float(1.0),
                   'BetaV':float(1.0),
                   'GammaT':float(1.352643115),
                   'GammaV':float(1.018702573),
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
                   
[tau,delta,Tc,rho_c]=scale_tau_and_delta_kw(T,density/h2_ch4_mix_params['M_amu'],x_h2,x_ch4,h2_params['Tc'],ch4_params['Tc'],h2_params['rho_c']/h2_params['M_amu'],ch4_params['rho_c']/ch4_params['M_amu'],h2_ch4_mix_params['BetaT'],h2_ch4_mix_params['BetaV'],h2_ch4_mix_params['GammaT'],h2_ch4_mix_params['GammaV']) 
P=pressure_vector(density,T,tau,delta,'mix',h2_params,ch4_params,h2_ch4_mix_params,x_h2,x_ch4,rho_c,Tc)
h=enthalpy_vector(density,T,tau,delta,Tc,rho_c,'mix',h2_params,ch4_params,h2_ch4_mix_params,x_h2,x_ch4) 
Excess_h=Excess_Enthalpy_vector(density,T,tau,delta,Tc,rho_c,h2_params,ch4_params,h2_ch4_mix_params,x_h2,x_ch4)

B=Virial_B_vector(tau,h2_params,ch4_params,h2_ch4_mix_params,x_h2,x_ch4,rho_c,Tc,'mix')
C=Virial_C_vector(tau,h2_params,ch4_params,h2_ch4_mix_params,x_h2,x_ch4,rho_c,Tc,'mix')
print Excess_h #Excess Enthalpy in J/mol
print P-1000*P_refprop #Pressure in kPa
print h-M_amu*enthalpy_refprop #Enthalpy in J/mol M_amu used to convert from kJ/kg
print B
print C
print tau,delta
print P
#...and BEHOLD! IT WORKS!                                                        