from gas_library import get_h2_params
from gas_library import get_ch4_params
from gas_library import get_he_params
from gas_library import get_h2o_params
from gas_library import get_h2_ch4_mix
from gas_library import get_h2_h2o_mix
from scale_tau_and_delta_kw_4comp import scale_tau_and_delta_kw_4comp
from thermodynamic_functions.pressure_vector4 import pressure_vector4
from numpy import asarray,zeros

h2_params=get_h2_params()
ch4_params=get_ch4_params()
he_params=get_he_params()
h2o_params=get_h2o_params()
h2_ch4_mix_params=get_h2_ch4_mix(0.5,0.5)
h2_h2o_mix_params=get_h2_h2o_mix(0.5,0.5)
density=asarray([0.1,0.1])
x_h2=asarray([0.99,0.99])
x_ch4=asarray([0.01,0.01])
x_h2o=asarray([0.00,0.00])
x_he=asarray([0.0,0.0])
T=asarray([500.0,500.0])


M_amu_mix=x_h2*h2_params['M_amu']+x_ch4*ch4_params['M_amu']+x_h2o*h2o_params['M_amu']+x_he*he_params['M_amu']
[tau,delta,Tc_mix,rho_c_mix]=scale_tau_and_delta_kw_4comp(T,density/M_amu_mix,x_h2,x_ch4,x_h2o,x_he,h2_params['Tc'],ch4_params['Tc'],h2o_params['Tc'],he_params['Tc'],h2_params['rho_c']/h2_params['M_amu'],ch4_params['rho_c']/ch4_params['M_amu'],h2o_params['rho_c']/h2o_params['M_amu'],he_params['rho_c']/he_params['M_amu'],h2_ch4_mix_params['BetaT'],h2_ch4_mix_params['BetaV'],h2_ch4_mix_params['GammaT'],h2_ch4_mix_params['GammaV'],h2_h2o_mix_params['BetaT'],h2_h2o_mix_params['BetaV'],h2_h2o_mix_params['GammaT'],h2_h2o_mix_params['GammaV'],1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)
print 'tau,delta',tau,delta
P=pressure_vector4(density,T,tau,delta,h2_params,ch4_params,h2o_params,he_params,h2_ch4_mix_params,h2_h2o_mix_params,x_h2,x_ch4,x_h2o,x_he,rho_c_mix,Tc_mix)
print P

