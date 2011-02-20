from gas_library import get_h2_params
from gas_library import get_ch4_params
from gas_library import get_he_params
from gas_library import get_h2o_params
from gas_library import get_h2_ch4_mix
from gas_library import get_h2_h2o_mix
from thermodynamic_functions.C_p_vector_mix import C_p_vector_mix
from scale_tau_and_delta_kw_4comp import scale_tau_and_delta_kw_4comp
from thermodynamic_functions.pressure_vector4 import pressure_vector4
from thermodynamic_functions.C_p_vector_mix4 import C_p_vector_mix4

from scipy.optimize import fsolve
from numpy import asarray,zeros

f=open('python_compressibility/calc_Cp/output_T_P.txt','r')
data=f.readline()
f.close()
T_i,P_h2_i,P_he_i,P_ch4_i,P_h2o_i=data.split("\t")

parameters_H2=get_h2_params()
parameters_h2o=get_h2o_params()
parameters_ch4=get_ch4_params()
parameters_he=get_he_params()


Bars_to_kPa=100
T=asarray([float(T_i),float(T_i)])
P_h2=asarray([float(P_h2_i),float(P_h2_i)])*Bars_to_kPa
P_he=asarray([float(P_he_i),float(P_he_i)])*Bars_to_kPa
P_ch4=asarray([float(P_ch4_i),float(P_ch4_i)])*Bars_to_kPa
P_h2o=asarray([float(P_h2o_i),float(P_h2o_i)])*Bars_to_kPa

                                             

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

density_ch4=P_ch4/((parameters_ch4['R']/parameters_ch4['M_amu'])*T)
#density_ch4=fsolve(residual_pressure_pure_substance,density_guess_ch4,args=(P_ch4,T,parameters_ch4))
delta_ch4=density_ch4/parameters_ch4['rho_c']


density=density_h2+density_h2o+density_ch4+density_he
#note we're calculating x_h2 and x_ch4 according to the mixture of just H2+CH4, not the whole thing He+h2o+etc...
molar_density_h2=density_h2/parameters_H2['M_amu']
molar_density_he=density_he/parameters_he['M_amu']
molar_density_ch4=density_ch4/parameters_ch4['M_amu']
molar_density_h2o=density_h2o/parameters_h2o['M_amu']

x_h2=molar_density_h2/(molar_density_h2+molar_density_ch4+molar_density_he+molar_density_h2o)
x_ch4=molar_density_ch4/(molar_density_h2+molar_density_ch4+molar_density_he+molar_density_h2o)
x_he=molar_density_he/(molar_density_h2+molar_density_ch4+molar_density_he+molar_density_h2o)
x_h2o=molar_density_h2o/(molar_density_h2+molar_density_ch4+molar_density_he+molar_density_h2o)

h2_ch4_mix_params=get_h2_ch4_mix(x_h2,x_ch4)
h2_h2o_mix_params=get_h2_h2o_mix(x_h2,x_h2o)

M_amu_mix=x_h2*parameters_H2['M_amu']+x_he*parameters_he['M_amu']+x_ch4*parameters_ch4['M_amu']+x_h2o*parameters_h2o['M_amu']
[tau,delta,Tc_mix,rho_c_mix]=scale_tau_and_delta_kw_4comp(T,density/M_amu_mix,x_h2,x_ch4,x_h2o,x_he,parameters_H2['Tc'],parameters_ch4['Tc'],parameters_h2o['Tc'],parameters_he['Tc'],parameters_H2['rho_c']/parameters_H2['M_amu'],parameters_ch4['rho_c']/parameters_ch4['M_amu'],parameters_h2o['rho_c']/parameters_h2o['M_amu'],parameters_he['rho_c']/parameters_he['M_amu'],h2_ch4_mix_params['BetaT'],h2_ch4_mix_params['BetaV'],h2_ch4_mix_params['GammaT'],h2_ch4_mix_params['GammaV'],h2_h2o_mix_params['BetaT'],h2_h2o_mix_params['BetaV'],h2_h2o_mix_params['GammaT'],h2_h2o_mix_params['GammaV'],1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0)
        
Pcalc=pressure_vector4(density,T,tau,delta,parameters_H2,parameters_ch4,parameters_h2o,parameters_he,h2_ch4_mix_params,h2_h2o_mix_params,x_h2,x_ch4,x_h2o,x_he,rho_c_mix,Tc_mix)

CP_mix=C_p_vector_mix4(density,T,delta,tau,parameters_H2,parameters_ch4,parameters_he,parameters_h2o,h2_ch4_mix_params,h2_h2o_mix_params,x_h2,x_ch4,x_he,x_h2o)
Pcalc=Pcalc/Bars_to_kPa
try:
   f=open('python_compressibility/calc_Cp/input_Cp.txt','w')
   f.write(str(CP_mix[0])+' '+str(Pcalc[0]))
   f.close()
   print 'success',str(Pcalc[0])
except:
   f=open('python_compressibility/calc_Cp/fail.txt','w')
   f.write('failed to write input_Cp.txt');
   f.close()
   print 'fail'
