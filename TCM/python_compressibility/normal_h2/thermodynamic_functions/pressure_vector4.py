def pressure_vector4(rho,T,tau,delta,parameter_set1,parameter_set2,parameter_set3,parameter_set4,parameter_set12,parameter_set13,x_1,x_2,x_3,x_4,rho_c_mix,Tc_mix):
        #Pressure is returned in kPa (kilo-Pascals)
        #Units will depend on critical density and Temperature (assumed to be in units of kg/m**3 and Kelvin respectively
        from helmholtz_functions.d_alpha_d_delta_vector import d_alpha_d_delta_vector
        from numpy import zeros
        dalpha_delta=zeros(len(tau))
        
        Rnum=x_1*parameter_set1['R']+x_2*parameter_set2['R']+x_3*parameter_set3['R']+x_4*parameter_set4['R']
        Mden=x_1*parameter_set1['M_amu']+x_2*parameter_set2['M_amu']+x_3*parameter_set3['M_amu']+x_4*parameter_set4['M_amu']
        Rspecific=Rnum/Mden
        dalpha_delta1=d_alpha_d_delta_vector(tau,delta,parameter_set1)
        dalpha_delta2=d_alpha_d_delta_vector(tau,delta,parameter_set2)
        dalpha_delta3=d_alpha_d_delta_vector(tau,delta,parameter_set3)
        dalpha_delta4=d_alpha_d_delta_vector(tau,delta,parameter_set4)
        
        dalpha_delta12=d_alpha_d_delta_vector(tau,delta,parameter_set12)
        dalpha_delta13=d_alpha_d_delta_vector(tau,delta,parameter_set13)
        
        dalpha_delta=x_1*dalpha_delta1+x_2*dalpha_delta2+x_3*dalpha_delta3+x_4*dalpha_delta4+x_1*x_2*dalpha_delta12+x_1*x_3*dalpha_delta13
       
        Pressure=(rho)*Rspecific*T*(delta*dalpha_delta+float(1.0)) #units! rho-->kg/m**3
                                                                   #        R (ideal) --> m**3 Pa/mol
                                                                   #        M --> g/mol
                                                                   #        T--Kelvin
                                                                   #        (kg/m**3)*(m**3 Pa/mol)*(K)=kg*Pa/mol
                                                                   #        kg*Pa/mol * (1 mol/M (g)) * (1000g/kg)=1000Pa=1kPa 
                                                                                                                                        
        return Pressure 
