def Excess_Enthalpy_vector(rho,T,tau,delta,Tc,rho_c,parameter_set1,parameter_set2,parameter_set12,x_1,x_2):
        from thermodynamic_functions.enthalpy_vector import enthalpy_vector
        from thermodynamic_functions.pressure_vector import pressure_vector
        from residual_pressure_pure_substance import residual_pressure_pure_substance
        from scipy.optimize import leastsq
        
        Enthalpy_of_mixture=enthalpy_vector(rho,T,tau,delta,Tc,rho_c,'mix',parameter_set1,parameter_set2,parameter_set12,x_1,x_2)
        Pressure_of_mixture=pressure_vector(rho,T,tau,delta,'mix',parameter_set1,parameter_set2,parameter_set12,x_1,x_2,Tc,rho_c)
        
        density_guess_component_1=Pressure_of_mixture/((parameter_set1['R']/parameter_set1['M_amu'])*T)
        density_guess_component_2=Pressure_of_mixture/((parameter_set2['R']/parameter_set2['M_amu'])*T)
        
        [density_component_1,message_1]=leastsq(residual_pressure_pure_substance,density_guess_component_1,args=(Pressure_of_mixture,T,parameter_set1))
        [density_component_2,message_2]=leastsq(residual_pressure_pure_substance,density_guess_component_2,args=(Pressure_of_mixture,T,parameter_set2))
        
        Tc1=parameter_set1['Tc']
        Tc2=parameter_set2['Tc']
        
        rho_c1=parameter_set1['rho_c']
        rho_c2=parameter_set2['rho_c']
        
        tau_1=Tc1/T
        tau_2=Tc2/T
        
        delta_1=density_component_1/parameter_set1['rho_c']
        delta_2=density_component_2/parameter_set2['rho_c']
         
        Enthalpy_component_1=enthalpy_vector(density_component_1,T,tau_1,delta_1,Tc1,rho_c1,'pure_substance',parameter_set1,parameter_set1,parameter_set1,1,1)
        Enthalpy_component_2=enthalpy_vector(density_component_2,T,tau_2,delta_2,Tc2,rho_c2,'pure_substance',parameter_set2,parameter_set2,parameter_set2,1,1)
        
        Excess_Enthalpy=Enthalpy_of_mixture-x_1*Enthalpy_component_1-x_2*Enthalpy_component_2
       
        return Excess_Enthalpy
