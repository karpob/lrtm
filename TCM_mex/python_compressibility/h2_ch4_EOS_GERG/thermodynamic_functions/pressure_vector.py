def pressure_vector(rho,T,tau,delta,eqn_type,parameter_set1,parameter_set2,parameter_set12,x_1,x_2,rho_c12,Tc12):
        #Pressure is returned in kPa (kilo-Pascals)
        #Units will depend on critical density and Temperature (assumed to be in units of kg/m**3 and Kelvin respectively
        from helmholtz_functions.d_alpha_d_delta_vector import d_alpha_d_delta_vector
        from pylab import zeros
        dalpha_delta=zeros(len(tau))
        #rho=rho_c12*delta # get density back
        #T=Tc12/tau # get Temperature back
        Rspecific=parameter_set12['R']/parameter_set12['M_amu']
        
        
        if(eqn_type=='pure_substance'):
                dalpha_delta=d_alpha_d_delta_vector(tau,delta,parameter_set1)
        else:        
                dalpha_delta1=d_alpha_d_delta_vector(tau,delta,parameter_set1)
                dalpha_delta2=d_alpha_d_delta_vector(tau,delta,parameter_set2)
                dalpha_delta12=d_alpha_d_delta_vector(tau,delta,parameter_set12)
               
                dalpha_delta=x_1*dalpha_delta1+x_2*dalpha_delta2+x_1*x_2*dalpha_delta12
        
        Pressure=(rho)*Rspecific*T*(delta*dalpha_delta+float(1.0)) #units! rho-->kg/m**3
                                                                   #        R (ideal) --> m**3 Pa/mol
                                                                   #        M --> g/mol
                                                                   #        T--Kelvin
                                                                   #        (kg/m**3)*(m**3 Pa/mol)*(K)=kg*Pa/mol
                                                                   #        kg*Pa/mol * (1 mol/M (g)) * (1000g/kg)=1000Pa=1kPa                                                                         
        return Pressure 
