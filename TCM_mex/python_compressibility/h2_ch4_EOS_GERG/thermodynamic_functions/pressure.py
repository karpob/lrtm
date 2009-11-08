def pressure(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms_wo_exp,n_power_terms_w_exp,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A,R,Tc,rho_c,M_amu):
        #Pressure is returned in kPa (kilo-Pascals)
        #Units will depend on critical density and Temperature (assumed to be in units of kg/m**3 and Kelvin respectively
        from helmholtz_functions.d_alpha_d_delta_vector import d_alpha_d_delta_vector
        from pylab import zeros
        dalpha_delta=zeros(len(tau))
        rho=rho_c*delta # get density back
        T=Tc/tau # get Temperature back
        Rspecific=R/M_amu
        
        dalpha_delta=d_alpha_d_delta_vector(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms_wo_exp,n_power_terms_w_exp,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
        Pressure=(rho)*Rspecific*T*(delta*dalpha_delta+float(1.0)) #units! rho-->kg/m**3
                                                                   #        R (ideal) --> m**3 Pa/mol
                                                                   #        M --> g/mol
                                                                   #        T--Kelvin
                                                                   #        (kg/m**3)*(m**3 Pa/mol)*(K)=kg*Pa/mol
                                                                   #        kg*Pa/mol * (1 mol/M (g)) * (1000g/kg)=1000Pa=1kPa                                                                         
        return Pressure 
