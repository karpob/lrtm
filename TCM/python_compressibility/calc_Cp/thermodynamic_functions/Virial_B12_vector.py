def Virial_B12_vector(tau,parameter_set1,parameter_set2,parameter_set12,x_1,x_2,rho_c12,Tc12):
        # this computes the fugacity for a mixture using finitite difference approximation
        # using "central difference"
        rho_in_limit==0.00000001 # taken as zero limit in REFPROP
        from helmholtz_functions.d_alpha_d_delta_vector import d_alpha_d_delta_vector
        from numpy import zeros
        dalpha_delta=zeros(len(tau))
        #rho=rho_c12*delta # get density back
        T=Tc12/tau # get Temperature back
        Rspecific=parameter_set12['R']/parameter_set12['M_amu']
        tau=Tc12/T
        delta=rho_in_limit/rho_c12
        dalpha_d_delta1=d_alpha_d_delta_vector(tau,delta,parameter_set1)
        dalpha_d_delta2=d_alpha_d_delta_vector(tau,delta,parameter_set2)                    
        dalpha_d_delta12=d_alpha_d_delta_vector(tau,delta,parameter_set12)
        dalpha_d_delta=dalpha_d_delta1+dalpha_d_delta2+dalpha_d_delta12
        Bmix=dalpha_d_delta/rho_in_limit
        
        tau1=parameter_set1['Tc']/T
        delta1=rho_in_limit/parameter_set1['rho_c']
        dalpha_1=d_alpha_d_delta_vector(tau1,delta1,parameter_set1)
        B1=dalpha_1/rho_in_limit
        
        tau2=parameter_set2['Tc']/T
        delta2=rho_in_limit/parameter_set2['rho_c']
        dalpha_2=d_alpha_d_delta_vector(tau1,delta1,parameter_set1)
        B2=dalpha_2/rho_in_limit
        
        B12=(Bmix-x_1*x_1*B1-x_2*x_2*B2)/2.0/x_1/x_2 
        return B12
