def Virial_B_vector(tau,parameter_set1,parameter_set2,parameter_set12,x_1,x_2,rho_c12,Tc12,mix_or_pure):
        
        rho_in_limit=1e-8 # taken as zero limit in REFPROP
        from helmholtz_functions.d_alpha_d_delta_vector import d_alpha_d_delta_vector
        from numpy import zeros,shape,ones
        
        delta=rho_in_limit*ones(shape(tau))/rho_c12
        
        if (mix_or_pure=='mix'):
                dalpha_d_delta1=d_alpha_d_delta_vector(tau,delta,parameter_set1)
                dalpha_d_delta2=d_alpha_d_delta_vector(tau,delta,parameter_set2)                    
                dalpha_d_delta12=d_alpha_d_delta_vector(tau,delta,parameter_set12)
                dalpha_d_delta=x_1*dalpha_d_delta1+x_2*dalpha_d_delta2+x_1*x_2*dalpha_d_delta12
        if(mix_or_pure=='pure'):
                dalpha_d_delta=d_alpha_d_delta_vector(tau,delta,parameter_set1)        
        
        B=dalpha_d_delta/rho_c12  # returns Second Virial coefficient which is consistent with 
                                  # REFERENCED NOTATION, in place of REFPROP KLUUUUDGE!
        
        return B 
