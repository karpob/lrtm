def residual_virial_C(density,Cin,T,h2_params,h2o_params,mix_params,x_h2,x_h2o):
        from thermodynamic_functions.Virial_C_vector import Virial_C_vector
        from scale_tau_and_delta_kw import scale_tau_and_delta_kw
        
        [tau,delta,Tc,rho_c]=scale_tau_and_delta_kw(T,density/mix_params['M_amu'],x_h2,x_h2o,h2_params['Tc'],h2o_params['Tc'],h2_params['rho_c']/h2_params['M_amu'],h2o_params['rho_c']/h2o_params['M_amu'],mix_params['BetaT'],mix_params['BetaV'],mix_params['GammaT'],mix_params['GammaV'])                           
        
        Ccalc=Virial_C_vector(tau,h2_params,h2o_params,mix_params,x_h2,x_h2o,rho_c,Tc,'mix')
       
        
        
        err=Cin-Ccalc
        
        return err
