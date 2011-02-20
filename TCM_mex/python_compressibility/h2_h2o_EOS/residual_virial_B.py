def residual_virial_B(density,Bin,T,h2_params,h2o_params,h2o_h2_mix_params,x_h2,x_h2o):
        from thermodynamic_functions.Virial_B_vector import Virial_B_vector
        from scale_tau_and_delta_kw import scale_tau_and_delta_kw
        
        [tau,delta,Tc,rho_c]=scale_tau_and_delta_kw(T,density/mix_params['M_amu'],x_h2,x_h2o,h2_params['Tc'],h2o_params['Tc'],h2_params['rho_c']/h2_params['M_amu'],h2o_params['rho_c']/h2o_params['M_amu'],mix_params['BetaT'],mix_params['BetaV'],mix_params['GammaT'],mix_params['GammaV'])                   
        #Pcalc=pressure_vector(tau,delta,'mix',h2_params,h2o_params,h2o_h2_mix_params,x_h2,x_h2o,rho_c,Tc)
        #print tau#,delta,rho_c,Tc
        Bcalc=Virial_B_vector(tau,h2_params,h2o_params,h2o_h2_mix_params,x_h2,x_h2o,rho_c,Tc,'mix')
       
        
        
        err=Bin-Bcalc
        
        return err
