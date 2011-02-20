def residual_pressure(density,Pin,T,h2_params,h2o_params,mix_params,x_h2,x_h2o):
        from thermodynamic_functions.pressure_vector import pressure_vector
        from thermodynamic_functions.enthalpy import enthalpy
        from scale_tau_and_delta_kw import scale_tau_and_delta_kw
        from pylab import shape
        
        
        [tau,delta,Tc,rho_c]=scale_tau_and_delta_kw(T,density/mix_params['M_amu'],x_h2,x_h2o,h2_params['Tc'],h2o_params['Tc'],h2_params['rho_c']/h2_params['M_amu'],h2o_params['rho_c']/h2o_params['M_amu'],mix_params['BetaT'],mix_params['BetaV'],mix_params['GammaT'],mix_params['GammaV'])       
        
        Pcalc=pressure_vector(density,T,tau,delta,'mix',h2_params,h2o_params,mix_params,x_h2,x_h2o,rho_c,Tc)
  
        err=Pin-Pcalc
        
        return err
