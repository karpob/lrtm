def residual_pressure_pure_substance(density,Pin,T,params):
        from thermodynamic_functions.pressure_vector import pressure_vector
        from numpy import shape
        
        tau=params['Tc']/T
        delta=density/params['rho_c']
        
        Pcalc=pressure_vector(density,T,tau,delta,'pure_substance',params,params,params,0,0,0,0)
       
        
        
        err=Pin-Pcalc
        
        return err
