def residual_pressure_pure_substance(density,Pin,T,params):
        from thermodynamic_functions.pressure_vector import pressure_vector
        from pylab import shape
        
        tau=params['Tc']/T
        delta=density/params['rho_c']
        Pcalc=pressure_vector(density,T,tau,delta,'pure_substance',params,0,params,0,0,params['rho_c'],params['Tc'])
       
        
        
        err=Pin-Pcalc
        
        return err
