def residual_enthalpy(p,actual_enthalpy,P,T,x_h2,x_h2o):
        from pylab import shape
        from enthalpy_h2_h2o_mixture import enthalpy_h2_h2o_mixture        
        F_ij=p[0]
        beta_ij=p[1]
        phi_ij=p[2]
        sigma_ij=p[3]
        xi_ij=p[4]
        calculated_enthalpy=enthalpy_h2_h2o_mixture(T,P,x_h2,x_h2o,F_ij,beta_ij,phi_ij,sigma_ij,xi_ij)
        err=actual_enthalpy-calculated_enthalpy
        print err
        return err
