def scale_tau_and_delta(T,rho,x1,x2,T_c1,T_c2,rho_c1,rho_c2,beta_ij,phi_ij,sigma_ij,xi_ij):
        #parameters for fitting xi_ij
        #                       sigma_ij
        #                       beta_ij
        #                       phi_ij
        
        #T_c_scaled=0.0
        #rho_c_scaled=0.0
        #for i in range(0,len(T_c)):
        #        T_c_scaled=x[i]*T_c[i]+T_c_scaled
        #        rho_c_scaled=x[i]/rho_c+rho_c_scaled
        T_c_scaled=x1*T_c1+x2*T_c2
        T_c_scaled=T_c_scaled+x1**(beta_ij)*x2**(phi_ij)*sigma_ij        
         
        #for i in range(0,len(T_c)-1):
        #        for j in range(1,len(T_c)):
        #                T_c_scaled=x[i]**(beta)*x[j]**(phi)*sigma_ij+T_c_scaled
        #                rho_c_scaled=x[i]*x[j]*xi_ij+rho_scaled
        
        rho_c_scaled=x1/rho_c1+x2/rho_c2+x1*x2*xi_ij
        
        rho_c_scaled=rho_c_scaled**(-1.0)
        
        tau=T_c_scaled/T
        delta=rho/rho_c_scaled
        return tau,delta,T_c_scaled,rho_c_scaled        
