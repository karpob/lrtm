def scale_tau_and_delta_kw(T,rho,x1,x2,T_c1,T_c2,rho_c1,rho_c2,BetaT,BetaV,GammaT,GammaV):
        from numpy import power
        #parameters for fitting xi_ij
        #                       sigma_ij
        #                       beta_ij
        #                       phi_ij
        
         
        front_matter_T_scaled=2.0*x1*x2*BetaT*GammaT
        front_matter_rho_scaled=2.0*x1*x2*BetaV*GammaV
        
        FracT=(x1+x2)/(BetaT*BetaT*x1+BetaT*BetaT*x2)
        FracRho=(x1+x2)/(BetaV*BetaV*x1+BetaV*BetaV*x2)
        
        T_c_scaled=x1*x1*T_c1+x2*x2*T_c2+front_matter_T_scaled*FracT*power((T_c1*T_c2),0.5)        
        one_third_pow_1=power(rho_c1,-1.0/3.0)
        one_third_pow_2=power(rho_c2,-1.0/3.0)
        inv_rho_c=x1*x1*(1.0/rho_c1)+x2*x2*(1.0/rho_c2)+front_matter_rho_scaled*FracRho*power((one_third_pow_1+one_third_pow_2),3.0)/8.0
        
        rho_c_scaled=1.0/inv_rho_c
        T_c_scaled=T_c_scaled
        tau=T_c_scaled/T
        delta=rho/(rho_c_scaled)
        
        T_c_scaled[x1==0.0]=0.0
        rho_c_scaled[x1==0.0]=0.0
        T_c_scaled[x2==0.0]=0.0
        rho_c_scaled[x2==0.0]=0.0
        
        T_c_scaled[x1==1.0]=T_c1
        T_c_scaled[x2==1.0]=T_c2
        
        rho_c_scaled[x1==1.0]=rho_c1
        rho_c_scaled[x2==1.0]=rho_c2
        return tau,delta,T_c_scaled,rho_c_scaled        
