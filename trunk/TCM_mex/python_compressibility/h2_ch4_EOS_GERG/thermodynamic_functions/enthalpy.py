def enthalpy(tau,delta,ni,ti,vi,ui,R,ho,so,rho_c,Tc,n_ideal_gas_terms_pow,n_ideal_gas_terms_exp,To,rho_o,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A,M_amu,ideal_n,ideal_gamma,ideal_eqn_type):       
        from helmholtz_functions.ideal_helmholtz_energy_dtau_from_Cp_o_over_R import ideal_helmholtz_energy_dtau_from_Cp_o_over_R
        from helmholtz_functions.d_alpha_d_tau import d_alpha_d_tau
        from helmholtz_functions.d_alpha_d_delta import d_alpha_d_delta
        from helmholtz_functions.ideal_helmholtz_from_coef_dtau import ideal_helmholtz_from_coef_dtau
        from pylab import zeros
        dalpha_tau=zeros(len(tau))
        dalpha_delta=zeros(len(tau))
        dalpha_o=zeros(len(tau))
        #M_H2=2.01594
        T=Tc/tau
        if(ideal_eqn_type=='Cp'):
                for i in range(0,len(tau)):
                        dalpha_o[i]=ideal_helmholtz_energy_dtau_from_Cp_o_over_R(ni,ti,vi,ui,R,ho,so,rho_c/M_amu,Tc,tau[i],delta[i],n_ideal_gas_terms_pow,n_ideal_gas_terms_exp,To,rho_o)
        if(ideal_eqn_type=='Coef'):
                for i in range(0,len(tau)):
                        dalpha_o[i]=ideal_helmholtz_from_coef_dtau(delta[i],tau[i],ideal_n,ideal_gamma)
        if(ideal_eqn_type=='none'):
                dalpha_o[1]==0.0        
        
        for i in range(0,len(tau)):                
                dalpha_tau[i]=d_alpha_d_tau(tau[i],delta[i],N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
                dalpha_delta[i]=d_alpha_d_delta(tau[i],delta[i],N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
        R_specific=R/M_amu
        h1=R_specific*T*1.0
	h2=R_specific*T*tau*(dalpha_o+dalpha_tau)
	h3=R_specific*T*delta*dalpha_delta
	h=h1+h2+h3
	return h
	
