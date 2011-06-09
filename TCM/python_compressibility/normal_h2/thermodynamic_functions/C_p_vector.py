def C_p_vector(delta,tau,component_params):
        from helmholtz_functions.d_alpha_d_delta_vector import d_alpha_d_delta_vector
        from helmholtz_functions.d_alpha_d_delta_d_tau_vector import d_alpha_d_delta_d_tau_vector
        from helmholtz_functions.d_alpha_d_tau_d_tau_vector import d_alpha_d_tau_d_tau_vector
        from helmholtz_functions.ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R_vector import ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R_vector
        from helmholtz_functions.d_alpha_d_delta_d_delta_vector import d_alpha_d_delta_d_delta_vector
        from helmholtz_functions.ideal_helmholtz_from_coef_dtau_dtau_vector import ideal_helmholtz_from_coef_dtau_dtau_vector
        from numpy import power
        
        R_specific=component_params['R']/component_params['M_amu']
        if(component_params['ideal_eqn_type']=='Cp'):
                
                dalpha_o_tau_tau=ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R_vector(component_params['ni'],component_params['ti'],component_params['vi'],component_params['ui'],component_params['R'],component_params['ho'],component_params['so'],component_params['rho_c']/component_params['M_amu'],component_params['Tc'],tau,delta,component_params['n_ideal_gas_terms_pow'],component_params['n_ideal_gas_terms_exp'],component_params['To'],component_params['rho_o'])
        else:
                dalpha_o_tau_tau=ideal_helmholtz_from_coef_dtau_dtau_vector(delta,tau,component_params['ideal_n'],component_params['ideal_gamma'])
        
        dalpha_tau_tau=d_alpha_d_tau_d_tau_vector(tau,delta,component_params)
        dalpha_delta=d_alpha_d_delta_vector(tau,delta,component_params)
        dalpha_delta_tau=d_alpha_d_delta_d_tau_vector(tau,delta,component_params)      
        dalpha_delta_delta=d_alpha_d_delta_d_delta_vector(tau,delta,component_params)        
        cv=-1.0*R_specific*(tau*tau*(dalpha_o_tau_tau+dalpha_tau_tau))
        cp=cv+R_specific*power((1.0+delta*dalpha_delta-delta*tau*dalpha_delta_tau),2)/(1.0+2.0*delta*dalpha_delta+delta*delta*dalpha_delta_delta)
        
        #default is in kJ/kg-K convert to erg/K/mol for DeBoer's TCM
	#  kJ    1000J    1e7 ergs    1 kg    M_amu g      erg
	#----- * ----- * --------- * ------ * --------  = ------- 
	# kg*K    kJ      1J         1000 g    mol         K mol
	
	#cp*M_amu*1e7
        return cp
