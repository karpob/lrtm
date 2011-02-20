def enthalpy_vector(rho,T,tau,delta,Tc,rho_c,eqn_type,parameter_set1,parameter_set2,parameter_set12,x_1,x_2):       
        from helmholtz_functions.ideal_helmholtz_energy_dtau_from_Cp_o_over_R_vector import ideal_helmholtz_energy_dtau_from_Cp_o_over_R_vector
        from helmholtz_functions.d_alpha_d_tau_vector import d_alpha_d_tau_vector
        from helmholtz_functions.d_alpha_d_delta_vector import d_alpha_d_delta_vector
        from helmholtz_functions.ideal_helmholtz_from_coef_dtau_vector import ideal_helmholtz_from_coef_dtau_vector
        from numpy import zeros,shape
        
        dalpha_tau=zeros(len(tau))
        dalpha_delta=zeros(len(tau))
        dalpha_o=zeros(len(tau))
  
        
        if(eqn_type=='mix'):
                rho=rho/parameter_set12['M_amu'] #transform density to mol/L from kg/m^3
                if(parameter_set1['ideal_eqn_type']=='Cp'):
                        dalpha_o1=ideal_helmholtz_energy_dtau_from_Cp_o_over_R_vector(parameter_set1['ni'],parameter_set1['ti'],parameter_set1['vi'],parameter_set1['ui'],parameter_set1['R'],parameter_set1['ho'],parameter_set1['so'],parameter_set1['rho_c']/parameter_set1['M_amu'],parameter_set1['Tc'],parameter_set1['Tc']/T,rho/(parameter_set1['rho_c']/parameter_set1['M_amu']),parameter_set1['n_ideal_gas_terms_pow'],parameter_set1['n_ideal_gas_terms_exp'],parameter_set1['To'],parameter_set1['rho_o'])
                if(parameter_set1['ideal_eqn_type']=='Coef'):
                        dalpha_o1=ideal_helmholtz_from_coef_dtau_vector(rho/(parameter_set1['rho_c']/parameter_set1['M_amu']),parameter_set1['Tc']/T,parameter_set1['ideal_n'],parameter_set1['ideal_gamma'])
                if(parameter_set2['ideal_eqn_type']=='Cp'):
                        dalpha_o2=ideal_helmholtz_energy_dtau_from_Cp_o_over_R_vector(parameter_set2['ni'],parameter_set2['ti'],parameter_set2['vi'],parameter_set2['ui'],parameter_set2['R'],parameter_set2['ho'],parameter_set2['so'],parameter_set2['rho_c']/parameter_set2['M_amu'],parameter_set2['Tc'],parameter_set2['Tc']/T,rho/(parameter_set2['rho_c']/parameter_set2['M_amu']),parameter_set2['n_ideal_gas_terms_pow'],parameter_set2['n_ideal_gas_terms_exp'],parameter_set2['To'],parameter_set2['rho_o'])
                if(parameter_set2['ideal_eqn_type']=='Coef'):
                        dalpha_o2=ideal_helmholtz_from_coef_dtau_vector(rho/(parameter_set2['rho_c']/parameter_set2['M_amu']),parameter_set2['Tc']/T,parameter_set2['ideal_n'],parameter_set2['ideal_gamma'])
                
                dalpha_tau1=d_alpha_d_tau_vector(tau,delta,parameter_set1)
                dalpha_tau2=d_alpha_d_tau_vector(tau,delta,parameter_set2)
                dalpha_tau12=d_alpha_d_tau_vector(tau,delta,parameter_set12)
                
                dalpha_delta1=d_alpha_d_delta_vector(tau,delta,parameter_set1)
                dalpha_delta2=d_alpha_d_delta_vector(tau,delta,parameter_set2)
                dalpha_delta12=d_alpha_d_delta_vector(tau,delta,parameter_set12)
                
                R_tau_dalpha_o=parameter_set1['R']*x_1*(dalpha_o1)*parameter_set1['Tc']/T+parameter_set2['R']*x_2*(dalpha_o2)*parameter_set2['Tc']/T
                
                dalpha_tau=x_1*dalpha_tau1+x_2*dalpha_tau2+x_1*x_2*dalpha_tau12
                
                dalpha_delta=x_1*dalpha_delta1+x_2*dalpha_delta2+x_1*x_2*dalpha_delta12
                
                h2=tau*dalpha_tau
                h3=delta*dalpha_delta
                
                h=parameter_set12['R']*T*(1.0+h2+h3)+R_tau_dalpha_o*T      
                 
        if(eqn_type=='pure_substance'):
                if(parameter_set1['ideal_eqn_type']=='Cp'):
                        dalpha_o=ideal_helmholtz_energy_dtau_from_Cp_o_over_R_vector(parameter_set1['ni'],parameter_set1['ti'],parameter_set1['vi'],parameter_set1['ui'],parameter_set1['R'],parameter_set1['ho'],parameter_set1['so'],parameter_set1['rho_c']/parameter_set1['M_amu'],parameter_set1['Tc'],parameter_set1['Tc']/T,rho/parameter_set1['rho_c'],parameter_set1['n_ideal_gas_terms_pow'],parameter_set1['n_ideal_gas_terms_exp'],parameter_set1['To'],parameter_set1['rho_o'])
                if(parameter_set1['ideal_eqn_type']=='Coef'):
                        dalpha_o=ideal_helmholtz_from_coef_dtau_vector(delta,tau,parameter_set1['ideal_n'],parameter_set1['ideal_gamma'])
                dalpha_tau=d_alpha_d_tau_vector(tau,delta,parameter_set1)
                dalpha_delta=d_alpha_d_delta_vector(tau,delta,parameter_set1)
                R=parameter_set1['R']#/parameter_set1['M_amu']
                h1=1.0
	        h2=tau*(dalpha_o+dalpha_tau)
	        h3=delta*dalpha_delta
	        h=R*T*(h1+h2+h3)   #h in J/mol?
	        
	        
	return h
	
