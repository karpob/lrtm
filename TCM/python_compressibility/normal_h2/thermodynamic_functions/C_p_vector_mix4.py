def C_p_vector_mix4(density,T,delta,tau,component_params1,component_params2,component_params3,component_params4,component_params12,component_params14,x_1,x_2,x_3,x_4):
        from helmholtz_functions.d_alpha_d_delta_vector import d_alpha_d_delta_vector
        from helmholtz_functions.d_alpha_d_delta_d_tau_vector import d_alpha_d_delta_d_tau_vector
        from helmholtz_functions.d_alpha_d_tau_d_tau_vector import d_alpha_d_tau_d_tau_vector
        from helmholtz_functions.ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R_vector import ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R_vector
        from helmholtz_functions.d_alpha_d_delta_d_delta_vector import d_alpha_d_delta_d_delta_vector
        from helmholtz_functions.ideal_helmholtz_from_coef_dtau_dtau_vector import ideal_helmholtz_from_coef_dtau_dtau_vector
        from numpy import power
        
        
        if(component_params1['ideal_eqn_type']=='Cp'):
         
                dalpha_o_tau_tau1=ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R_vector(component_params1['ni'],component_params1['ti'],component_params1['vi'],component_params1['ui'],component_params1['R'],component_params1['ho'],component_params1['so'],component_params1['rho_c']/component_params1['M_amu'],component_params1['Tc'],component_params1['Tc']/T,density/component_params1['rho_c'],component_params1['n_ideal_gas_terms_pow'],component_params1['n_ideal_gas_terms_exp'],component_params1['To'],component_params1['rho_o'])
        else:
               
                dalpha_o_tau_tau1=ideal_helmholtz_from_coef_dtau_dtau_vector(density/component_params1['rho_c'],component_params1['Tc']/T,component_params1['ideal_n'],component_params1['ideal_gamma'])
        
        if(component_params2['ideal_eqn_type']=='Cp'):
               
                
                dalpha_o_tau_tau2=ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R_vector(component_params2['ni'],component_params2['ti'],component_params2['vi'],component_params2['ui'],component_params2['R'],component_params2['ho'],component_params2['so'],component_params2['rho_c'],component_params2['Tc'],component_params2['Tc']/T,density/component_params2['rho_c'],component_params2['n_ideal_gas_terms_pow'],component_params2['n_ideal_gas_terms_exp'],component_params2['To'],component_params2['rho_o'])
        else:
               
                dalpha_o_tau_tau2=ideal_helmholtz_from_coef_dtau_dtau_vector(density/component_params2['rho_c'],component_params2['Tc']/T,component_params2['ideal_n'],component_params2['ideal_gamma'])
        
        if(component_params3['ideal_eqn_type']=='Cp'):
               
                
                dalpha_o_tau_tau3=ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R_vector(component_params3['ni'],component_params3['ti'],component_params3['vi'],component_params3['ui'],component_params3['R'],component_params3['ho'],component_params3['so'],component_params3['rho_c'],component_params3['Tc'],component_params3['Tc']/T,density/component_params3['rho_c'],component_params3['n_ideal_gas_terms_pow'],component_params3['n_ideal_gas_terms_exp'],component_params3['To'],component_params3['rho_o'])
        else:
               
                dalpha_o_tau_tau3=ideal_helmholtz_from_coef_dtau_dtau_vector(density/component_params3['rho_c'],component_params3['Tc']/T,component_params3['ideal_n'],component_params3['ideal_gamma'])
        
        if(component_params4['ideal_eqn_type']=='Cp'):
               
                
                dalpha_o_tau_tau4=ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R_vector(component_params4['ni'],component_params4['ti'],component_params4['vi'],component_params4['ui'],component_params4['R'],component_params4['ho'],component_params4['so'],component_params4['rho_c'],component_params4['Tc'],component_params4['Tc']/T,density/component_params4['rho_c'],component_params4['n_ideal_gas_terms_pow'],component_params4['n_ideal_gas_terms_exp'],component_params4['To'],component_params4['rho_o'])
        else:
               
                dalpha_o_tau_tau4=ideal_helmholtz_from_coef_dtau_dtau_vector(density/component_params4['rho_c'],component_params4['Tc']/T,component_params4['ideal_n'],component_params4['ideal_gamma'])
                
        
        dalpha_tau_tau1=d_alpha_d_tau_d_tau_vector(tau,delta,component_params1)
        dalpha_delta1=d_alpha_d_delta_vector(tau,delta,component_params1)
        dalpha_delta_tau1=d_alpha_d_delta_d_tau_vector(tau,delta,component_params1)      
        dalpha_delta_delta1=d_alpha_d_delta_d_delta_vector(tau,delta,component_params1)
        
        dalpha_tau_tau2=d_alpha_d_tau_d_tau_vector(tau,delta,component_params2)
        dalpha_delta2=d_alpha_d_delta_vector(tau,delta,component_params2)
        dalpha_delta_tau2=d_alpha_d_delta_d_tau_vector(tau,delta,component_params2)      
        dalpha_delta_delta2=d_alpha_d_delta_d_delta_vector(tau,delta,component_params2)
        
        dalpha_tau_tau3=d_alpha_d_tau_d_tau_vector(tau,delta,component_params3)
        dalpha_delta3=d_alpha_d_delta_vector(tau,delta,component_params3)
        dalpha_delta_tau3=d_alpha_d_delta_d_tau_vector(tau,delta,component_params3)      
        dalpha_delta_delta3=d_alpha_d_delta_d_delta_vector(tau,delta,component_params3)
        
        dalpha_tau_tau4=d_alpha_d_tau_d_tau_vector(tau,delta,component_params4)
        dalpha_delta4=d_alpha_d_delta_vector(tau,delta,component_params4)
        dalpha_delta_tau4=d_alpha_d_delta_d_tau_vector(tau,delta,component_params4)      
        dalpha_delta_delta4=d_alpha_d_delta_d_delta_vector(tau,delta,component_params4)
        
        
        dalpha_tau_tau12=d_alpha_d_tau_d_tau_vector(tau,delta,component_params12)
        dalpha_delta12=d_alpha_d_delta_vector(tau,delta,component_params12)
        dalpha_delta_tau12=d_alpha_d_delta_d_tau_vector(tau,delta,component_params12)      
        dalpha_delta_delta12=d_alpha_d_delta_d_delta_vector(tau,delta,component_params12)
        
        dalpha_tau_tau14=d_alpha_d_tau_d_tau_vector(tau,delta,component_params14)
        dalpha_delta14=d_alpha_d_delta_vector(tau,delta,component_params14)
        dalpha_delta_tau14=d_alpha_d_delta_d_tau_vector(tau,delta,component_params14)      
        dalpha_delta_delta14=d_alpha_d_delta_d_delta_vector(tau,delta,component_params14)
        
        
        
        dalpha_tau_tau=x_1*dalpha_tau_tau1+x_2*dalpha_tau_tau2 +x_3*dalpha_tau_tau3 +x_4*dalpha_tau_tau4 +x_1*x_2*dalpha_tau_tau12 + x_1*x_4*dalpha_tau_tau14
        dalpha_delta=x_1*dalpha_delta1+x_2*dalpha_delta2+x_3*dalpha_delta3+x_4*dalpha_delta4+x_1*x_2*dalpha_delta12+x_1*x_4*dalpha_delta14
        dalpha_delta_tau=x_1*dalpha_delta_tau1+x_2*dalpha_delta_tau2+x_3*dalpha_delta_tau3+x_4*dalpha_delta_tau4 +x_1*x_2*dalpha_delta_tau12+ x_1*x_4*dalpha_delta_tau14
        dalpha_delta_delta=x_1*dalpha_delta_delta1+x_2*dalpha_delta_delta2+ x_3*dalpha_delta_delta3 + x_4*dalpha_delta_delta4  +  x_1*x_2*dalpha_delta_delta12 +x_1*x_4*dalpha_delta_delta14
        
        tau1=component_params1['Tc']/T
        tau2=component_params2['Tc']/T
        tau3=component_params3['Tc']/T
        tau4=component_params4['Tc']/T
        #unlike the pure component, its just easier to do things in terms of moles,
        # then convert to per/kg, then ergs for DeBoers stuff...
        
        R_ideal12=component_params12['R']#/component_params12['M_amu']
        Ri1=component_params1['R']#/component_params1['M_amu']
        Ri2=component_params2['R']#/component_params2['M_amu']
        Ri3=component_params3['R']
        Ri4=component_params4['R']
        ave_R=(x_1*Ri1+x_2*Ri2 +x_3*Ri3+x_4*Ri4)/(x_1+x_2+x_3+x_4)      
        cv_r=-1.0*ave_R*(tau*tau*(dalpha_tau_tau))
        cv_o=x_1*(-1.0)*Ri1*(tau1*tau1*dalpha_o_tau_tau1)+x_2*(-1.0)*Ri2*(tau2*tau2*dalpha_o_tau_tau2)+x_3*(-1.0)*Ri3*(tau3*tau3*dalpha_o_tau_tau3)+x_4*(-1.0)*Ri4*(tau4*tau4*dalpha_o_tau_tau4)
       
        cp_mol=cv_r+cv_o+ave_R*power((1.0+delta*dalpha_delta-delta*tau*dalpha_delta_tau),2)/(1.0+2.0*delta*dalpha_delta+delta*delta*dalpha_delta_delta)
        #(the above gives J/(mol K)
        # now convert to kJ/(kg-K)
        M_ave=x_1*component_params1['M_amu']+x_2*component_params2['M_amu']+x_3*component_params3['M_amu']+x_4*component_params4['M_amu']
        cp=cp_mol/M_ave#component_params12['M_amu']
        #default is in kJ/kg-K convert to erg/K/mol for DeBoer's TCM
	#  kJ    1000J    1e7 ergs    1 kg    M_amu g      erg
	#----- * ----- * --------- * ------ * --------  = ------- 
	# kg*K    kJ      1J         1000 g    mol         K mol
	
	cp=cp*M_ave*1e7
	
        return cp
