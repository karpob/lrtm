def residual_pressure_from_enthalpy(density_arr,Pin,h2_N_i,h2_t_i,h2_d_i,h2_p_i,h2_phi_i,h2_Beta_i,h2_gamma_i,h2_D_i,h2_n_power_terms,h2_n_gaussian_terms,h2_n_critical_terms,h2_RES_a,h2_RES_b,h2_RES_B,h2_RES_C,h2_RES_D,h2_RES_A,h2_R,h2_Tc,h2_rho_c,h2_M_amu,h2o_N_i,h2o_t_i,h2o_d_i,h2o_p_i,h2o_phi_i,h2o_Beta_i,h2o_gamma_i,h2o_D_i,h2o_n_power_terms,h2o_n_gaussian_terms,h2o_n_critical_terms,h2o_RES_a,h2o_RES_b,h2o_RES_B,h2o_RES_C,h2o_RES_D,h2o_RES_A,h2o_R,h2o_Tc,h2o_rho_c,h2o_M_amu,h2o_h2_mix_N_i,h2o_h2_mix_t_i,h2o_h2_mix_d_i,h2o_h2_mix_p_i,h2o_h2_mix_phi_i,h2o_h2_mix_Beta_i,h2o_h2_mix_gamma_i,h2o_h2_mix_D_i,h2o_h2_mix_n_power_terms,h2o_h2_mix_n_gaussian_terms,h2o_h2_mix_n_critical_terms,h2o_h2_mix_RES_a,h2o_h2_mix_RES_b,h2o_h2_mix_RES_B,h2o_h2_mix_RES_C,h2o_h2_mix_RES_D,h2o_h2_mix_RES_A,h2o_h2_mix_R,h2o_h2_mix_Tc,h2o_h2_mix_rho_c,h2o_h2_mix_M_amu,beta_ij,phi_ij,sigma_ij,xi_ij,x_h2,x_h2o,T):
        from thermodynamic_functions.pressure import pressure
        from thermodynamic_functions.enthalpy import enthalpy
        from scale_tau_and_delta import scale_tau_and_delta
        from pylab import shape
        density=density_arr        
        h2_delta=(h2_M_amu/((x_h2o*h2o_M_amu)+(x_h2*h2_M_amu)))*(x_h2*density)/h2_rho_c # Need to scale by molecular weight?
        h2_tau=h2_Tc/T
        h2o_delta=(h2o_M_amu/((x_h2o*h2o_M_amu)+(x_h2*h2_M_amu)))*(x_h2o*density)/h2o_rho_c #Need to scale by molecular weight?
        h2o_tau=h2o_Tc/T
        [h2o_h2_mix_tau,h2o_h2_mix_delta,h2o_h2_mix_Tc,h2o_h2_mix_rho_c]=scale_tau_and_delta(T,density,x_h2,x_h2o,h2_Tc,h2o_Tc,h2_rho_c,h2o_rho_c,beta_ij,phi_ij,sigma_ij,xi_ij)
        
        pressure_h2_EOS_contribution=pressure(h2_tau,h2_delta,h2_N_i,h2_t_i,h2_d_i,h2_p_i,h2_phi_i,h2_Beta_i,h2_gamma_i,h2_D_i,h2_n_power_terms,h2_n_gaussian_terms,h2_n_critical_terms,h2_RES_a,h2_RES_b,h2_RES_B,h2_RES_C,h2_RES_D,h2_RES_A,h2_R,h2_Tc,h2_rho_c,h2_M_amu)
        pressure_h2o_EOS_contribution=pressure(h2o_tau,h2o_delta,h2o_N_i,h2o_t_i,h2o_d_i,h2o_p_i,h2o_phi_i,h2o_Beta_i,h2o_gamma_i,h2o_D_i,h2o_n_power_terms,h2o_n_gaussian_terms,h2o_n_critical_terms,h2o_RES_a,h2o_RES_b,h2o_RES_B,h2o_RES_C,h2o_RES_D,h2o_RES_A,h2o_R,h2o_Tc,h2o_rho_c,h2o_M_amu)                
        
        pressure_delta_h2o_h2_contribution=pressure(h2o_h2_mix_tau,h2o_h2_mix_delta,h2o_h2_mix_N_i,h2o_h2_mix_t_i,h2o_h2_mix_d_i,h2o_h2_mix_p_i,h2o_h2_mix_phi_i,h2o_h2_mix_Beta_i,h2o_h2_mix_gamma_i,h2o_h2_mix_D_i,h2o_h2_mix_n_power_terms,h2o_h2_mix_n_gaussian_terms,h2o_h2_mix_n_critical_terms,h2o_h2_mix_RES_a,h2o_h2_mix_RES_b,h2o_h2_mix_RES_B,h2o_h2_mix_RES_C,h2o_h2_mix_RES_D,h2o_h2_mix_RES_A,h2o_h2_mix_R,h2o_h2_mix_Tc,h2o_h2_mix_rho_c,h2o_h2_mix_M_amu)
        
        #enthalpy_delta_h2o_h2=enthalpy(h2o_h2_mix_tau,h2o_h2_mix_delta,h2o_h2_mix_ni,h2o_h2_mix_ti,h2o_h2_mix_vi,h2o_h2_mix_ui,h2o_h2_mix_R,h2o_h2_mix_ho,h2o_h2_mix_so,h2o_h2_mix_rho_c,h2o_h2_mix_Tc,h2o_h2_mix_n_ideal_gas_terms_pow,h2o_h2_mix_n_ideal_gas_terms_exp,h2o_h2_mix_To,h2o_h2_mix_rho_o,h2o_h2_mix_N_i,h2o_h2_mix_t_i,h2o_h2_mix_d_i,h2o_h2_mix_p_i,h2o_h2_mix_phi_i,h2o_h2_mix_Beta_i,h2o_h2_mix_gamma_i,h2o_h2_mix_D_i,h2o_h2_mix_n_power_terms,h2o_h2_mix_n_gaussian_terms,h2o_h2_mix_n_critical_terms,h2o_h2_mix_RES_a,h2o_h2_mix_RES_b,h2o_h2_mix_RES_B,h2o_h2_mix_RES_C,h2o_h2_mix_RES_D,h2o_h2_mix_RES_A,h2o_h2_mix_M_amu,h2o_h2_mix_ideal_n,h2o_h2_mix_ideal_gamma,'none')
                
        #enthalpy_excess=enthalpy_delta_h2o_h2
        Pcalc=pressure_h2_EOS_contribution+pressure_h2o_EOS_contribution+pressure_delta_h2o_h2_contribution
        err=Pin-Pcalc
        
        return err
