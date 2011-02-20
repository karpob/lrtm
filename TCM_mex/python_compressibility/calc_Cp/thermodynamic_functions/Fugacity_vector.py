def Fugacity_vector(tau,delta,eqn_type,parameter_set1,parameter_set2,parameter_set12,x_1,x_2,rho_c12,Tc12):
        # this computes the fugacity for a mixture using finitite difference approximation
        # using "central difference"
        step_in_moles=1e-4 # taken as mol step in REFPROP
        from helmholtz_functions.helmholtz_energy_residual_vector import helmholtz_energy_residual_vector
        from numpy import zeros
        dalpha_delta=zeros(len(tau))
        rho=rho_c12*delta # get density back
        T=Tc12/tau # get Temperature back
        Rspecific=parameter_set12['R']/parameter_set12['M_amu']
        for i in range(0,len(x_1)):
                if(x_1[i]<step_in_moles):
                        delta_neg[i]=-x_1[i]/2.0
                        delta_pos[i]=-delta_neg[i]
                else:
                        delta_neg[i]=-step_in_moles
                        delta_pos[i]=step_in_moles                
        delta_plus_1=1.0/(1.0+delta_pos)
        delta_minus_1=1.0/(1.0+delta_neg)
        
        xplus=(x_1*delta_plus_1+delta_pos)*delta_plus_1
        xminus=(x_1*delta_minus_1+delta_neg)*delta_minus_1
        
        xplus2=(x_2*delta_plus_1+delta_pos)*delta_plus_1
        xminus2=(x_2*delta_minus_1+delta_neg)*delta_minus_1
        
        Dplus=rho*(1.0+delta_pos)
        Dminus=rho*(1.0+delta_neg)
        
        # call reducing function here with xplus
        [tau,delta,Tc,rho_c]=scale_tau_and_delta_kw(T,Dplus,xplus,xplus2,parameter_set1['Tc'],parameter_set2['Tc'],parameter_set1['rho_c'],parameter_set2['rho_c'],parameter_set12['BetaT'],parameter_set12['BetaV'],parameter_set12['GammaT'],parameter_set12['GammaV'])
        
        
        A_r1=helmholtz_energy_residual_vector(tau,delta,parameter_set1)
        A_r2=helmholtz_energy_residual_vector(tau,delta,parameter_set2)
        A_r12=helmholtz_energy_residual_vector(tau,delta,parameter_set12)
        Helm_plus=A_r1+A_r2+A_r12
        
        #call reducing function here with xminus
        [tau,delta,Tc,rho_c]=scale_tau_and_delta_kw(T,Dminus,xminus,xminus2,parameter_set1['Tc'],parameter_set2['Tc'],parameter_set1['rho_c'],parameter_set2['rho_c'],parameter_set12['BetaT'],parameter_set12['BetaV'],parameter_set12['GammaT'],parameter_set12['GammaV'])
        tau=Tc_mix/T
        delta=Dminus/rho_c_mix
        A_r1=helmholtz_energy_residual_vector(tau,delta,parameter_set1)
        A_r2=helmholtz_energy_residual_vector(tau,delta,parameter_set2)
        A_r12=helmholtz_energy_residual_vector(tau,delta,parameter_set12)
        Helm_minus=A_r1+A_r2+A_r12
        
        da_dn=((1.0+delta_pos)*Helm_plus-(1.0+delta_neg)*Helm_minus)/(delta_pos-delta_neg)
        fugacity=x_1*Rspecific*T*rho*exp(da_dn) #return in kPa                                           
        return fugacity 
