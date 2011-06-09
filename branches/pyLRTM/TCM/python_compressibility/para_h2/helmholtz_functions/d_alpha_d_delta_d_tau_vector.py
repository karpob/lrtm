def d_alpha_d_delta_d_tau_vector(tau,delta,parameters):
        from numpy import exp, power,zeros,shape,meshgrid
        t1=zeros([len(parameters['N_i']),len(tau)])
        t1_v=zeros(len(tau))
        t2=zeros([len(parameters['N_i']),len(tau)])
        t2_v=zeros(len(tau))
        t3=zeros([len(parameters['N_i']),len(tau)])
        t3_v=zeros(len(tau))
        t4_v=zeros(len(tau))
        power_terms_wo_exp=range(0,parameters['n_power_terms_wo_exp'])
        power_terms_w_exp=range(parameters['n_power_terms_wo_exp'],parameters['n_power_terms_wo_exp']+parameters['n_power_terms_w_exp'])
        n_power_terms=parameters['n_power_terms_wo_exp']+parameters['n_power_terms_w_exp']
        gaussian_terms=range(n_power_terms,n_power_terms+parameters['n_gaussian_terms'])
        critical_terms=range(n_power_terms+parameters['n_gaussian_terms'],n_power_terms+parameters['n_gaussian_terms']+parameters['n_critical_terms'])
        #Resize vectors to matricies
        tau_m,N_i_m=meshgrid(tau,parameters['N_i'])
        delta_m,dummy2=meshgrid(delta,parameters['N_i'])
        dummy,t_m=meshgrid(delta,parameters['t_i'])
        dummy,d_m=meshgrid(delta,parameters['d_i'])
        dummy,c_m=meshgrid(delta,parameters['p_i'])
        
        t1 = N_i_m[power_terms_wo_exp][:]*d_m[power_terms_wo_exp][:]*t_m[power_terms_wo_exp][:]*power(delta,(d_m[power_terms_wo_exp][:] - 1.0))*power(tau_m[power_terms_wo_exp][:],(t_m[power_terms_wo_exp][:] - 1.0))
        t1_v=t1.sum(axis=0)
        
        p1  = (d_m[power_terms_w_exp][:] - c_m[power_terms_w_exp][:]*power(delta_m[power_terms_w_exp][:],c_m[power_terms_w_exp][:]))*exp((-power(delta_m[power_terms_w_exp][:],c_m[power_terms_w_exp][:])))
        t2 = N_i_m[power_terms_w_exp][:]*t_m[power_terms_w_exp][:]*power(delta_m[power_terms_w_exp][:],d_m[power_terms_w_exp][:] - 1.0)*power(tau_m[power_terms_w_exp][:],t_m[power_terms_w_exp][:] - 1.0)*p1
        t2_v=t2.sum(axis=0)
        if(parameters['n_gaussian_terms']>0):
                m,n=shape(t1)
                p1=zeros([m,n])
                p2=zeros([m,n])
                dummy,alpha_m=meshgrid(delta,parameters['phi_i'])
	        dummy,beta_m=meshgrid(delta,parameters['Beta_i'])
	        dummy,gamma_m=meshgrid(delta,parameters['gamma_i'])
	        dummy,eps_m=meshgrid(delta,parameters['D_i'])
                d1=delta_m[gaussian_terms][:]-eps_m[gaussian_terms][:]
                t1=tau_m[gaussian_terms][:]-gamma_m[gaussian_terms][:]
                e1 = -1.0*alpha_m[gaussian_terms][:]*d1*d1 - beta_m[gaussian_terms][:]*t1*t1
                
		
                f1 = t_m[gaussian_terms][:] - 2*beta_m[gaussian_terms][:]*tau_m[gaussian_terms][:]*(tau_m[gaussian_terms][:] - gamma_m[gaussian_terms][:])
		g1 = d_m[gaussian_terms][:] - 2*alpha_m[gaussian_terms][:]*delta_m[gaussian_terms][:]*(delta_m[gaussian_terms][:] - eps_m[gaussian_terms][:])

		sumz = N_i_m[gaussian_terms][:] * f1 * power(tau_m[gaussian_terms][:],t_m[gaussian_terms][:]-1) * g1 * power(delta_m[gaussian_terms][:],d_m[gaussian_terms][:]-1) * exp(e1)
		t3_v=sumz.sum(axis=0)
                
     
        if(parameters['n_critical_terms']>0):
                dummy,a_m=meshgrid(delta,parameters['RES_a'])
	        dummy,b_m=meshgrid(delta,parameters['RES_b'])
	        dummy,B_m=meshgrid(delta,parameters['RES_B'])
	        dummy,C_m=meshgrid(delta,parameters['RES_C'])
	        #dummy,C_m=meshgrid(tau,parameters['RES_D'])
	        dummy,A_m=meshgrid(delta,parameters['RES_A'])
	        dummy,D_m=meshgrid(delta,parameters['RES_D'])
                d1 = delta_m[critical_terms][:] - 1.0
                t1 = tau_m[critical_terms][:] - 1.0
                psi = exp(-C_m[critical_terms][:]*d1*d1 - D_m[critical_terms][:]*t1*t1)
                theta = 1 - tau_m[critical_terms][:] + A_m[critical_terms][:]*power(d1*d1, (1/(2*beta_m[critical_terms][:])))
                Delta = theta*theta + B_m[critical_terms][:]*power(d1*d1, (a_m[critical_terms][:]))
                dDeltabidtau = -2*theta*b_m[critical_terms][:]*power(Delta, (b_m[critical_terms][:]-1))
                dpsidtau = -2*D_m[critical_terms][:]*t1*psi
                dDeltaddelta = d1*(A_m[critical_terms][:]*theta*2/beta_m[critical_terms][:]*pow(d1*d1, (1/(2*beta_m[critical_terms][:]) - 1)) + 2*B_m[critical_terms][:]*a_m[critical_terms][:]*power(d1*d1, (a_m[critical_terms][:] - 1)))
                d2psiddeltadtau = 4*C_m[critical_terms][:]*D_m[critical_terms][:]*d1*t1*psi
                dDeltabiddelta = b_m[critical_terms][:]*power(Delta, (b_m[critical_terms][:]-1))*dDeltaddelta
                dpsiddelta = -2*C_m[critical_terms][:]*d1*psi
                d2Deltabiddeltadtau = -A_m[critical_terms][:]*b_m[critical_terms][:]*2/beta_m[critical_terms][:]*power(Delta, (b_m[critical_terms][:]-1))*d1*power(d1*d1, (1/(2*beta_m[critical_terms][:])-1)) -2*theta*b_m[critical_terms][:]*(b_m[critical_terms][:]-1)*power(Delta, (b_m[critical_terms][:]-2))*dDeltaddelta
                t4 = N_i_m[critical_terms][:]*(power(Delta, b_m[critical_terms][:])*(dpsidtau + delta_m[critical_terms][:]*d2psiddeltadtau) + delta_m[critical_terms][:]*dDeltabiddelta*dpsidtau +dDeltabidtau*(psi + delta_m[critical_terms][:]*dpsiddelta) +d2Deltabiddeltadtau*delta_m[critical_terms][:]*psi)
                t4_v=t4.sum(axis=0)        
    
        return t1_v+t2_v+t3_v+t4_v
