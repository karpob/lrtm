def d_alpha_d_delta_d_delta_vector(tau,delta,parameters):
    #from math import exp, pow
    from numpy import exp,power,zeros,meshgrid,shape
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
    
    t1=N_i_m[ power_terms_wo_exp][:]*d_m[power_terms_wo_exp][:]*(d_m[power_terms_wo_exp][:]-1.0)*power(delta_m[power_terms_wo_exp][:],d_m[power_terms_wo_exp][:]-2.0)*power(tau_m[power_terms_wo_exp][:],t_m[power_terms_wo_exp][:])
    t1_v=t1.sum(axis=0)
    
    p1=power(delta_m[power_terms_w_exp][:],d_m[power_terms_w_exp][:]-2.0)*power(tau_m[power_terms_w_exp][:],t_m[power_terms_w_exp][:])
   
    p2=(d_m[power_terms_w_exp][:]-c_m[power_terms_w_exp][:]*power(delta_m[power_terms_w_exp][:],c_m[power_terms_w_exp][:]))*(d_m[power_terms_w_exp][:]-1.0-(c_m[power_terms_w_exp][:])*power(delta_m[power_terms_w_exp][:],c_m[power_terms_w_exp][:]))
   
    t2=N_i_m[power_terms_w_exp][:]*exp(-power(delta_m[power_terms_w_exp][:],c_m[power_terms_w_exp][:]))*(p1*(p2-(c_m[power_terms_w_exp][:]*c_m[power_terms_w_exp][:])*power(delta_m[power_terms_w_exp][:],c_m[power_terms_w_exp][:])))
    t2_v=t2.sum(axis=0)
   
              
    if(parameters['n_gaussian_terms']>0):
                m,n=shape(t1)
                p1=zeros([m,n])
                p2=zeros([m,n])
                dummy,alpha_m=meshgrid(delta,parameters['phi_i'])
	        dummy,beta_m=meshgrid(delta,parameters['Beta_i'])
	        dummy,gamma_m=meshgrid(delta,parameters['gamma_i'])
	        dummy,eps_m=meshgrid(delta,parameters['D_i'])
                # p1 = exp(-RES_alpha[i]*((delta - RES_eps[i])*(delta - RES_eps[i])) - RES_beta[i]*((tau-RES_gamma[i])*(tau-RES_gamma[i])))
                pp1= exp(-alpha_m[gaussian_terms][:]*((delta_m[gaussian_terms][:] - eps_m[gaussian_terms][:])*(delta_m[gaussian_terms][:] - eps_m[gaussian_terms][:])) - beta_m[gaussian_terms][:]*((tau_m[gaussian_terms][:]-gamma_m[gaussian_terms][:])*(tau_m[gaussian_terms][:]-gamma_m[gaussian_terms][:])))
                pp2=2.0*alpha_m[gaussian_terms][:]*power(delta_m[gaussian_terms][:],d_m[gaussian_terms][:])
                
                pp3=4.0*(alpha_m[gaussian_terms][:]*alpha_m[gaussian_terms][:])*power(delta_m[gaussian_terms][:],d_m[gaussian_terms][:])*((delta_m[gaussian_terms][:]-eps_m[gaussian_terms][:])*(delta_m[gaussian_terms][:]-eps_m[gaussian_terms][:]))
                pp4=4.0*d_m[gaussian_terms][:]*alpha_m[gaussian_terms][:]*power(delta_m[gaussian_terms][:],d_m[gaussian_terms][:] -1.0)*(delta_m[gaussian_terms][:]-eps_m[gaussian_terms][:])
                pp5=d_m[gaussian_terms][:]*(d_m[gaussian_terms][:]-1.0)*power(delta_m[gaussian_terms][:],d_m[gaussian_terms][:]-2.0)
                t3=N_i_m[gaussian_terms][:]*power(tau_m[gaussian_terms][:],t_m[gaussian_terms][:])*pp1*(-pp2+pp3-pp4+pp5)
                t3_v=t3.sum(axis=0)
    
    if(parameters['n_critical_terms']>0):
                dummy,a_m=meshgrid(delta,parameters['RES_a'])
	        dummy,b_m=meshgrid(delta,parameters['RES_b'])
	        dummy,B_m=meshgrid(delta,parameters['RES_B'])
	        dummy,C_m=meshgrid(delta,parameters['RES_C'])
	        #dummy,C_m=meshgrid(tau,parameters['RES_D'])
	        dummy,A_m=meshgrid(delta,parameters['RES_A'])
	        dummy,D_m=meshgrid(delta,parameters['RES_D'])
                THETA=(1.0-tau_m[critical_terms][:])+A_m[critical_terms][:]*power((delta_m[critical_terms][:]-1.0)*(delta_m[critical_terms][:]-1.0),0.5/beta_m[critical_terms][:])
                
                DELTA = THETA*THETA + B_m[critical_terms][:]*power((delta_m[critical_terms][:] - 1.0)*(delta_m[critical_terms][:] - 1.0),a_m[critical_terms][:])
                #PSI   = exp(-RES_C[i]*((delta - 1.0)*(delta - 1.0)) - RES_D[i]*((tau - 1.0)*(tau - 1.0)))
                PSI   = exp(-C_m[critical_terms][:]*((delta_m[critical_terms][:] - 1.0)*(delta_m[critical_terms][:] - 1.0)) - D_m[critical_terms][:]*((tau_m[critical_terms][:] - 1.0)*(tau_m[critical_terms][:] - 1.0)))
                #PSI   = exp(-C_m[critical_terms][:]*((delta_m[critical_terms][:] - 1.0)*(delta_m[critical_terms][:] - 1.0)) - D_m[critical_terms][:]*((tau_m[critical_terms][:] - 1.0)*(tau_m[critical_terms][:] - 1.0)))
                #print PSI
                D_DELTA_delta = (delta_m[critical_terms][:] - 1.0)*(A_m[critical_terms][:]*THETA*(2.0/beta_m[critical_terms][:])*power((delta_m[critical_terms][:] - 1.0)*(delta_m[critical_terms][:] - 1.0),0.5/beta_m[critical_terms][:] - 1.0) + 2.0*B_m[critical_terms][:]*a_m[critical_terms][:]*power((delta_m[critical_terms][:] - 1.0)*(delta_m[critical_terms][:] - 1.0),a_m[critical_terms][:]-1.0))
                D_DELTA_bi_delta = b_m[critical_terms][:]*power(DELTA,b_m[critical_terms][:]-1.0)*D_DELTA_delta
                D_PSI_delta = -2.0*C_m[critical_terms][:]*(delta_m[critical_terms][:] - 1.0)*PSI
                D2_PSI_delta = 2.0*C_m[critical_terms][:]*PSI*(2.0*C_m[critical_terms][:]*(delta_m[critical_terms][:] - 1.0)*(delta_m[critical_terms][:] - 1.0) - 1.0)
                ppp1 = 4.0*B_m[critical_terms][:]*a_m[critical_terms][:]*(a_m[critical_terms][:] - 1.0)*power((delta_m[critical_terms][:] - 1.0)*(delta_m[critical_terms][:] - 1.0),a_m[critical_terms][:] - 2.0)
                ppp2 = 2.0*(A_m[critical_terms][:]*A_m[critical_terms][:])*power(beta_m[critical_terms][:],-2.0)*power(power((delta_m[critical_terms][:] - 1.0)*(delta_m[critical_terms][:] - 1.0),0.5/beta_m[critical_terms][:] - 1.0),2)
                ppp3 = A_m[critical_terms][:]*THETA*(4.0/beta_m[critical_terms][:])*(0.5/beta_m[critical_terms][:] - 1.0)*pow((delta_m[critical_terms][:] - 1.0)*(delta_m[critical_terms][:] - 1.0),0.5/beta_m[critical_terms][:] - 2.0)
                D2_DELTA_delta = (D_DELTA_delta/(delta_m[critical_terms][:] - 1.0)) + ((delta_m[critical_terms][:] - 1.0)*(delta_m[critical_terms][:] - 1.0))*(ppp1 + ppp2 + ppp3)
                D2_DELTA_bi_delta = b_m[critical_terms][:]*(pow(DELTA,b_m[critical_terms][:] - 1.0)*D2_DELTA_delta + (b_m[critical_terms][:] - 1.0)*power(DELTA,b_m[critical_terms][:] - 2.0)*power(D_DELTA_delta,2.0))        
                pppp1 = pow(DELTA,b_m[critical_terms][:])*(2.0*D_PSI_delta + delta*D2_PSI_delta)       
                pppp2 = 2.0*D_DELTA_bi_delta*(PSI + delta_m[critical_terms][:]*D_PSI_delta)
                pppp3 = D2_DELTA_bi_delta*delta_m[critical_terms][:]*PSI 
      
                t4 = N_i_m[critical_terms][:]*(pppp1 + pppp2 + pppp3)
                t4_v=t4.sum(axis=0)
                
    return t1_v+t2_v+t3_v+t4_v
