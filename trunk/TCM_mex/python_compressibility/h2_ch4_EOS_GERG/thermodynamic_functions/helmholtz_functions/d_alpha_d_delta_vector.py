def d_alpha_d_delta_vector(tau,delta,parameters):
	from numpy import exp, power,meshgrid,sum,zeros,shape
	
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
	
	t1=zeros([len(parameters['N_i']),len(tau)])
	t1_v=zeros(len(tau))

        t1=(N_i_m[power_terms_wo_exp][:]*d_m[power_terms_wo_exp][:]*power(delta_m[power_terms_wo_exp][:],d_m[power_terms_wo_exp][:]-1.0)*power(tau_m[power_terms_wo_exp][:],t_m[power_terms_wo_exp][:]))
        t1_v=t1.sum(axis=0)
	
        #print shape(N_i_m),power_terms_w_exp,shape(c_m)
	t2=N_i_m[power_terms_w_exp][:]*exp(-power(delta_m[power_terms_w_exp][:],c_m[power_terms_w_exp][:]))*(power(delta_m[power_terms_w_exp][:],d_m[power_terms_w_exp][:]-1.0)*power(tau_m[power_terms_w_exp][:],t_m[power_terms_w_exp][:])*(d_m[power_terms_w_exp][:]-c_m[power_terms_w_exp][:]*power(delta_m[power_terms_w_exp][:],c_m[power_terms_w_exp][:])))          
        t2_v=t2.sum(axis=0)

        if(parameters['n_gaussian_terms']>0):
                m,n=shape(t1)
                p1=zeros([m,n])
                p2=zeros([m,n])
                dummy,alpha_m=meshgrid(delta,parameters['phi_i'])
	        dummy,beta_m=meshgrid(delta,parameters['Beta_i'])
	        dummy,gamma_m=meshgrid(delta,parameters['gamma_i'])
	        dummy,eps_m=meshgrid(delta,parameters['D_i'])

                p1=exp(-alpha_m[gaussian_terms][:]*((delta_m[gaussian_terms][:]-eps_m[gaussian_terms][:])*(delta_m[gaussian_terms][:]-eps_m[gaussian_terms][:])) - beta_m[gaussian_terms][:]*((tau_m[gaussian_terms][:]-gamma_m[gaussian_terms][:])*(tau_m[gaussian_terms][:]-gamma_m[gaussian_terms][:])))
                 
                p2=(d_m[gaussian_terms][:]/delta_m[gaussian_terms][:]) -2.0*alpha_m[gaussian_terms][:]*(delta_m[gaussian_terms][:]-eps_m[gaussian_terms][:])
                t3=N_i_m[gaussian_terms][:]*power(delta,d_m[gaussian_terms][:])*power(tau_m[gaussian_terms][:],t_m[gaussian_terms][:])*p1*p2
                t3_v=t3.sum(axis=0)
        else:
                t3_v[0]=0.0        
        if(parameters['n_critical_terms']>0):
                dummy,a_m=meshgrid(delta,parameters['RES_a'])
	        dummy,b_m=meshgrid(delta,parameters['RES_b'])
	        dummy,B_m=meshgrid(delta,parameters['RES_B'])
	        dummy,C_m=meshgrid(delta,parameters['RES_C'])
	        #dummy,C_m=meshgrid(tau,parameters['RES_D'])
	        dummy,A_m=meshgrid(delta,parameters['RES_A'])
	        dummy,D_m=meshgrid(delta,parameters['RES_D'])
	
                THETA = (1.0 - tau_m[critical_terms][:]) + A_m[critical_terms][:]*power((delta_m[critical_terms][:] - 1.0)*(delta_m[critical_terms][:] - 1.0),0.5/beta_m[critical_terms][:])
                DELTA = THETA*THETA + B_m[critical_terms][:]*power((delta_m[critical_terms][:] - 1.0)*(delta_m[critical_terms][:] - 1.0),a_m[critical_terms][:])
                PSI   = exp(-C_m[critical_terms][:]*((delta_m[critical_terms][:] - 1.0)*(delta_m[critical_terms][:] - 1.0)) - D_m[critical_terms][:]*((tau_m[critical_terms][:] - 1.0)*(tau_m[critical_terms][:] - 1.0)))
                D_DELTA_delta = (delta_m[critical_terms][:] - 1.0)*(A_m[critical_terms][:]*THETA*(2.0/beta_m[critical_terms][:])*power((delta_m[critical_terms][:] - 1.0)*(delta_m[critical_terms][:] - 1.0),0.5/beta_m[critical_terms][:] - 1.0) + 2.0*B_m[critical_terms][:]*a_m[critical_terms][:]*power((delta_m[critical_terms][:] - 1.0)*(delta_m[critical_terms][:] - 1.0),a_m[critical_terms][:]-1.0))
                D_DELTA_bi_delta = b_m[critical_terms][:]*power(DELTA,b_m[critical_terms][:]-1.0)*D_DELTA_delta
                D_PSI_delta = -2.0*C_m[critical_terms][:]*(delta_m[critical_terms][:] - 1.0)*PSI
       
                t4 = N_i_m[critical_terms][:]*(power(DELTA,b_m[critical_terms][:])*(PSI + delta_m[critical_terms][:]*D_PSI_delta) + D_DELTA_bi_delta*delta_m[critical_terms][:]*PSI)
        
                t4_v=t4.sum(axis=0)
        else:
                t4_v[0]=0.0
                                     
        return t1_v+t2_v+t3_v+t4_v             

   
