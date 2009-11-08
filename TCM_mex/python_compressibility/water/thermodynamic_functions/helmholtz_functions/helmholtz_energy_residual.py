def helmholtz_energy_residual(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A):
	from math import exp, pow,sqrt
	from pylab import shape,zeros
	alpha_1=0.0
	alpha_2=0.0
	alpha_3=0.0
        alpha_4=0.0
        t4=0.0

	#alpha_mat=zeros(len(N_i))
        for i in range(0,n_power_terms):
		if(abs(float(p_i[i]))!=0.0):
#			print 'exp fun',N_i[i],i+1,d_i[i],t_i[i],p_i[i]
			alpha_1=alpha_1+N_i[i]*pow(delta,float(d_i[i]))*pow(tau,t_i[i])*exp(-pow(delta,p_i[i]))
			
		else:
#			print 'power fun', N_i[i],i+1,d_i[i],t_i[i]
			alpha_2=alpha_2+N_i[i]*pow(delta,float(d_i[i]))*pow(tau,t_i[i])
			
			#print 'powerfun', i+1,N_i[i]*pow(delta,float(d_i[i]))*pow(tau,t_i[i])
	for i in range(n_power_terms,n_power_terms+n_gaussian_terms):
	       #print d_i[i]
		alpha_3=alpha_3+N_i[i]*pow(delta,float(d_i[i]))*pow(tau,float(t_i[i]))*exp(-phi_i[i]*pow(delta-D_i[i],2)-Beta_i[i]*pow(tau-gamma_i[i],2))
		
        for i in range(n_power_terms+n_gaussian_terms,n_power_terms+n_gaussian_terms+n_critical_terms):
                THETA = (1.0 - tau) + RES_A[i]*pow((delta - 1.0)*(delta - 1.0),0.5/Beta_i[i])
                DELTA = THETA*THETA + RES_B[i]*pow((delta - 1.0)*(delta - 1.0),RES_a[i])
                PSI   = exp(-RES_C[i]*((delta - 1.0)*(delta - 1.0)) - RES_D[i]*((tau - 1.0)*(tau - 1.0)))
                t4 += N_i[i]*pow(DELTA,RES_b[i])*delta*PSI
        alpha_4=t4
        #print 'my first alpha',alpha_mat[0]+alpha_mat[1]+alpha_mat[2]+alpha_mat[3]+alpha_mat[4]      
	return alpha_1+alpha_2+alpha_3+alpha_4
