def d_alpha_d_delta(tau,delta,RES_n,RES_t,RES_d,RES_c,RES_alpha,RES_beta,RES_gamma,RES_eps,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A):
	from math import exp, pow
	t1=0.0
	t2=0.0
	
	for i in range(0,n_power_terms):
	        
	        if(abs(float(RES_c[i]))>0.0):         
	                t2 +=RES_n[i]*exp(-pow(delta,RES_c[i]))*(pow(delta,RES_d[i]-1.0)*pow(tau,RES_t[i])*(RES_d[i] - RES_c[i]*pow(delta,RES_c[i])))#RES_n[i]*pow(tau,RES_t[i])*exp(-RES_gamma[i]*pow(delta,RES_c[i]))*(pow(delta,RES_d[i]-1)*RES_d[i]-pow(delta,RES_d[i]+RES_c[i]-1)*RES_gamma[i]*RES_c[i]) #(RES_n[i]*exp(-RES_gamma[i]*pow(delta,RES_c[i]))*(pow(delta,RES_d[i]-1.0)*pow(tau,RES_t[i])*(RES_d[i] - RES_gamma[i]*RES_c[i]*pow(delta,RES_c[i]))))
	             #  print 'exp term',i+1,RES_n[i]*exp(-pow(delta,RES_c[i]))*(pow(delta,RES_d[i]-1.0)*pow(tau,RES_t[i])*(RES_d[i] - RES_c[i]*pow(delta,RES_c[i])))
	                
                else:
                       # if(abs(RES_n[i])>0.0):
                        t1 +=(RES_n[i]*RES_d[i]*pow(delta,RES_d[i]-1.0)*pow(tau,RES_t[i]))
                        #print t1
	                 #if(float(RES_t[i])<float(0.0)):
	                #print 'power_term',i+1,RES_t[i],RES_d[i],RES_n[i]#(RES_n[i]*RES_d[i]*pow(delta,RES_d[i]-1.0)*pow(tau,RES_t[i]))#*RES_d[i]*pow(delta,RES_d[i]-1.0)*pow(tau,RES_t[i])
	t3=0.0                   
        for i in range(n_power_terms,n_power_terms+n_gaussian_terms):
                p1 = exp(-RES_alpha[i]*((delta - RES_eps[i])*(delta - RES_eps[i])) - RES_beta[i]*((tau - RES_gamma[i])*(tau - RES_gamma[i])))
                p2 = (RES_d[i]/delta) - 2.0*RES_alpha[i]*(delta - RES_eps[i])
                t3 += RES_n[i]*pow(delta,RES_d[i])*pow(tau,RES_t[i])*p1*p2
                
                
        t4=0.0
        for i in range(n_power_terms+n_gaussian_terms,n_power_terms+n_gaussian_terms+n_critical_terms):
                THETA = (1.0 - tau) + RES_A[i]*pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i])
                DELTA = THETA*THETA + RES_B[i]*pow((delta - 1.0)*(delta - 1.0),RES_a[i])
                PSI   = exp(-RES_C[i]*((delta - 1.0)*(delta - 1.0)) - RES_D[i]*((tau - 1.0)*(tau - 1.0)))

                D_DELTA_delta = (delta - 1.0)*(RES_A[i]*THETA*(2.0/RES_beta[i])*pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i] - 1.0) + 2.0*RES_B[i]*RES_a[i]*pow((delta - 1.0)*(delta - 1.0),RES_a[i]-1.0))
                D_DELTA_bi_delta = RES_b[i]*pow(DELTA,RES_b[i]-1.0)*D_DELTA_delta
                D_PSI_delta = -2.0*RES_C[i]*(delta - 1.0)*PSI
        
                t4 += RES_n[i]*(pow(DELTA,RES_b[i])*(PSI + delta*D_PSI_delta) + D_DELTA_bi_delta*delta*PSI)
        #print t1,t2,t3,t4               
        return t1+t2+t3+t4             

   
