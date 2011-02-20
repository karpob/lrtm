def d_alpha_d_delta_d_delta(tau,delta,RES_n,RES_t,RES_d,RES_c,RES_alpha,RES_beta,RES_gamma,RES_eps,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A):
    from math import exp, pow
    t1 = 0.0
    
    for i in range(0,n_power_terms):
        if(abs(float(RES_c[i]))>0.0):
                p1 = pow(delta,RES_d[i] - 2.0)*pow(tau,RES_t[i])
                p2 = (RES_d[i] - RES_c[i]*pow(delta,RES_c[i]))*(RES_d[i] - 1.0 - (RES_c[i])*pow(delta,RES_c[i]))
                t1 += RES_n[i]*exp(-pow(delta,RES_c[i]))*(p1*(p2 - (RES_c[i]*RES_c[i])*pow(delta,RES_c[i])))
                
        else:
                t1 += RES_n[i]*RES_d[i]*(RES_d[i]-1.0)*pow(delta,RES_d[i]-2.0)*pow(tau,RES_t[i])
        
    t3 = 0.0
    #print t1
    for i in range(n_power_terms,n_power_terms+n_gaussian_terms):
        p1 = exp(-RES_alpha[i]*((delta - RES_eps[i])*(delta - RES_eps[i])) - RES_beta[i]*((tau-RES_gamma[i])*(tau-RES_gamma[i])))
        p2 = 2.0*RES_alpha[i]*pow(delta,RES_d[i])
        p3 = 4.0*(RES_alpha[i]*RES_alpha[i])*pow(delta,RES_d[i])*((delta - RES_eps[i])*(delta - RES_eps[i]))
        p4 = 4.0*RES_d[i]*RES_alpha[i]*pow(delta,RES_d[i] - 1.0)*(delta - RES_eps[i])
        p5 = RES_d[i]*(RES_d[i] - 1.0)*pow(delta,RES_d[i] - 2.0)
        t1 += RES_n[i]*pow(tau,RES_t[i])*p1*(-p2 + p3 - p4 + p5)
    
    t4 = 0.0
    for i in range(n_power_terms+n_gaussian_terms,n_power_terms+n_gaussian_terms+n_critical_terms):
        THETA = (1.0 - tau) + RES_A[i]*pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i])
        
        DELTA = THETA*THETA + RES_B[i]*pow((delta - 1.0)*(delta - 1.0),RES_a[i])
        
        PSI   = exp(-RES_C[i]*((delta - 1.0)*(delta - 1.0)) - RES_D[i]*((tau - 1.0)*(tau - 1.0)))
        
        D_DELTA_delta = (delta - 1.0)*(RES_A[i]*THETA*(2.0/RES_beta[i])*pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i] - 1.0) + 2.0*RES_B[i]*RES_a[i]*pow((delta - 1.0)*(delta - 1.0),RES_a[i]-1.0))
        D_DELTA_bi_delta = RES_b[i]*pow(DELTA,RES_b[i]-1.0)*D_DELTA_delta
        D_PSI_delta = -2.0*RES_C[i]*(delta - 1.0)*PSI
        
        D2_PSI_delta = 2.0*RES_C[i]*PSI*(2.0*RES_C[i]*(delta - 1.0)*(delta - 1.0) - 1.0)

        p1 = 4.0*RES_B[i]*RES_a[i]*(RES_a[i] - 1.0)*pow((delta - 1.0)*(delta - 1.0),RES_a[i] - 2.0)
        p2 = 2.0*(RES_A[i]*RES_A[i])*pow(RES_beta[i],-2.0)*pow(pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i] - 1.0),2)
        p3 = RES_A[i]*THETA*(4.0/RES_beta[i])*(0.5/RES_beta[i] - 1.0)*pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i] - 2.0)

        D2_DELTA_delta = (D_DELTA_delta/(delta - 1.0)) + ((delta - 1.0)*(delta - 1.0))*(p1 + p2 + p3)
        D2_DELTA_bi_delta = RES_b[i]*(pow(DELTA,RES_b[i] - 1.0)*D2_DELTA_delta + (RES_b[i] - 1.0)*pow(DELTA,RES_b[i] - 2.0)*pow(D_DELTA_delta,2.0))
        
        #p1 - p5 can be discarded now
        p1 = pow(DELTA,RES_b[i])*(2.0*D_PSI_delta + delta*D2_PSI_delta)
        p2 = 2.0*D_DELTA_bi_delta*(PSI + delta*D_PSI_delta)
        p3 = D2_DELTA_bi_delta*delta*PSI
        
        #finally
        t1 += RES_n[i]*(p1 + p2 + p3)
        #print t1
    return t1
