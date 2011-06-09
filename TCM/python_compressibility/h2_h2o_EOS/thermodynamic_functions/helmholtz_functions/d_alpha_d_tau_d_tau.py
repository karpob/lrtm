def d_alpha_d_tau_d_tau(tau,delta,RES_n,RES_t,RES_d,RES_c,RES_alpha,RES_beta,RES_gamma,RES_eps,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A):
    from math import exp, pow
    t1 = 0.0
    t2=0.0
    for i in range(0,n_power_terms):
        if(abs(float(RES_c[i]))>0.0):
                t2 += RES_n[i]*RES_t[i]*(RES_t[i] - 1.0)*pow(delta,RES_d[i])*pow(tau,RES_t[i] - 2.0)*exp(-pow(delta,RES_c[i]))
        else:
                t1 += RES_n[i]*RES_t[i]*(RES_t[i] - 1.0)*pow(delta,RES_d[i])*pow(tau,RES_t[i] - 2.0)
    t3 = 0.0
    for i in range(n_power_terms,n_power_terms+n_gaussian_terms):
        p1 = -RES_alpha[i]*((delta - RES_eps[i])*(delta - RES_eps[i])) - RES_beta[i]*((tau - RES_gamma[i])*(tau - RES_gamma[i]))
        p2 = pow((RES_t[i]/tau - 2.0*RES_beta[i]*(tau - RES_gamma[i])),2.0) - RES_t[i]/(tau*tau) - 2.0*RES_beta[i]
        t3 += RES_n[i]*pow(delta,RES_d[i])*pow(tau,RES_t[i])*exp(p1)*p2

    t4 = 0.0
    for i in range(n_power_terms+n_gaussian_terms,n_power_terms+n_gaussian_terms+n_critical_terms):
        THETA = (1.0 - tau) + RES_A[i]*pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i])
        
        DELTA = THETA*THETA + RES_B[i]*pow((delta - 1.0)*(delta - 1.0),RES_a[i])
        
        PSI   = exp(-RES_C[i]*((delta - 1.0)*(delta - 1.0)) - RES_D[i]*((tau - 1.0)*(tau - 1.0)))
        
        D_DELTA_bi_tau = -2.0*THETA*RES_b[i]*pow(DELTA,RES_b[i] - 1.0)
        
        D_PSI_tau = -2.0*RES_D[i]*(tau - 1.0)*PSI
         
        D2_DELTA_bi_tau = 2.0*RES_b[i]*pow(DELTA,RES_b[i] - 1.0) + 4.0*(THETA*THETA)*RES_b[i]*(RES_b[i] - 1.0)*pow(DELTA,RES_b[i] - 2.0)
        
        D2_PSI_tau = (2*RES_D[i]*((tau - 1.0)*(tau - 1.0)) - 1.0)*2.0*RES_D[i]*PSI
        
        t4 += RES_n[i]*delta*(D2_DELTA_bi_tau*PSI + 2.0*D_DELTA_bi_tau*D_PSI_tau + pow(DELTA,RES_b[i])*D2_PSI_tau)
    return t1+t2+t3+t4
