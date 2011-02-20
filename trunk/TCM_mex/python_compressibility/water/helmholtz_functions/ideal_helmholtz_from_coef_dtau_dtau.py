def ideal_helmholtz_from_coef_dtau_dtau(delta,tau,ideal_n,ideal_gamma):
    from math import exp,pow
    t1 = 0.0
    for i in range(3,len(ideal_gamma)):
        t1 += ideal_n[i]*(ideal_gamma[i]**2)*exp(-ideal_gamma[i]*tau)*pow(1.0 - exp(-ideal_gamma[i]*tau),-2.0)
    return -ideal_n[2]/(tau*tau) - t1
