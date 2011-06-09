def ideal_helmholtz_from_coef_dtau(delta,tau,ideal_n,ideal_gamma):
    from math import pow,exp
    t1 = 0.0
    for i in range(3,len(ideal_gamma)):
        t1 += ideal_n[i]*ideal_gamma[i]*(pow(1.0 - exp(-ideal_gamma[i]*tau),-1.0) - 1.0)
    return ideal_n[1] + ideal_n[2]/tau + t1
