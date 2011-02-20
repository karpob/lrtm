def ideal_helmholtz_from_coef(delta,tau,ideal_n,ideal_gamma):
    t1 = 0.0
    from math import log,exp
    for i in range(3,len(ideal_gamma)):
        t1 += ideal_n[i]*log(1.0 - exp(-ideal_gamma[i]*tau))
        
    return log(delta) + ideal_n[0] + ideal_n[1]*tau + ideal_n[2]*log(tau) + t1
