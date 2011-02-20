def ideal_helmholtz_from_coef_dtau_vector(delta,tau,ideal_n,ideal_gamma):
    from numpy import power,exp,meshgrid
    tau_m,ideal_n_m=meshgrid(tau,ideal_n)
    delta_m,ideal_gamma_m=meshgrid(delta,ideal_gamma)
    #for i in range(3,len(ideal_gamma)):
    #    t1 += ideal_n[i]*ideal_gamma[i]*(pow(1.0 - exp(-ideal_gamma[i]*tau),-1.0) - 1.0)
    #return ideal_n[1] + ideal_n[2]/tau + t1
    t1=ideal_n_m[3:len(ideal_n)][:]*ideal_gamma_m[3:len(ideal_gamma)][:]*(power(1.0-exp(-ideal_gamma_m[3:len(ideal_gamma)][:]*tau_m[3:len(ideal_gamma)][:]),-1.0)-1.0)
    t1_v=t1.sum(axis=0) 
    return ideal_n_m[1][:]+ideal_n_m[2][:]/tau_m[2][:]+t1_v
