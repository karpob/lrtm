def ideal_helmholtz_from_coef_vector(delta,tau,ideal_n,ideal_gamma):
    from numpy import log,exp,meshgrid
    tau_m,ideal_n_m=meshgrid(tau,ideal_n)
    delta_m,ideal_gamma_m=meshgrid(delta,ideal_gamma)
    
    t1=ideal_n_m[3:len(ideal_gamma)][:]*log(1.0-exp(-ideal_gamma_m[3:len(ideal_gamma)][:]*tau_m[3:len(ideal_gamma)][:]))    
    t1_v=t1.sum(axis=0)
    return log(delta) + ideal_n_m[0][:] + ideal_n_m[1][:]*tau_m[2][:] + ideal_n_m[2]*log(tau_m[2][:]) + t1_v
    
  
