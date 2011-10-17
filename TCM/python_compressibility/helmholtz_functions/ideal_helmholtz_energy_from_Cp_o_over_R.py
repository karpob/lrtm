def ideal_helmholtz_energy_from_Cp_o_over_R(ni,ti,vi,ui,R,ho,so,dc,Tc,tau,delta,n_pow,n_exp,To,do):
	from math import pow, exp, log
	from pylab import zeros
	#To=273.15        
        #do=((0.001/(R*To))*1000)
        c0=0
        c1=0
	betas_ao=zeros(len(ui))
	out2=0
        #'Add up all values that contribute to the tau^0 "constant" term (the constants)
	c0 = -so / R - 1 - log(do/dc*To/Tc)
	c0 = c0 + ni + ni * log(To) - ni * log(Tc)
        
	#'Add up all values that contribute to the "linear" tau term
	c1 = ho / R / Tc - ni * To/Tc
        
        if (n_exp>0):
        	for i in range(n_pow,n_exp):
			uu=-1.0*ui[i]
               		Tstar_on_tau=Tc/tau
			out2=out2+vi[i]*log(1-exp(uu/Tstar_on_tau))
                	v=ui[i]/To
			c0 = c0 + vi[i] * v * exp(v) / (exp(v) - 1) - vi[i] * log(exp(v) - 1)
			c1 = c1 - vi[i] * ui[i] / Tc / (exp(ui[i]/To) - 1)
        out=out2+c0+c1*tau-log(tau)+log(delta)+ni*log(tau)
	#print c1
	return out

