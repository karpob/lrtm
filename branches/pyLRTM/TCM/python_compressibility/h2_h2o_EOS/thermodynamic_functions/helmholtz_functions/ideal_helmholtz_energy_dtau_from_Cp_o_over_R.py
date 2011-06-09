def ideal_helmholtz_energy_dtau_from_Cp_o_over_R(ni,ti,vi,ui,R,ho,so,dc,Tc,tau,delta,n_power,n_exp,To,do):
	from math import pow, exp, log

        c0=0
        c1=0
	
	out2=0
        #'Add up all values that contribute to the tau^0 term
	c0 = -so / R - 1 - log(do/dc*To/Tc)
	c0 = c0 + ni + ni * log(To) - ni * log(Tc)
        
	#'Add up all values that contribute to the "linear" tau term
	c1 = ho / R / Tc - ni * To/Tc
        
        term=0.0
	if(n_exp>0):
        	for i in range(n_power,n_exp):
                	beta=ui[i]
               
                	Tstar_on_tau=Tc/tau
			
			expo = exp(-beta/Tstar_on_tau)
			term += vi[i] * (beta /Tstar_on_tau)/tau* expo / (1.0-expo)
			c1 = c1 - vi[i] * ui[i] / Tc / (exp(ui[i]/To) - 1)
        
        out=(c1)+(ni/tau)+term+(-1.0/tau)

        return out  
	
