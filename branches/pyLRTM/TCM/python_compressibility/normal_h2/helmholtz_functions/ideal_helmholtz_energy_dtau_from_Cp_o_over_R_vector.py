def ideal_helmholtz_energy_dtau_from_Cp_o_over_R_vector(ni,ti,vi,ui,R,ho,so,dc,Tc,tau,delta,n_power,n_exp,To,do):
	from numpy import power, exp, log,meshgrid,sum

        tau_m,ui_m=meshgrid(tau,ui)
        delta_m,vi_m=meshgrid(delta,vi)
        
        #'Add up all values that contribute to the tau^0 term 
        #note this won't work for anything with multiple power terrms!
        # Designed to work with H2, and He only!
        #  cross your fingers replacing 2.5 with ni
	c0 = -so / R - 1 - log(do/dc*To/Tc)
	c0 = c0 + ni + ni * log(To) - ni * log(Tc)
        
	#'Add up all values that contribute to the "linear" tau term
	c1 = ho / R / Tc - ni * To/Tc
        
        #term=0.0
	#if(n_exp>0):
        #	for i in range(n_power,n_exp):
        #        	beta=ui[i]
               
        #        	Tstar_on_tau=Tc/tau
			
	#		expo = exp(-beta/Tstar_on_tau)
	#		term += vi[i] * (beta /Tstar_on_tau)/tau* expo / (1.0-expo)
	#		c1 = c1 - vi[i] * ui[i] / Tc / (exp(ui[i]/To) - 1)
        beta=ui_m[n_power:n_exp][:]
        Tstar_on_tau=Tc/tau_m[n_power:n_exp][:]
        expo=exp(-beta/Tstar_on_tau)
        term=vi_m[n_power:n_exp][:]*(beta/Tstar_on_tau)/tau*expo/(1.0-expo)
        cminus=vi_m[n_power:n_exp][:]*ui_m[n_power:n_exp]/Tc/(exp(ui_m[n_power:n_exp]/To)-1.0)
        cminus_v=cminus.sum(axis=0)
        term_v=term.sum(axis=0)
        cterm=c1-cminus_v
        
        out=cterm+(ni/tau)+term_v+(-1.0/tau)

        return out  
	
