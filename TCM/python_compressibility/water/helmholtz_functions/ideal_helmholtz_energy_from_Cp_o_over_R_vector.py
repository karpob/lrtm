def ideal_helmholtz_energy_from_Cp_o_over_R_vector(ni,ti,vi,ui,R,ho,so,dc,Tc,tau,delta,n_pow,n_exp,To,do):
	from numpy import power, exp, log,meshgrid,sum

        tau_m,ui_m=meshgrid(tau,ui)
        delta_m,vi_m=meshgrid(delta,vi)
        
	#To=273.15        
        #do=((0.001/(R*To))*1000)
        #c0=0
        #c1=0
	#betas_ao=zeros(len(ui))
	#out2=0
        #'Add up all values that contribute to the tau^0 "constant" term (the constants)
        # replaced 2.5 with ni, cross your fingers and hope this works!
	c0 = -so / R - 1 - log(do/dc*To/Tc)
	c0 = c0 + ni + ni * log(To) - ni * log(Tc)
        
	#'Add up all values that contribute to the "linear" tau term
	c1 = ho / R / Tc - ni * To/Tc
        
        #if (n_exp>0):
        #	for i in range(n_pow,n_exp):
	#		uu=-1.0*ui[i]
        #       		Tstar_on_tau=Tc/tau
	#		out2=out2+vi[i]*log(1-exp(uu/Tstar_on_tau))
        #        	v=ui[i]/To
	#		c0 = c0 + vi[i] * v * exp(v) / (exp(v) - 1) - vi[i] * log(exp(v) - 1)
	#		c1 = c1 - vi[i] * ui[i] / Tc / (exp(ui[i]/To) - 1)
	if(n_exp>0):
	        uu=-1.0*ui_m[n_pow:n_exp][:]
	        Tstar_on_tau=Tc/tau_m[n_pow:n_exp][:]
	        out2=vi_m[n_pow:n_exp][:]*log(1-exp(uu/Tstar_on_tau))
	        out2_v=out2.sum(axis=0)
	        v=ui_m[n_pow:n_exp][:]/To
	        CC0=vi_m[n_pow:n_exp][:]*v*exp(v)/ (exp(v)-1.0) -vi_m[n_pow:n_exp][:]*log(exp(v)-1.0)
	        CC0_v=CC0.sum(axis=0)+c0
	        CC1=vi_m[n_pow:n_exp][:]*ui_m[n_pow:n_exp][:]/Tc/(exp(ui_m[n_pow:n_exp][:]/To)-1)
	        CC1_v=c1-CC1.sum(axis=0)
	        
	#print CC1_v
	out=out2_v+CC0_v+CC1_v*tau_m[1][:]-log(tau_m[1][:])+log(delta_m[1][:])+ni*log(tau_m[1][:])#+c0+c1	
        #out=out2+c0+c1*tau-log(tau)+log(delta)+2.5*log(tau)
	#vi_m[n_pow:n_exp][:]*ui_m[n_pow:n_exp][:]#*exp(ui_m[n_pow:n_exp][:])
	return out

