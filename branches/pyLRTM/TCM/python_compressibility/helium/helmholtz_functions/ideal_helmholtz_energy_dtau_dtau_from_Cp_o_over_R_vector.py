def ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R_vector(ni,ti,vi,ui,R,ho,so,dc,Tc,tau,delta,n_power,n_exp,To,do):
	#from math import pow, exp, log
	from numpy import power, exp, log,meshgrid,sum

        tau_m,ui_m=meshgrid(tau,ui)
        delta_m,vi_m=meshgrid(delta,vi)
	#To=273.15        
        #do=((0.001/(R*To))*1000)
        #c0=0
        #c1=0
	
	#out2=0
        #'Add up all values that contribute to the tau^0 "constant" term (the constants)
        # cross your fingers replacing 2.5 with ni
	c0 = -so / R - 1 - log(do/dc*To/Tc)
	c0 = c0 + ni + ni * log(To) - ni * log(Tc)
        
	#'Add up all values that contribute to the "linear" tau term
	c1 = ho / R / Tc - ni * To/Tc
        
        #term=0.0
        #sumz=2.5 # 2.5 is added here as the only power term T^0
        T=Tc/tau_m[n_power:n_exp][:]
        #for i in range(n_power,n_exp):
        #        x=ui[i]/T
        #        ea=exp(-x)
        #        d=(1.0-ea)*(1.0-ea)
        #        term=vi[i]*x*x*ea/d
        #        sumz +=term
        x=ui_m[n_power:n_exp][:]/T
        ea=exp(-x)
        d=(1.0-ea)*(1.0-ea)
        term=vi_m[n_power:n_exp][:]*x*x*ea/d
        sumz=term.sum(axis=0)+ni
        
        out=(1-sumz)/power(tau_m[1][:],2)
        
        return out  
        
        

