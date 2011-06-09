def ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R(ni,ti,vi,ui,R,ho,so,dc,Tc,tau,delta,n_power,n_exp,To,do):
	from math import pow, exp, log
	#To=273.15        
        #do=((0.001/(R*To))*1000)
        c0=0
        c1=0
	
	out2=0
        #'Add up all values that contribute to the tau^0 "constant" term (the constants)
	c0 = -so / R - 1 - log(do/dc*To/Tc)
	c0 = c0 + ni + ni * log(To) - ni * log(Tc)
        
	#'Add up all values that contribute to the "linear" tau term
	c1 = ho / R / Tc - ni * To/Tc
        
        term=0.0
        sumz=ni # 2.5 is added here as the only power term T^0
        T=Tc/tau
        for i in range(n_power,n_exp):
                x=ui[i]/T
                ea=exp(-x)
                d=(1.0-ea)*(1.0-ea)
                term=vi[i]*x*x*ea/d
                sumz +=term
        
        out=(1-sumz)/pow(tau,2)
        
        return out  
        
        

