def enthalpy_helium(T,density):
        from math import sqrt,exp
        from numpy import sort,asarray,shape,concatenate,power
	from helmholtz_functions.ideal_helmholtz_energy_from_Cp_o_over_R import ideal_helmholtz_energy_from_Cp_o_over_R
	from helmholtz_functions.ideal_helmholtz_energy_dtau_from_Cp_o_over_R import ideal_helmholtz_energy_dtau_from_Cp_o_over_R
	from helmholtz_functions.ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R import ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R	
	from helmholtz_functions.mBWR_to_Helmholtz import mBWR_to_Helmholtz
	from helmholtz_functions.helmholtz_energy_residual import  helmholtz_energy_residual
	from helmholtz_functions.d_alpha_d_tau import d_alpha_d_tau
        from helmholtz_functions.d_alpha_d_tau_d_tau import d_alpha_d_tau_d_tau
	from helmholtz_functions.d_alpha_d_delta import d_alpha_d_delta
	from helmholtz_functions.d_alpha_d_delta_d_delta import d_alpha_d_delta_d_delta
	from helmholtz_functions.d_alpha_d_delta_d_tau import d_alpha_d_delta_d_tau

	M_He=4.0026#02 g/mol
        R=8.314310 #

        R_He=R/M_He
        To=4.230359714841141  #Kelvin
	rho_o=31.163394763964778 # mol/L  
	ho=108.78863197310453    # J/mol
        so=3.6929233790579463    # J/mol/K
	
        
        Tc=5.19530# Kelvin
        rho_c=17.3990 #(mol/L)
                
        mBWR_Coef=[0.4558980227431e-4,#1 
                   0.1260692007853e-2,#2
                   -0.7139657549318e-2,#3
                   0.9728903861441e-2,#4
                   -0.1589302471562e-1,#5
                   0.1454229259623e-5,#6
                   -0.4708238429298e-4,#7
                   0.1132915232587e-2,#8
                   0.2410763742104e-2,#9
                   -0.5093547838381e-8,#10
                   0.2699726927900e-5,#11
                   -0.3954146691114e-4,#12`
                   0.1551961438127e-8,#13
                   0.1050712335785e-7,#14
                   -0.5501158366750e-7,#15
                   -0.1037673478521e-9,#16
                   0.6446881346448e-12,#17
                   0.3298960057071e-10,#18
                   -0.3555585738784e-12,#19
                   -0.6885401367690e-2,#20
                   0.9166109232806e-2,#21
                   -0.6544314242937e-5,#22
                   -0.3315398880031e-4,#23
                   -0.2067693644676e-7,#24
                   0.3850153114958e-7,#25
                   -0.1399040626999e-10,#26
                   -0.1888462892389e-11,#27
                   -0.4595138561035e-14,#28
                   0.6872567403738e-14,#29
                   -0.6097223119177e-18,#30
                   -0.7636186157005e-17,#31
                   0.3848665703556e-17]#32
        rho=density/M_He    #mol/L       
#       delta=rho/rho_c      
#       tau=Tc/T
#	tau_o=Tc/To
#	delta_o=rho_o/rho_c
        
	ni=2.5
	ti=0.0
	vi=[0.0]
	ui=[0.0]
	npower=1
	n_exp=0
	tau=Tc/T
	delta=rho/rho_c
	n_power=1
        #phi_o=ideal_helmholtz_energy_from_Cp_o_over_R(ni,ti,vi,ui,R,ho,so,rho_c,Tc,tau,delta,n_power,n_exp,To,rho_o)

        #N=mBWR_Coef
        #abwr=ABWR(N,T,rho,rho_c)
	#residual=abwr/(R*T)
	[N_i,d_i,t_i,p_i]=mBWR_to_Helmholtz(mBWR_Coef,T,rho,rho_c,Tc,R)
	
	
	#print shape(sorted_array)
	#print d_i
	#print t_i
	#N_i=sorted_array[3,:]
	#d_i=sorted_array[2,:]
	#t_i=sorted_array[1,:]
	#p_i=sorted_array[0,:]
	#print 'heres p_i',p_i,p_ii
	Beta_i=0
	gamma_i=0
	D_i=0
	n_power_terms=len(N_i)
	n_gaussian_terms=0
	n_critical_terms=0
	RES_a=0
	RES_b=0
	RES_B=0
	RES_C=0
	RES_D=0
	RES_A=0
	phi_i=0
        #residual2=helmholtz_energy_residual(tau,delta,ni,ti,di,pi,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	#helmholtz=residual2#-residual#residual2#(ideal+residual)/M_He#(ideal+residual)/M_He
        #print residual,residual2

	#Calculate Ideal Terms
	ideal=ideal_helmholtz_energy_from_Cp_o_over_R(ni,ti,vi,ui,R,ho,so,rho_c,Tc,tau,delta,n_power,n_exp,To,rho_o)
        
        dalpha_o=ideal_helmholtz_energy_dtau_from_Cp_o_over_R(ni,ti,vi,ui,R,ho,so,rho_c,Tc,tau,delta,n_power,n_exp,To,rho_o)
	
	dalpha_o_tau_tau=ideal_helmholtz_energy_dtau_dtau_from_Cp_o_over_R(ni,ti,vi,ui,R,ho,so,rho_c,Tc,tau,delta,n_power,n_exp,To,rho_o)
	
	#Calculate Residual Terms
	residual=helmholtz_energy_residual(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	
	dalpha_tau=d_alpha_d_tau(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	
	dalpha_delta=d_alpha_d_delta(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
	
	dalpha_tau_tau=d_alpha_d_tau_d_tau(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
        
        dalpha_delta_delta=d_alpha_d_delta_d_delta(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)
        
        dalpha_delta_tau=d_alpha_d_delta_d_tau(tau,delta,N_i,t_i,d_i,p_i,phi_i,Beta_i,gamma_i,D_i,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A)

	#compressibility Density<-->Pressure
	Z=delta*dalpha_delta+float(1.0)
        #Enthalpy
	h1=R_He*T*1.0
	h2=R_He*T*tau*(dalpha_o+dalpha_tau)
	h3=R_He*T*delta*dalpha_delta
	h=h1+h2+h3
	#Entropy
	s=R_He*(tau*dalpha_o+tau*dalpha_tau-ideal-residual)
	#Absolute Helmholtz Energy
	absolute_helmholtz_energy=R_He*T*(ideal+residual)
	#Isochoric heat capacity
	cv=-1.0*R_He*(tau*tau*(dalpha_o_tau_tau+dalpha_tau_tau))
	#Isobaric heat capacity
	cp=cv+R_He*pow((1.0+delta*dalpha_delta-delta*tau*dalpha_delta_tau),2)/(1.0+2.0*delta*dalpha_delta+delta*delta*dalpha_delta_delta)
	
	# Speed of Sound?
	
	a1 = 1 + delta*dalpha_delta - delta*tau*dalpha_delta_tau;
        b1 = tau*tau*(dalpha_o_tau_tau + dalpha_tau_tau);
        w = 1 + 2*delta*dalpha_delta + delta*delta*dalpha_delta_delta - a1*a1/b1;
        speed_o_sound=power(R_He*T*w*1000,0.5)


        return Z,h,s,absolute_helmholtz_energy,cv,cp,speed_o_sound,ideal,residual
