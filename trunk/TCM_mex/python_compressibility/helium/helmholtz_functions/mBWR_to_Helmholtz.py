def mBWR_to_Helmholtz(b,t,rho,rho_c,T_c,R):
        from numpy import zeros,multiply
        from math import sqrt,exp

	thou_over_R=1000.0/(R) #factor to normalize helmholtz energy, extra T^-1 is added to temperature powers

	temperature_powers=zeros(32)
	ti=zeros(32)
	temperature_powers=[float(1.0),float(0.5),float(0.0),float(-1),float(-2),   #2
			float(1),float(0.0),float(-1),float(-2),           #3
			float(1),float(0.0),float(-1),              #4
                        float(0.0),                   #5
			float(-1),float(-2),               #6
               		float(-1),                  #7
			float(-1),float(-2),               #8
			float(-2),                  #9
			float(-2),float(-3),               #10
			float(-2),float(-4),               #11
			float(-2),float(-3),               #12
			float(-2),float(-4),               #13
			float(-2),float(-3),               #14
			float(-2),float(-3),float(-4)]            #15

        temperature_powers=multiply(1.0,temperature_powers)
	temperature_powers=temperature_powers-1.0  # to account for the extra T^-1 to 
                                                 #  normalize helmholtz energy (thou_over_R takes care of the rest)
	ti=multiply(-1.0,temperature_powers)
	ti[0]=0.0
	ti[5]=0.0
	ti[9]=0.0
	density_powers=[float(1.0),float(1.0),float(1.0),float(1.0),float(1.0),   #2
			float(2),float(2),float(2),float(2),           #3
			float(3),float(3),float(3),              #4
                        float(4),                   #5
			float(5),float(5),               #6
              		float(6),                  #7
			float(7),float(7),               #8
			float(8),                  #9
			float(2),float(2),               #10
			float(4),float(4),               #11
			float(6),float(6),               #12
			float(8),float(8),               #13
			float(10),float(10),               #14
			float(12),float(12),float(12)]            #15
	di=density_powers           # Note, these change  with A10-15 these are best 
                                    # interpreted as powers for cricitcal density
				    # Powers after #9 are sometimes replaced with various "extra" polynomial terms, "extra" delta terms, etc 
	# Bloody Hell, this gets big 80 Coefficients from an initial 32! Thanks Bennedict, Webb, and Ruben!

	total_length=80

	ni=zeros(total_length)
	di_out=zeros(total_length)
	ti_out=zeros(total_length)
	pi=zeros(total_length)
	ni_out=zeros(total_length)
	
	for i in range(0, len(ti)):
		ni[i]=b[i]*pow(rho_c,di[i])*pow(T_c,temperature_powers[i])   #Note, rho_c^density_power is taken care of here
		
	for i in range(0,19):
		ni_out[i]=thou_over_R*ni[i]/(di[i]) # see McLinden B6, or Span
		pi[i]=0.0
		di_out[i]=di[i]
		ti_out[i]=ti[i]

#############################################
#    Calculate values associated with term A10 
#############################################
#        ar=ar-0.5*an[10]*rhocn*(expdel-1.0)	
        ni_out[19]=(-0.5*ni[19])*thou_over_R
	ni_out[20]=(-0.5*ni[20])*thou_over_R
	ti_out[19]=ti[19]
	ti_out[20]=ti[20]
	di_out[19]=0
	di_out[20]=0
	pi[19]=2.0
	pi[20]=2.0

	# zero density term extra
	ni_out[21]=thou_over_R*(0.5*ni[19])
	ni_out[22]=thou_over_R*(0.5*ni[20])
	ti_out[21]=ti[19]
	ti_out[22]=ti[20]
	pi[21]=0
	pi[22]=0
	di_out[21]=0
	di_out[22]=0


#########################################
# Calculate terms associated with A11-A15
#########################################

#     polynomial coefficients from McLinden B6 used to avoid recursion.
	A11_vector=[-1,1,1]
	A12_vector=[-2,2,2,1]
	A13_vector=[-6,6,6,3,1]
	A14_vector=[-24,24,24,12,4,1]
	A15_vector=[-120,120,120,60,20,5,1]

	pi_vector=[0.0, 2.0, 2.0, 2.0,2.0,2.0,2.0]
	di_vector=[0.0, 0.0,2.0,4.0,6.0,8.0,10.0]

	Coef_A11_1=multiply(-0.5*thou_over_R*ni[21],A11_vector)
	Coef_A11_2=multiply(-0.5*thou_over_R*ni[22],A11_vector)

	Coef_A12_1=multiply(-0.5*thou_over_R*ni[23],A12_vector)
	Coef_A12_2=multiply(-0.5*thou_over_R*ni[24],A12_vector)

	Coef_A13_1=multiply(-0.5*thou_over_R*ni[25],A13_vector)
	Coef_A13_2=multiply(-0.5*thou_over_R*ni[26],A13_vector)

	Coef_A14_1=multiply(-0.5*thou_over_R*ni[27],A14_vector)
	Coef_A14_2=multiply(-0.5*thou_over_R*ni[28],A14_vector)

	Coef_A15_1=multiply(-0.5*thou_over_R*ni[29],A15_vector)
	Coef_A15_2=multiply(-0.5*thou_over_R*ni[30],A15_vector)
	Coef_A15_3=multiply(-0.5*thou_over_R*ni[31],A15_vector)

#Still kinda Kludgey, but you can read it, and RECURSION.EQ.EVIL
	for i in range(11,16):
		######################################
		# Calculate terms Associated with A11
		######################################	
		if(i==11):
			for j in range(0,len(Coef_A11_1)):				
				ni_out[j+23]=Coef_A11_1[j]
				ti_out[j+23]=ti[21]
				pi[j+23]=pi_vector[j]
				di_out[j+23]=di_vector[j]
			current_length=len(Coef_A11_1)+23
			ti_offset=21+1
			for j in range(0,len(Coef_A11_2)):
				ni_out[j+current_length]=Coef_A11_2[j]
				ti_out[j+current_length]=ti[ti_offset]
				pi[j+current_length]=pi_vector[j]
				di_out[j+current_length]=di_vector[j]
			current_length=len(Coef_A11_2)+current_length
			ti_offset=ti_offset+1
		######################################
		# Calculate terms Associated with A12
		######################################
		if(i==12):
			for j in range(0,len(Coef_A12_1)):
				ni_out[j+current_length]=Coef_A12_1[j]
				ti_out[j+current_length]=ti[ti_offset]
				pi[j+current_length]=pi_vector[j]
				di_out[j+current_length]=di_vector[j]
			current_length=len(Coef_A12_1)+current_length
			ti_offset=ti_offset+1
			for j in range(0,len(Coef_A12_2)):
				ni_out[j+current_length]=Coef_A12_2[j]
				ti_out[j+current_length]=ti[ti_offset]
				pi[j+current_length]=pi_vector[j]
				di_out[j+current_length]=di_vector[j]
			current_length=len(Coef_A12_2)+current_length
			ti_offset=ti_offset+1
		######################################
		# Calculate terms Associated with A13
		######################################	
		if(i==13):
			for j in range(0,len(Coef_A13_1)):
				ni_out[j+current_length]=Coef_A13_1[j]
				ti_out[j+current_length]=ti[ti_offset]
				pi[j+current_length]=pi_vector[j]
				di_out[j+current_length]=di_vector[j]
			current_length=len(Coef_A13_1)+current_length
			ti_offset=ti_offset+1
			for j in range(0,len(Coef_A13_2)):
				ni_out[j+current_length]=Coef_A13_2[j]
				ti_out[j+current_length]=ti[ti_offset]
				pi[j+current_length]=pi_vector[j]
				di_out[j+current_length]=di_vector[j]
			current_length=len(Coef_A13_2)+current_length
			ti_offset=ti_offset+1		
		######################################
		# Calculate terms Associated with A14
		######################################
		if(i==14):
			for j in range(0,len(Coef_A14_1)):
				ni_out[j+current_length]=Coef_A14_1[j]
				ti_out[j+current_length]=ti[ti_offset]
				pi[j+current_length]=pi_vector[j]
				di_out[j+current_length]=di_vector[j]
			current_length=len(Coef_A14_1)+current_length
			ti_offset=ti_offset+1
			for j in range(0,len(Coef_A14_2)):
				ni_out[j+current_length]=Coef_A14_2[j]
				ti_out[j+current_length]=ti[ti_offset]
				pi[j+current_length]=pi_vector[j]
				di_out[j+current_length]=di_vector[j]
			current_length=len(Coef_A14_2)+current_length
			ti_offset=ti_offset+1
		######################################
		# Calculate terms Associated with A15
		#######################################		
		if(i==15):
			for j in range(0,len(Coef_A15_1)):
				ni_out[j+current_length]=Coef_A15_1[j]
				ti_out[j+current_length]=ti[ti_offset]
				pi[j+current_length]=pi_vector[j]
				di_out[j+current_length]=di_vector[j]
			current_length=len(Coef_A15_1)+current_length
			ti_offset=ti_offset+1
			for j in range(0,len(Coef_A15_2)):
				ni_out[j+current_length]=Coef_A15_2[j]
				ti_out[j+current_length]=ti[ti_offset]
				pi[j+current_length]=pi_vector[j]
				di_out[j+current_length]=di_vector[j]
			current_length=len(Coef_A15_2)+current_length
			ti_offset=ti_offset+1		
			for j in range(0,len(Coef_A15_3)):
				ni_out[j+current_length]=Coef_A15_3[j]
				ti_out[j+current_length]=ti[ti_offset]
				pi[j+current_length]=pi_vector[j]
				di_out[j+current_length]=di_vector[j]
	sorted_index=pi.argsort(axis=0)
	ni_out=ni_out[sorted_index]
	di_out=di_out[sorted_index]
	ti_out=ti_out[sorted_index]
	pi=pi[sorted_index]
	count_zero=0
	count_two=0
	                      		
        return ni_out,di_out,ti_out,pi
