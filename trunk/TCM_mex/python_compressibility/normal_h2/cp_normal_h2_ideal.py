def cp_normal_h2_ideal(T):
	from math import exp, pow,log
	ni=2.5   
	ti=0.0
        vi=[1.616,-0.4117,-0.792,0.758,1.217] 
        ui=[531,751,1989,2484,6859]
	term1=0
	term2=0

	#for i in range(0,len(ti)):
	term1=ni*pow(T,ti)
	for i in range(0,len(ui)):
		term2=term2+vi[i]*pow((ui[i]/T),2)*(exp(-ui[i]/T))/pow(1-exp(-ui[i]/T),2)
        return term1+term2
