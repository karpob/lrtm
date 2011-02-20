def d_alpha_o_d_tau(tau,c0,c1,c2,betas_ao,vi):
	from math import pow, exp, log
	from pylab import zeros
	out2=0.0
	out1=c1+c2/tau
	
	for i in range(0,len(betas_ao)):
		out2=out2-vi[i]*(betas_ao[i]/tau)*(pow(1.0-exp(betas_ao[i]),-1.0)-1.0)
        out=out1+out2
	return out
	
#	def IG_phi_tau(delta,tau):
#    t1 = 0.0
#    for i in range(4,9):
#        t1 += IG_n[i]*IG_gamma[i]*(pow(1.0 - exp(-IG_gamma[i]*tau),-1.0) - 1.0)
#    return IG_n[2] + IG_n[3]/tau + t1
