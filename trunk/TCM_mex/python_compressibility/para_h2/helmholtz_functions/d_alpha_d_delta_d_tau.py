def d_alpha_d_delta_d_tau(tau,delta,RES_n,RES_t,RES_d,RES_c,RES_alpha,RES_beta,RES_gamma,RES_eps,n_power_terms,n_gaussian_terms,n_critical_terms,RES_a,RES_b,RES_B,RES_C,RES_D,RES_A):
        from math import exp, pow

        tt1=0.0
        tt2=0.0
        tt3=0.0
        tt4=0.0
        for i in range(0,n_power_terms):
                if(abs(float(RES_c[i]))>0.0):
                        p1  = (RES_d[i] - RES_c[i]*pow(delta,RES_c[i]))*exp((-pow(delta,RES_c[i])))
                        tt2 += RES_n[i]*RES_t[i]*pow(delta,RES_d[i] - 1.0)*pow(tau,RES_t[i] - 1.0)*p1
                else:
                        tt1 += RES_n[i]*RES_d[i]*RES_t[i]*pow(delta,(RES_d[i] - 1.0))*pow(tau,(RES_t[i] - 1.0))
        
        for i in range(n_power_terms,n_power_terms+n_gaussian_terms):
                d1 = delta - RES_eps[i]
		t1 = tau - RES_gamma[i]
		e1 = -1.0*RES_alpha[i]*d1*d1 - RES_beta[i]*t1*t1
                f1 = RES_t[i] - 2*RES_beta[i]*tau*(tau - RES_gamma[i])
		g1 = RES_d[i] - 2*RES_alpha[i]*delta*(delta - RES_eps[i])

		sumz = RES_n[i] * f1 * pow(tau,RES_t[i]-1) * g1 * pow(delta,RES_d[i]-1) * exp(e1)
		tt3 += sumz
                
     
        for i in range(n_power_terms+n_gaussian_terms,n_power_terms+n_gaussian_terms+n_critical_terms):
                d1 = delta - 1.0
                t1 = tau - 1.0
                psi = exp(-RES_C[i]*d1*d1 - RES_D[i]*t1*t1)
                theta = 1 - tau + RES_A[i]*pow(d1*d1, (1/(2*RES_beta[i])))
                Delta = theta*theta + RES_B[i]*pow(d1*d1, (RES_a[i]))
                dDeltabidtau = -2*theta*RES_b[i]*pow(Delta, (RES_b[i]-1))
                dpsidtau = -2*RES_D[i]*t1*psi
                dDeltaddelta = d1*(RES_A[i]*theta*2/RES_beta[i]*pow(d1*d1, (1/(2*RES_beta[i]) - 1)) + 2*RES_B[i]*RES_a[i]*pow(d1*d1, (RES_a[i] - 1)))
                d2psiddeltadtau = 4*RES_C[i]*RES_D[i]*d1*t1*psi
                dDeltabiddelta = RES_b[i]*pow(Delta, (RES_b[i]-1))*dDeltaddelta
                dpsiddelta = -2*RES_C[i]*d1*psi
                d2Deltabiddeltadtau = -RES_A[i]*RES_b[i]*2/RES_beta[i]*pow(Delta, (RES_b[i]-1))*d1*pow(d1*d1, (1/(2*RES_beta[i])-1)) -2*theta*RES_b[i]*(RES_b[i]-1)*pow(Delta, (RES_b[i]-2))*dDeltaddelta
                tt4 += RES_n[i]*(pow(Delta, RES_b[i])*(dpsidtau + delta*d2psiddeltadtau) + delta*dDeltabiddelta*dpsidtau +dDeltabidtau*(psi + delta*dpsiddelta) +d2Deltabiddeltadtau*delta*psi)
                        
    
        return tt1+tt2+tt3+tt4
