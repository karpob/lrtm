#!/usr/bin/env python

# Copyright (c) 2008, Kiran Pashikanti 
# 
# Permission to use, copy, modify, and/or distribute this software for any 
# purpose with or without fee is hereby granted, provided that the above 
# copyright notice and this permission notice appear in all copies. 
# 
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES 
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF 
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY 
# SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES 
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN 
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR 
# IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. 
# 

"""
IAPWS95.py

Implementation of IAPWS-95 specification for water and steam properties 
"""

from math import exp,log,pow,sqrt

CRITICAL_TEMPERATURE  = 647.096     #K
CRITICAL_DENSITY      = 322.0       #kg m^-3
SPECIFIC_GAS_CONSTANT = 0.46151805  #kJ kg^-1 K^-1

IG_n = [0,-8.32044648201,6.6832105268,3.00632,0.012436,0.97315,1.27950,0.96956,0.24873]
IG_gamma = [0,0,0,0,1.28728967,3.53734222,7.74073708,9.24437796,27.5075105]

def IG_phi(delta,tau):
    t1 = 0.0
    for i in range(4,9):
        t1 += IG_n[i]*log(1.0 - exp(-IG_gamma[i]*tau))
        
    return log(delta) + IG_n[1] + IG_n[2]*tau + IG_n[3]*log(tau) + t1
    

def IG_phi_delta(delta,tau):
    return 1.0/delta

def IG_phi_delta_delta(delta,tau):
    return -1.0/(delta*delta)

def IG_phi_tau(delta,tau):
    t1 = 0.0
    for i in range(4,9):
        t1 += IG_n[i]*IG_gamma[i]*(pow(1.0 - exp(-IG_gamma[i]*tau),-1.0) - 1.0)
    return IG_n[2] + IG_n[3]/tau + t1

def IG_phi_tau_tau(delta,tau):
    t1 = 0.0
    for i in range(4,9):
        t1 += IG_n[i]*(IG_gamma[i]**2)*exp(-IG_gamma[i]*tau)*pow(1.0 - exp(-IG_gamma[i]*tau),-2.0)
    return -IG_n[3]/(tau*tau) - t1

def IG_phi_delta_tau(delta,tau):
    return 0.0

RES_n =  [0,
          0.12533547935523e-1,  0.78957634722828e+1, -0.87803203303561e+1,
          0.31802509345418e+0, -0.26145533859358e+0, -0.78199751687981e-2,
          0.88089493102134e-2, -0.66856572307965e+0,  0.20433810950965e+0,
         -0.66212605039687e-4, -0.19232721156002e+0, -0.25709043003438e+0,
          0.16074868486251e+0, -0.40092828925807e-1,  0.39343422603254e-6,
         -0.75941377088144e-5,  0.56250979351888e-3, -0.15608652257135e-4,
          0.11537996422951e-8,  0.36582165144204e-6, -0.13251180074668e-11,
         -0.62639586912454e-9, -0.10793600908932e+0,  0.17611491008752e-1,
          0.22132295167546e+0, -0.40247669763528e+0,  0.58083399985759e+0,
          0.49969146990806e-2, -0.31358700712549e-1, -0.74315929710341e+0,
          0.47807329915480e+0,  0.20527940895948e-1, -0.13636435110343e+0,
          0.14180634400617e-1,  0.83326504880713e-2, -0.29052336009585e-1,
          0.38615085574206e-1, -0.20393486513704e-1, -0.16554050063734e-2,
          0.19955571979541e-2,  0.15870308324157e-3, -0.16388568342530e-4,
          0.43613615723811e-1,  0.34994005463765e-1, -0.76788197844621e-1,
          0.22446277332006e-1, -0.62689710414685e-4, -0.55711118565645e-9,
         -0.19905718354408e+0,  0.31777497330738e+0, -0.11841182425981e+0,
         -0.31306260323435e+2,  0.31546140237781e+2, -0.25213154341695e+4,
         -0.14874640856724e+0,  0.31806110878444e+0]

RES_c =  [0,
          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
          1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0,
          2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
          2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 3.0, 4.0, 6.0, 6.0, 6.0, 6.0]
RES_d =  [0,
          1.0, 1.0, 1.0, 2.0,  2.0,  3.0,  4.0,  1.0,  1.0, 1.0, 2.0,  2.0,  3.0,  4.0,
          4.0, 5.0, 7.0, 9.0, 10.0, 11.0, 13.0, 15.0,  1.0, 2.0, 2.0,  2.0,  3.0,  4.0,
          4.0, 4.0, 5.0, 6.0,  6.0,  7.0,  9.0,  9.0,  9.0, 9.0, 9.0, 10.0, 10.0, 12.0,
          3.0, 4.0, 4.0, 5.0, 14.0,  3.0,  6.0,  6.0,  6.0, 3.0, 3.0,  3.0]
RES_t =  [0,
          -0.5, 0.875,  1.0,  0.5,  0.75, 0.375,  1.0,  4.0,  6.0, 12.0,  1.0,
           5.0, 4.0  ,  2.0, 13.0,  9.0 , 3.0  ,  4.0, 11.0,  4.0, 13.0,  1.0,
           7.0, 1.0  ,  9.0, 10.0, 10.0 , 3.0  ,  7.0, 10.0, 10.0,  6.0, 10.0,
          10.0, 1.0  ,  2.0,  3.0,  4.0 , 8.0  ,  6.0,  9.0,  8.0, 16.0, 22.0,
          23.0,23.0  , 10.0, 50.0, 44.0, 46.0  , 50.0,  0.0,  1.0,  4.0]

#offset constants, no point in creating giant nearly empty lists
RES_alpha = {52:20.0,53:20.0,54:20.0}
RES_beta  = {52:150.0,53:150.0,54:250.0,55:0.3,56:0.3}
RES_gamma = {52:1.21,53:1.21,54:1.25}
RES_eps   = {52:1.0, 53:1.0 ,54:1.0 }
RES_a     = {55:3.5 ,56:3.5}
RES_b     = {55:0.85,56:0.95}
RES_B     = {55:0.2,56:0.2}
RES_C     = {55:28.0,56:32.0}
RES_D     = {55:700.0,56:800.0}
RES_A     = {55:0.32,56:0.32}

def RES_phi(delta,tau):
    t1 = 0.0
    for i in range(1,8):
        t1 += RES_n[i]*(pow(delta,RES_d[i]))*(pow(tau,RES_t[i]))

    t2 = 0.0
    for i in range(8,52):
        t2 += RES_n[i]*(pow(delta,RES_d[i]))*(pow(tau,RES_t[i]))*exp(-pow(delta,RES_c[i]))

    t3 = 0.0
    for i in range(52,55):
        p1  = (delta-RES_eps[i])*(delta-RES_eps[i])
        p2  = (tau-RES_gamma[i])*(tau-RES_gamma[i])
        t3 += RES_n[i]*(pow(delta,RES_d[i]))*(pow(tau,RES_t[i]))*exp(-RES_alpha[i]*p1 - RES_beta[i]*p2)
        
    t4 = 0.0
    for i in range(55,57):
        THETA = (1.0 - tau) + RES_A[i]*pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i])
        DELTA = THETA*THETA + RES_B[i]*pow((delta - 1.0)*(delta - 1.0),RES_a[i])
        PSI   = exp(-RES_C[i]*((delta - 1.0)*(delta - 1.0)) - RES_D[i]*((tau - 1.0)*(tau - 1.0)))
        t4 += RES_n[i]*pow(DELTA,RES_b[i])*delta*PSI

    return t1 + t2 + t3 + t4

def RES_phi_delta(delta,tau):
    t1 = 0.0
    for i in range(1,8):
        t1 += RES_n[i]*RES_d[i]*pow(delta,RES_d[i]-1.0)*pow(tau,RES_t[i])
    
    t2 = 0.0
    for i in range(8,52):
        t2 += RES_n[i]*exp(-pow(delta,RES_c[i]))*(pow(delta,RES_d[i]-1.0)*pow(tau,RES_t[i])*(RES_d[i] - RES_c[i]*pow(delta,RES_c[i])))

    t3 = 0.0
    for i in range(52,55):
        p1 = exp(-RES_alpha[i]*((delta - RES_eps[i])*(delta - RES_eps[i])) - RES_beta[i]*((tau - RES_gamma[i])*(tau - RES_gamma[i])))
        p2 = (RES_d[i]/delta) - 2.0*RES_alpha[i]*(delta - RES_eps[i])
        t3 += RES_n[i]*pow(delta,RES_d[i])*pow(tau,RES_t[i])*p1*p2
    
    t4 = 0.0
    for i in range(55,57):
        THETA = (1.0 - tau) + RES_A[i]*pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i])
        DELTA = THETA*THETA + RES_B[i]*pow((delta - 1.0)*(delta - 1.0),RES_a[i])
        PSI   = exp(-RES_C[i]*((delta - 1.0)*(delta - 1.0)) - RES_D[i]*((tau - 1.0)*(tau - 1.0)))

        D_DELTA_delta = (delta - 1.0)*(RES_A[i]*THETA*(2.0/RES_beta[i])*pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i] - 1.0) + 2.0*RES_B[i]*RES_a[i]*pow((delta - 1.0)*(delta - 1.0),RES_a[i]-1.0))
        D_DELTA_bi_delta = RES_b[i]*pow(DELTA,RES_b[i]-1.0)*D_DELTA_delta
        D_PSI_delta = -2.0*RES_C[i]*(delta - 1.0)*PSI
        
        t4 += RES_n[i]*(pow(DELTA,RES_b[i])*(PSI + delta*D_PSI_delta) + D_DELTA_bi_delta*delta*PSI)
        
    return t1 + t2 + t3 + t4

def RES_phi_delta_delta(delta,tau):
    t1 = 0.0
    for i in range(1,8):
        t1 += RES_n[i]*RES_d[i]*(RES_d[i]-1.0)*pow(delta,RES_d[i]-2.0)*pow(tau,RES_t[i])
    
    t2 = 0.0
    for i in range(8,52):
        p1 = pow(delta,RES_d[i] - 2.0)*pow(tau,RES_t[i])
        p2 = (RES_d[i] - RES_c[i]*pow(delta,RES_c[i]))*(RES_d[i] - 1.0 - (RES_c[i])*pow(delta,RES_c[i]))
        t2 += RES_n[i]*exp(-pow(delta,RES_c[i]))*(p1*(p2 - (RES_c[i]*RES_c[i])*pow(delta,RES_c[i])))
        
    t3 = 0.0
    for i in range(52,55):
        p1 = exp(-RES_alpha[i]*((delta - RES_eps[i])*(delta - RES_eps[i])) - RES_beta[i]*((tau-RES_gamma[i])*(tau-RES_gamma[i])))
        p2 = 2.0*RES_alpha[i]*pow(delta,RES_d[i])
        p3 = 4.0*(RES_alpha[i]*RES_alpha[i])*pow(delta,RES_d[i])*((delta - RES_eps[i])*(delta - RES_eps[i]))
        p4 = 4.0*RES_d[i]*RES_alpha[i]*pow(delta,RES_d[i] - 1.0)*(delta - RES_eps[i])
        p5 = RES_d[i]*(RES_d[i] - 1.0)*pow(delta,RES_d[i] - 2.0)
        t3 += RES_n[i]*pow(tau,RES_t[i])*p1*(-p2 + p3 - p4 + p5)

    t4 = 0.0
    for i in range(55,57):
        THETA = (1.0 - tau) + RES_A[i]*pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i])
        DELTA = THETA*THETA + RES_B[i]*pow((delta - 1.0)*(delta - 1.0),RES_a[i])
        PSI   = exp(-RES_C[i]*((delta - 1.0)*(delta - 1.0)) - RES_D[i]*((tau - 1.0)*(tau - 1.0)))

        D_DELTA_delta = (delta - 1.0)*(RES_A[i]*THETA*(2.0/RES_beta[i])*pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i] - 1.0) + 2.0*RES_B[i]*RES_a[i]*pow((delta - 1.0)*(delta - 1.0),RES_a[i]-1.0))
        D_DELTA_bi_delta = RES_b[i]*pow(DELTA,RES_b[i]-1.0)*D_DELTA_delta
        D_PSI_delta = -2.0*RES_C[i]*(delta - 1.0)*PSI
        
        D2_PSI_delta = 2.0*RES_C[i]*PSI*(2.0*RES_C[i]*(delta - 1.0)*(delta - 1.0) - 1.0)

        p1 = 4.0*RES_B[i]*RES_a[i]*(RES_a[i] - 1.0)*pow((delta - 1.0)*(delta - 1.0),RES_a[i] - 2.0)
        p2 = 2.0*(RES_A[i]*RES_A[i])*pow(RES_beta[i],-2.0)*pow(pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i] - 1.0),2)
        p3 = RES_A[i]*THETA*(4.0/RES_beta[i])*(0.5/RES_beta[i] - 1.0)*pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i] - 2.0)

        D2_DELTA_delta = (D_DELTA_delta/(delta - 1.0)) + ((delta - 1.0)*(delta - 1.0))*(p1 + p2 + p3)
        D2_DELTA_bi_delta = RES_b[i]*(pow(DELTA,RES_b[i] - 1.0)*D2_DELTA_delta + (RES_b[i] - 1.0)*pow(DELTA,RES_b[i] - 2.0)*pow(D_DELTA_delta,2.0))
        
        #p1 - p5 can be discarded now
        p1 = pow(DELTA,RES_b[i])*(2.0*D_PSI_delta + delta*D2_PSI_delta)
        p2 = 2.0*(D_DELTA_bi_delta)*(PSI + delta*D_PSI_delta)
        p3 = D2_DELTA_bi_delta*delta*PSI
        
        #finally
        t4 += RES_n[i]*(p1 + p2 + p3)
        
        
    return t1 + t2 + t3 + t4


def RES_phi_tau(delta,tau):
    t1 = 0.0
    for i in range(1,8):
        t1 += RES_n[i]*RES_t[i]*pow(delta,RES_d[i])*pow(tau,RES_t[i] - 1.0)

    t2 = 0.0
    for i in range(8,52):
        t2 += RES_n[i]*RES_t[i]*pow(delta,RES_d[i])*pow(tau,RES_t[i] - 1.0)*exp(-pow(delta,RES_c[i]))

    t3 = 0.0
    for i in range(52,55):
        p1 = -RES_alpha[i]*((delta - RES_eps[i])*(delta - RES_eps[i])) - RES_beta[i]*((tau - RES_gamma[i])*(tau - RES_gamma[i]))
        p2 = RES_t[i]/tau - 2.0*RES_beta[i]*(tau - RES_gamma[i])
        t3 += RES_n[i]*pow(delta,RES_d[i])*pow(tau,RES_t[i])*exp(p1)*p2
    
    t4 = 0.0
    for i in range(55,57):
        THETA = (1.0 - tau) + RES_A[i]*pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i])
        DELTA = THETA*THETA + RES_B[i]*pow((delta - 1.0)*(delta - 1.0),RES_a[i])
        PSI   = exp(-RES_C[i]*((delta - 1.0)*(delta - 1.0)) - RES_D[i]*((tau - 1.0)*(tau - 1.0)))

        D_DELTA_bi_tau = -2.0*THETA*RES_b[i]*pow(DELTA,RES_b[i] - 1.0)
        D_PSI_tau = -2.0*RES_D[i]*(tau - 1.0)*PSI
        t4 += RES_n[i]*delta*(D_DELTA_bi_tau*PSI + pow(DELTA,RES_b[i])*D_PSI_tau)


    return t1 + t2 + t3 + t4


def RES_phi_tau_tau(delta,tau):
    t1 = 0.0
    for i in range(1,8):
        t1 += RES_n[i]*RES_t[i]*(RES_t[i] - 1.0)*pow(delta,RES_d[i])*pow(tau,RES_t[i] - 2.0)

    t2 = 0.0
    for i in range(8,52):
        t2 += RES_n[i]*RES_t[i]*(RES_t[i] - 1.0)*pow(delta,RES_d[i])*pow(tau,RES_t[i] - 2.0)*exp(-pow(delta,RES_c[i]))

    t3 = 0.0
    for i in range(52,55):
        p1 = -RES_alpha[i]*((delta - RES_eps[i])*(delta - RES_eps[i])) - RES_beta[i]*((tau - RES_gamma[i])*(tau - RES_gamma[i]))
        p2 = pow((RES_t[i]/tau - 2.0*RES_beta[i]*(tau - RES_gamma[i])),2.0) - RES_t[i]/(tau*tau) - 2.0*RES_beta[i]
        t3 += RES_n[i]*pow(delta,RES_d[i])*pow(tau,RES_t[i])*exp(p1)*p2

    t4 = 0.0
    for i in range(55,57):
        THETA = (1.0 - tau) + RES_A[i]*pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i])
        DELTA = THETA*THETA + RES_B[i]*pow((delta - 1.0)*(delta - 1.0),RES_a[i])
        PSI   = exp(-RES_C[i]*((delta - 1.0)*(delta - 1.0)) - RES_D[i]*((tau - 1.0)*(tau - 1.0)))

        D_DELTA_bi_tau = -2.0*THETA*RES_b[i]*pow(DELTA,RES_b[i] - 1.0)
        D_PSI_tau = -2.0*RES_D[i]*(tau - 1.0)*PSI

        D2_DELTA_bi_tau = 2.0*RES_b[i]*pow(DELTA,RES_b[i] - 1.0) + 4.0*(THETA*THETA)*RES_b[i]*(RES_b[i] - 1.0)*pow(DELTA,RES_b[i] - 2.0)
        D2_PSI_tau = (2*RES_D[i]*((tau - 1.0)*(tau - 1.0)) - 1.0)*2.0*RES_D[i]*PSI

        t4 += RES_n[i]*delta*(D2_DELTA_bi_tau*PSI + 2.0*D_DELTA_bi_tau*D_PSI_tau + pow(DELTA,RES_b[i])*D2_PSI_tau)


    return t1 + t2 + t3 + t4


def RES_phi_delta_tau(delta,tau):
    t1 = 0.0
    for i in range(1,8):
        t1 += RES_n[i]*RES_d[i]*RES_t[i]*pow(delta,RES_d[i] - 1.0)*pow(tau,RES_t[i] - 1.0)

    t2 = 0.0
    for i in range(8,52):
        p1  = (RES_d[i] - RES_c[i]*pow(delta,RES_c[i]))*exp(-pow(delta,RES_c[i]))
        t2 += RES_n[i]*RES_t[i]*pow(delta,RES_d[i] - 1.0)*pow(tau,RES_t[i] - 1.0)*p1

    t3 = 0.0
    for i in range(52,55):
        p1 = -RES_alpha[i]*((delta - RES_eps[i])*(delta - RES_eps[i])) - RES_beta[i]*((tau - RES_gamma[i])*(tau - RES_gamma[i]))
        p2 = (-RES_d[i]/delta - 2.0*RES_alpha[i]*(delta - RES_eps[i]))*(RES_t[i]/tau - 2.0*RES_beta[i]*(tau - RES_gamma[i]))
        t3 += RES_n[i]*pow(delta,RES_d[i])*pow(tau,RES_t[i])*exp(p1)*p2

    t4 = 0.0
    for i in range(55,57):
        THETA = (1.0 - tau) + RES_A[i]*pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i])
        DELTA = THETA*THETA + RES_B[i]*pow((delta - 1.0)*(delta - 1.0),RES_a[i])
        PSI   = exp(-RES_C[i]*((delta - 1.0)*(delta - 1.0)) - RES_D[i]*((tau - 1.0)*(tau - 1.0)))

        D_DELTA_bi_tau = -2.0*THETA*RES_b[i]*pow(DELTA,RES_b[i] - 1.0)
        D_PSI_tau = -2.0*RES_D[i]*(tau - 1.0)*PSI

        D2_DELTA_bi_tau = 2.0*RES_b[i]*pow(DELTA,RES_b[i] - 1.0) + 4.0*(THETA*THETA)*RES_b[i]*(RES_b[i] - 1.0)*pow(DELTA,RES_b[i] - 2.0)
        D2_PSI_tau = (2*RES_D[i]*((tau - 1.0)*(tau - 1.0)) - 1.0)*2.0*RES_D[i]*PSI

        D_DELTA_delta = (delta - 1.0)*(RES_A[i]*THETA*(2.0/RES_beta[i])*pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i] - 1.0) + 2.0*RES_B[i]*RES_a[i]*pow((delta - 1.0)*(delta - 1.0),RES_a[i]-1.0))
        D_DELTA_bi_delta = RES_b[i]*pow(DELTA,RES_b[i]-1.0)*D_DELTA_delta
        D_PSI_delta = -2.0*RES_C[i]*(delta - 1.0)*PSI

        D2_PSI_delta_tau = 4.0*RES_C[i]*RES_D[i]*(delta - 1.0)*(tau - 1.0)*PSI
        
        p1 = RES_A[i]*RES_b[i]*(2.0/RES_beta[i])*pow(DELTA,RES_b[i] - 1.0)*(delta - 1.0)*pow((delta - 1.0)*(delta - 1.0),0.5/RES_beta[i] - 1.0)
        p2 = 2.0*THETA*RES_b[i]*(RES_b[i] - 1.0)*pow(DELTA,RES_b[i] - 1.0)*D_DELTA_delta
        D2_DELTA_bi_delta_tau = -p1 - p2
        
        p1 = pow(DELTA,RES_b[i])*(D_PSI_tau + delta*D2_PSI_delta_tau)
        p2 = delta*D_DELTA_bi_tau*D_PSI_tau
        p3 = D_DELTA_bi_tau*(PSI + delta*D_PSI_delta)
        p4 = D2_DELTA_bi_delta_tau*delta*PSI
        
        t4 += RES_n[i]*(p1 + p2 + p3 + p4)

    return t1 + t2 + t3 + t4

def pressure(rho,T):
    delta = rho/CRITICAL_DENSITY
    tau = CRITICAL_TEMPERATURE/T
    return rho*SPECIFIC_GAS_CONSTANT*T*(1.0 + delta*RES_phi_delta(delta,tau))

def internal_energy(rho,T):
    delta = rho/CRITICAL_DENSITY
    tau = CRITICAL_TEMPERATURE/T
    return tau*SPECIFIC_GAS_CONSTANT*T*(IG_phi_tau(delta,tau) + RES_phi_tau(delta,tau))

def entropy(rho,T):
    delta = rho/CRITICAL_DENSITY
    tau = CRITICAL_TEMPERATURE/T
    return SPECIFIC_GAS_CONSTANT*(tau*(IG_phi_tau(delta,tau) + RES_phi_tau(delta,tau)) - IG_phi(delta,tau) - RES_phi(delta,tau))

def enthalpy(rho,T):
    delta = rho/CRITICAL_DENSITY
    tau = CRITICAL_TEMPERATURE/T
    return (SPECIFIC_GAS_CONSTANT*T)*(1.0 + tau*(IG_phi_tau(delta,tau) + RES_phi_tau(delta,tau)) + delta*RES_phi_delta(delta,tau))

def isochoric_heat_capacity(rho,T):
    delta = rho/CRITICAL_DENSITY
    tau = CRITICAL_TEMPERATURE/T
    return (-(tau*tau)*(IG_phi_tau_tau(delta,tau) + RES_phi_tau_tau(delta,tau)))*SPECIFIC_GAS_CONSTANT

def isobaric_heat_capacity(rho,T):
    delta = rho/CRITICAL_DENSITY
    tau = CRITICAL_TEMPERATURE/T
    p1 = -(tau*tau)*(IG_phi_tau_tau(delta,tau) + RES_phi_tau_tau(delta,tau))
    p2 = pow((1.0 + delta*RES_phi_delta(delta,tau) - delta*tau*RES_phi_delta_tau(delta,tau)),2.0)
    p3 = 1.0 + 2.0*delta*RES_phi_delta(delta,tau) + (delta*delta)*RES_phi_delta_delta(delta,tau)
    return (p1 + p2/p3)*SPECIFIC_GAS_CONSTANT

def speed_of_sound(rho,T):
    delta = rho/CRITICAL_DENSITY
    tau = CRITICAL_TEMPERATURE/T
    p1 = 1.0 + 2.0*delta*RES_phi_delta(delta,tau) + (delta*delta)*RES_phi_delta_delta(delta,tau)
    p2 = pow((1.0 + delta*RES_phi_delta(delta,tau) - delta*tau*RES_phi_delta_tau(delta,tau)),2.0)
    p3 = (tau*tau)*(IG_phi_tau_tau(delta,tau) + RES_phi_tau_tau(delta,tau))
    return sqrt((SPECIFIC_GAS_CONSTANT*T)*(p1 - p2/p3)*1000.0)

def joule_thompson_coefficient(rho,T):
    delta = rho/CRITICAL_DENSITY
    tau = CRITICAL_TEMPERATURE/T
    p1 = delta*RES_phi_delta(delta,tau) + (delta*delta)*RES_phi_delta_delta(delta,tau) + (delta*tau)*RES_phi_delta_tau(delta,tau)
    p2 = pow(1.0 + delta*RES_phi_delta(delta,tau) - delta*tau*RES_phi_delta_tau(delta,tau),2.0)
    p3 = IG_phi_tau_tau(delta,tau) + RES_phi_tau_tau(delta,tau)
    p4 = 1.0 + 2.0*delta*RES_phi_delta(delta,tau) + (delta*delta)*RES_phi_delta_delta(delta,tau)
    return (-p1/(p2 - (tau*tau)*p3*p4))/(SPECIFIC_GAS_CONSTANT*rho)

def isothermal_throttling_coefficient(rho,T):
    delta = rho/CRITICAL_DENSITY
    tau = CRITICAL_TEMPERATURE/T
    p1 = 1.0 + delta*RES_phi_delta(delta,tau)-delta*tau*RES_phi_delta_tau(delta,tau)
    p2 = 1.0 + 2.0*delta*RES_phi_delta(delta,tau) + (delta*delta)*RES_phi_delta_delta(delta,tau)
    return (1.0 - p1/p2)/rho

def isentropic_temperature_pressure_coefficent(rho,T):
    delta = rho/CRITICAL_DENSITY
    tau = CRITICAL_TEMPERATURE/T
    p1 = 1.0 + delta*RES_phi_delta(delta,tau) - delta*tau*RES_phi_delta_tau(delta,tau)
    p2 = pow(1.0 + delta*RES_phi_delta(delta,tau) - delta*tau*RES_phi_delta_tau(delta,tau),2.0)
    p3 = IG_phi_tau_tau(delta,tau) + RES_phi_tau_tau(delta,tau)
    p4 = 1.0 + 2.0*delta*RES_phi_delta(delta,tau) + (delta*delta)*RES_phi_delta_delta(delta,tau)
    return (p1/(p2 - (tau*tau)*p3*p4))/(SPECIFIC_GAS_CONSTANT*rho)

def phase_equilibrium(rho1,rho2,T,P):
    delta1 = rho1/CRITICAL_DENSITY
    delta2 = rho2/CRITICAL_DENSITY
    tau = CRITICAL_TEMPERATURE/T

    eqn1 = P/(SPECIFIC_GAS_CONSTANT*T*rho1) - 1.0 - delta1*RES_phi_delta(delta1,tau)
    eqn2 = P/(SPECIFIC_GAS_CONSTANT*T*rho2) - 1.0 - delta2*RES_phi_delta(delta2,tau)
    a = RES_phi(delta1,tau)
    b = RES_phi(delta2,tau)
    eqn3 = (P/(SPECIFIC_GAS_CONSTANT*T))*((1.0/rho2) - (1.0/rho1)) - log(rho1/rho2) - a + b
    return (eqn1,eqn2,eqn3)


if __name__ == '__main__':
    #Verification tests
    #Computed results must match the results in the relevant tables.
    
    delta = 838.025/CRITICAL_DENSITY
    tau = CRITICAL_TEMPERATURE/500.0

    print 'IAPWS95: TABLE 6 (IDEAL AND RESIDUAL)'
    print '-------------------------------'
    print 'PHI               %1.9G \t %1.9G' % (IG_phi(delta,tau),RES_phi(delta,tau))
    print 'PHI_delta         %1.9G \t %1.9G' % (IG_phi_delta(delta,tau),RES_phi_delta(delta,tau))
    print 'PHI_delta_delta  %1.9G \t  %1.9G' %  (IG_phi_delta_delta(delta,tau),RES_phi_delta_delta(delta,tau))
    print 'PHI_tau           %1.9G \t %1.9G' % (IG_phi_tau(delta,tau),RES_phi_tau(delta,tau))
    print 'PHI_tau_tau      %1.9G \t %1.9G' %  (IG_phi_tau_tau(delta,tau),RES_phi_tau_tau(delta,tau))
    print 'PHI_delta_tau     %1.9G \t %1.9G' % (IG_phi_delta_tau(delta,tau),RES_phi_delta_tau(delta,tau))
    
    print 'IAPWS95: TABLE 7 (VERIFICATION)'
    print '-------------------------------'
    T = 300.0
    rho = 0.996556e3
    print '%d %1.9G %1.9G %1.9G %1.9G %1.9G'% (T,rho,pressure(rho,T)*1e-3,isochoric_heat_capacity(rho,T),speed_of_sound(rho,T),entropy(rho,T))
    rho = 0.1005308e4
    print '%d %1.9G %1.9G %1.9G %1.9G %1.9G'% (T,rho,pressure(rho,T)*1e-3,isochoric_heat_capacity(rho,T),speed_of_sound(rho,T),entropy(rho,T))
    rho = 0.1188202e4
    print '%d %1.9G %1.9G %1.9G %1.9G %1.9G'% (T,rho,pressure(rho,T)*1e-3,isochoric_heat_capacity(rho,T),speed_of_sound(rho,T),entropy(rho,T))

    print '\n'
    
    T = 500.0
    rho = 0.435
    print '%d %1.9G %1.9G %1.9G %1.9G %1.9G'% (T,rho,pressure(rho,T)*1e-3,isochoric_heat_capacity(rho,T),speed_of_sound(rho,T),entropy(rho,T))
    rho = 0.4352e1
    print '%d %1.9G %1.9G %1.9G %1.9G %1.9G'% (T,rho,pressure(rho,T)*1e-3,isochoric_heat_capacity(rho,T),speed_of_sound(rho,T),entropy(rho,T))
    rho = 0.838025e3
    print '%d %1.9G %1.9G %1.9G %1.9G %1.9G'% (T,rho,pressure(rho,T)*1e-3,isochoric_heat_capacity(rho,T),speed_of_sound(rho,T),entropy(rho,T))
    rho = 0.1084564e4
    print '%d %1.9G %1.9G %1.9G %1.9G %1.9G'% (T,rho,pressure(rho,T)*1e-3,isochoric_heat_capacity(rho,T),speed_of_sound(rho,T),entropy(rho,T))

    print '\n'

    T = 647.0
    rho = 0.358e3
    print '%d %1.9G %1.9G %1.9G %1.9G %1.9G'% (T,rho,pressure(rho,T)*1e-3,isochoric_heat_capacity(rho,T),speed_of_sound(rho,T),entropy(rho,T))

    print '\n'

    T = 900.0
    rho = 0.241
    print '%d %1.9G %1.9G %1.9G %1.9G %1.9G'% (T,rho,pressure(rho,T)*1e-3,isochoric_heat_capacity(rho,T),speed_of_sound(rho,T),entropy(rho,T))
    rho = 0.526150e2
    print '%d %1.9G %1.9G %1.9G %1.9G %1.9G'% (T,rho,pressure(rho,T)*1e-3,isochoric_heat_capacity(rho,T),speed_of_sound(rho,T),entropy(rho,T))
    rho = 0.870769e3
    print '%d %1.9G %1.9G %1.9G %1.9G %1.9G'% (T,rho,pressure(rho,T)*1e-3,isochoric_heat_capacity(rho,T),speed_of_sound(rho,T),entropy(rho,T))

    print '\n'
    
    print 'IAPWS95: TABLE 8 (2-PHASE)'
    print '--------------------------'
    T = 275.0,450.0,625.0
    rhoL = 0.999887406e3,0.890341250e3,0.567090385e3
    rhoV = 0.550664919e-2,0.481200360e1,0.118290280e3
    
    print '%1.9G %1.9G %1.9G' % (pressure(rhoL[0],T[0]),pressure(rhoL[1],T[1]),pressure(rhoL[2],T[2]))
    print '%1.9G %1.9G %1.9G' % (pressure(rhoV[0],T[0]),pressure(rhoV[1],T[1]),pressure(rhoV[2],T[2]))
    print '%1.9G %1.9G %1.9G' % rhoL
    print '%1.9G %1.9G %1.9G' % rhoV
    print '%1.9G %1.9G %1.9G' % (enthalpy(rhoL[0],T[0]),enthalpy(rhoL[1],T[1]),enthalpy(rhoL[2],T[2]))
    print '%1.9G %1.9G %1.9G' % (enthalpy(rhoV[0],T[0]),enthalpy(rhoV[1],T[1]),enthalpy(rhoV[2],T[2]))
    print '%1.9G %1.9G %1.9G' % (entropy(rhoL[0],T[0]),entropy(rhoL[1],T[1]),entropy(rhoL[2],T[2]))
    print '%1.9G %1.9G %1.9G' % (entropy(rhoV[0],T[0]),entropy(rhoV[1],T[1]),entropy(rhoV[2],T[2]))

    print 'Phase equilibrium'
    print  phase_equilibrium(rhoL[0],rhoV[0],T[0],pressure(rhoL[0],T[0]))
    
