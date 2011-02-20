def scale_tau_and_delta_kw_4comp(T,rho,x1,x2,x3,x4,T_c1,T_c2,T_c3,T_c4,rho_c1,rho_c2,rho_c3,rho_c4,BetaT1,BetaV1,GammaT1,GammaV1,BetaT2,BetaV2,GammaT2,GammaV2,BetaT3,BetaV3,GammaT3,GammaV3,BetaT4,BetaV4,GammaT4,GammaV4,BetaT5,BetaV5,GammaT5,GammaV5,BetaT6,BetaV6,GammaT6,GammaV6):
        from scale_tau_and_delta_kw import scale_tau_and_delta_kw
        from numpy import inf,ones,shape
        
        #nominally these relate to correlations according toe GERG 2004, and Karpowicz's correlation:
        #hydrogen-methane
        [tau,delta,Tc12,rho12]=scale_tau_and_delta_kw(T,rho,x1,x2,T_c1,T_c2,rho_c1,rho_c2,BetaT1,BetaV1,GammaT1,GammaV1) 
        #hydrogen-water
        [tau,delta,Tc13,rho13]=scale_tau_and_delta_kw(T,rho,x1,x3,T_c1,T_c3,rho_c1,rho_c3,BetaT2,BetaV2,GammaT2,GammaV2) 
        
        #nominally there Betas and Gammas should be one (ideal mixing), unless you've worked on a correlation for them (I didn't)
        #hydrogen-helium
        [tau,delta,Tc14,rho14]=scale_tau_and_delta_kw(T,rho,x1,x4,T_c1,T_c4,rho_c1,rho_c4,BetaT3,BetaV3,GammaT3,GammaV3) 
        #methane-water
        [tau,delta,Tc23,rho23]=scale_tau_and_delta_kw(T,rho,x2,x3,T_c2,T_c3,rho_c2,rho_c3,BetaT4,BetaV4,GammaT4,GammaV4)
        #methane-helium
        [tau,delta,Tc24,rho24]=scale_tau_and_delta_kw(T,rho,x2,x4,T_c2,T_c4,rho_c2,rho_c4,BetaT5,BetaV5,GammaT5,GammaV5) 
        #water-helium
        [tau,delta,Tc34,rho34]=scale_tau_and_delta_kw(T,rho,x3,x4,T_c3,T_c4,rho_c3,rho_c4,BetaT6,BetaV6,GammaT6,GammaV6) 
        
        Tc=Tc12+Tc13+Tc14+Tc23+Tc24+Tc34
        rho_inv1=(1.0/rho12)
        rho_inv2=(1.0/rho13)
        rho_inv3=(1.0/rho14)
        rho_inv4=(1.0/rho23)
        rho_inv5=(1.0/rho24)
        rho_inv6=(1.0/rho34)
        
        #make sure things don't blow up in your face if you've got a 0 mole fraction of something
        rho_inv1[rho_inv1==inf]=0.0
        rho_inv2[rho_inv2==inf]=0.0
        rho_inv3[rho_inv3==inf]=0.0
        rho_inv4[rho_inv4==inf]=0.0
        rho_inv5[rho_inv5==inf]=0.0
        rho_inv6[rho_inv6==inf]=0.0
        
     
        rho_inv=rho_inv1+rho_inv2+rho_inv3+rho_inv4+rho_inv5+rho_inv6
        
        # logic here is a little kludgy, but it makes sense, and is marginally easier than modifying scale_tau_and_delta_kw
        # We have extra x_i*x_i*(1/rho_ci)'s from scale_tau_and_delta_kw depending upon whether or not we have 3 or 4 components
        # if we have 2...no problem, otherwise..we've got issues
        #
        # condition1 will apply if there is x1,x2,and x3 -->subtract off ONE x_i*x_i*(1/rho_ci)
        # condition2 will apply if there is x1,x2,x4--> subtract off  ONE x_i*x_i*(1/rho_ci)
        # now what happens if there is x1,x2,x3,x4? --> subtract off TWO x_i*x_i*(1/rho_ci) (or if you prefer apply condition1 and condition2 applies, so just do them)
        
        
        condition1=(x1!=0.0)&(x2!=0.0)&(x3!=0.0)
        condition2=(x1!=0.0)&(x2!=0.0)&(x4!=0.0)
        
        
        rho_inv[condition1]=rho_inv[condition1]+(-x1[condition1]*x1[condition1]*(1.0/rho_c1)-x2[condition1]*x2[condition1]*(1.0/rho_c2)-x4[condition1]*x4[condition1]*(1.0/rho_c4)-x3[condition1]*x3[condition1]*(1.0/rho_c3))        
        Tc[condition1]=Tc[condition1]+(-x4[condition1]*x4[condition1]*T_c4-x2[condition1]*x2[condition1]*T_c2-x1[condition1]*x1[condition1]*T_c1-x3[condition1]*x3[condition1]*T_c3)
        
        rho_inv[condition2]=rho_inv[condition2]+(-x1[condition2]*x1[condition2]*(1.0/rho_c1)-x2[condition2]*x2[condition2]*(1.0/rho_c2)-x4[condition2]*x4[condition2]*(1.0/rho_c4)-x3[condition2]*x3[condition2]*(1.0/rho_c3))        
        Tc[condition2]=Tc[condition2]+(-x4[condition2]*x4[condition2]*T_c4-x2[condition2]*x2[condition2]*T_c2-x1[condition2]*x1[condition2]*T_c1-x3[condition2]*x3[condition2]*T_c3)
        
        
        rho_c=1.0/rho_inv
        
        
        #Enforce x_i=1.0 for x1,x2,x3,x4
        rho_c[x1==1.0]=rho_c1
        Tc[x1==1.0]=T_c1
       
        rho_c[x2==1.0]=rho_c2
        Tc[x2==1.0]=T_c2
        
        rho_c[x3==1.0]=rho_c3
        Tc[x3==1.0]=T_c3
        
        rho_c[x4==1.0]=rho_c4
        Tc[x4==1.0]=T_c4

        
        #finally calc tau and delta
        tau=Tc/T
        delta=rho/rho_c
                                   
        return tau,delta,Tc,rho_c
