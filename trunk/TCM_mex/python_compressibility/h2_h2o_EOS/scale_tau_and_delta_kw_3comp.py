def scale_tau_and_delta_kw_4comp(T,rho,x1,x2,x3,x4,T_c1,T_c2,T_c3,T_c4,rho_c1,rho_c2,rho_c3,rho_c4,BetaT1,BetaV1,GammaT1,GammaV1,BetaT2,BetaV2,GammaT2,GammaV2,BetaT3,BetaV3,GammaT3,GammaV3,BetaT4,BetaV4,GammaT4,GammaV4,BetaT5,BetaV5,GammaT5,GammaV5,BetaT6,BetaV6,GammaT6,GammaV6):
        from scale_tau_and_delta_kw import scale_tau_and_delta_kw
        [tau,delta,Tc12,rho12]=scale_tau_and_delta_kw(T,rho,x1,x2,T_c1,T_c2,rho_c1,rho_c2,BetaT1,BetaV1,GammaT1,GammaV1)
        [tau,delta,Tc13,rho13]=scale_tau_and_delta_kw(T,rho,x1,x3,T_c1,T_c3,rho_c1,rho_c3,BetaT2,BetaV2,GammaT2,GammaV2)
        [tau,delta,Tc14,rho14]=scale_tau_and_delta_kw(T,rho,x1,x4,T_c1,T_c4,rho_c1,rho_c4,BetaT3,BetaV3,GammaT3,GammaV3)
        
        [tau,delta,Tc23,rho23]=scale_tau_and_delta_kw(T,rho,x1,x2,T_c2,T_c3,rho_c2,rho_c3,BetaT4,BetaV4,GammaT4,GammaV4)
        [tau,delta,Tc24,rho24]=scale_tau_and_delta_kw(T,rho,x1,x3,T_c2,T_c4,rho_c2,rho_c4,BetaT5,BetaV5,GammaT5,GammaV5)
        
        [tau,delta,Tc34,rho34]=scale_tau_and_delta_kw(T,rho,x1,x4,T_c3,T_c4,rho_c3,rho_c4,BetaT6,BetaV6,GammaT6,GammaV6)
        
        Tc=Tc12+Tc13+Tc14+Tc23+Tc24+Tc34
        rho_inv=(1.0/rho12)+(1.0/rho13)+(1.0/rho14)+(1.0/rho23)+(1.0/rho24)+(1.0/rho34)
        rho_c=1/rho_inv
        tau=Tc/T
        delta=rho/rho_c
        return tau,delta,Tc,rho_c
