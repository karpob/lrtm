def DeBoerTCM(TP_list,TP_force,XH2S_i,XNH3_i,XH2O_i,XCH4_i,
                                XPH3_i,XHe_i,XCO,P_temp,T_temp, g0_i,R0e_i,
                                P0_i,T_targ_i,P_targ_i,P_term_i,
                                use_lindal,SuperSatSelf_H2S,SuperSatSelf_NH3,
                                SuperSatSelf_PH3,SuperSatSelf_H2O,supersatNH3,
                                supersatH2S,AutoStep_constant,fp,dz,oblateness_factor,use_dz,
                                dP_init,dP_fine,P_fine_start,P_fine_stop,frain,select_ackerman):
        import shutil
        from TCM import TCM                        
        # function DeBoer_TCM is a matlab wrapper for the DeBoer Thermo chemical model. The function TCM is a mex
        # module which must be compiled. The source/mex module is located in TCM_mex. TCM module is compiled by:
        # >mex *.C -o TCM
        # Any warnings issued by the mex compiler may be ignored. If you get errors, then you need to worry.
        # 
        # Variable definitions:
        #        -->  INPUT
        #             ->        TP_list: a vector of strings which represent different temperature pressure profiles.
        #                               TP_list(1): Seiff Jupiter TP profile from Galileo proble
        #                               TP_list(2): Lindal Jupiter TP profile from Voyager
        #                               TP_list(3): Lindal Saturn TP profile from Voyager
        #                               TP_list(4): Lindal Uranus TP profile from Voyager
        #                               TP_list(5): Lindal Neptune TP profile from Voyager
        #                               TP_list(6): Load TP profile TP.TCM in the TCM_mex directory
        #
        #               ->       TP_force: a string which selects the TP profile from the TP_list vector.
        #               ->       XH2S_i: absolute deep abundance of hydrogen sulfide 
        #               ->       XNH3_i: absolute deep abundance of ammonia
        #               ->       XH2O_i: absolute deep abundance of water vapor
        #               ->       XCH4_i: absolute deep abundance of methane
        #               ->       XPH3_i: absolute deep abundance of phosphine
        #               ->       XHe_i: absolute deep abundance of helium
        #               ->       XCO_i: absolute deep abundance of carbon monoxide
        #               ->       P_temp: The deep pressure level (bottom)
        #               ->       T_temp: The Temperature at deep pressure level (bottom)
        #               ->       P_term_i: The lowest pressure level (top)
        #               ->       T_term_i: The temperature at lowest pressure level (top)
        #               ->       use_lindal: Flag to use TP profile as a guess/forcing
        #               ->       SuperSatSelf_H2S: abundance of hydrogen sulfide in supersaturation
        #               ->       SuperSatSelf_PH3: abundance of phosphine in supersaturation
        #               ->       SuperSatSelf_H2O: abundance of water vapor in supersaturation
        #               ->       SuperSatNH3: abundace of ammonia in supersaturation (?for solution cloud?)
        #               ->       SuperSatH2S: abundacne of hydrogen sulfide in supersaturation (?for solution cloud?)
        #               ->       Autostep_constant: Mysterious constant semi-empirically derived from scale height  
        #               ->       fp: ortho-para hydrogen selector 0.0=equilibrium,  -1.0=intermediate, 0.25 = normal
        #               ->       dz: altitude step, if you go with autostep, set this to 0
        #               ->       oblateness_factor: ratio of polar to equitorial radius
        #               ->       use_dz: Set this to 1 if you want to step using a step in altitude, if you want a step in
        #                                pressure set this to 0
        #               ->       dP_init: initial coarse pressure step
        #               ->       dP_fine: fine pressure step
        #               ->       P_fine_start: Pressure level to start stepping with fine dP
        #               ->       P_fine_stop: Pressure level to stop stepping with fine dP 
        #
        #  OUTPUT
        #              <-        me:length of the tcme array (number of pressure levels)
        #              <-        tcme: Thermo-chemical model output variables equitorial
        #                             tcme(1,1:me): Pressure levels (bars)
        #                             tcme(2,1:me): Temperature (Kelvin)
        #                             tcme(3,1:me): Altitude (km)
        #                             tcme(4,1:me): Hydrogen (H2) abundance
        #                             tcme(5,1:me): Helium abundance
        #                             tcme(6,1:me): hydrogen sulfide abundance
        #                             tcme(7,1:me): ammonia abundance 
        #                             tcme(8,1:me): water vapor abundance
        #                             tcme(9,1:me): methane abundance
        #                             tcme(10,1:me): phospine abundance
        #                             tcme(11,1:me): DeBoer cloud mask variable (long int) (ice phase, liquid, no cloud, etc)
        #                             tcme(12,1:me): Ammonium hydrosulfide cloud 
        #                             tcme(13,1:me): Hydrogen sulfide cloud density (g/cm^3)
        #                             tcme(14,1:me): Ammonia ice cloud density (g/cm^3)
        #                             tcme(15,1:me): Water ice cloud density (g/cm^3)
        #                             tcme(16,1:me): Methane cloud density (g/cm^3)
        #                             tcme(17,1:me): Phosphine cloud density (g/cm^3)
        #                             tcme(18,1:me): H2O-NH3 solution cloud density (g/cm^3)
        #                             tcme(19,1:me): Gravity (cm/s^2)
        #                             tcme(20,1:me): Average mass mu (Atomic Mass unit AMU)
        #                             tcme(21,1:me): refractivity with H2 and He only *isn't used in maintamone*
        #                             tcme(22,1:me): refractivity with H2, He, H2S, PH3, NH3, CH4, and H2O
        #                                             R = REFRACT_X*PX*(293/T) 
        #                                             n = (R/1E6) + 1 
        #               <-       tcmp: Thermo-chemical model output variables polar (same as tcme with scaled altitude)
        #               <-       DSOL_NH3: Density of Solution cloud which is ammonia in solution
        #                      *note* refractivity isn't used in maintamone it is a legacy from the DeBoer TCM model
        #

        if(use_lindal=='Y'):
                lindal_profile_switch=1
        else:
                lindal_profile_switch=0
        
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        # Copy over desired Temp Pressure profile
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        if(TP_list[0]==TP_force):
                shutil.copyfile('TCM/TP_Seiff_Jupiter.JUP', 'TP.TCM')
                f=open('TP.TCM','r')
                lines=f.readlines()
                n_lindal=len(lines)    
                f.close()
        elif(TP_list[1]==TP_force):
                shutil.copyfile('TCM/TP.JUP','TP.TCM')
                f=open('TP.TCM','r')
                lines=f.readlines()
                n_lindal=len(lines)
                f.close()
        elif(TP_list[2]==TP_force):
                shutil.copyfile('TCM/TP.SAT TP.TCM')
                f=open('TP.TCM','r')
                lines=f.readlines()
                n_lindal=len(lines)
                f.close()
        elif(TP_list[3]==TP_force):
                shutil.copyfile('TCM/TP.URN', 'TP.TCM')
                f=open('TP.TCM','r')
                lines=f.readlines()
                n_lindal=len(lines)
                f.close()
        elif(TP_list[4]==TP_force):
                shutil.copyfile('TCM/TP.NEP', 'TP.TCM')
                f=open('TP.TCM','r')
                lines=f.readlines()
                n_lindal=len(lines)
                f.close()
        elif(TP_list[5]==TP_force):
                TP_in_directory=load('TP.TCM')
                f=open('TP.TCM','r')
                lines=f.readlines()
                n_lindal=len(lines)
                f.close()
        else:
                print 'Invalid TP profile.'
        
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        # Run DeBoer TCM shared object (must be compiled first)
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

        #[P,T,XH2,XHe,XH2S,XNH3,XH20,XCH4,XPH3,
        #clouds,DNH4SH_i,DH2S_i,DNH3_i,DH2O_i,DCH4_i,DPH3_i,DSOL_i,
        #g,mu,ref_w_o,ref_w,z,DSOL_NH3,Preal]
        
        nlevels=TCM.intoTheVoid(dz,XHe_i,XH2S_i,XNH3_i,XH2O_i,XCH4_i,XPH3_i,XCO,
                              P_temp,T_temp,g0_i,R0e_i,P0_i,T_targ_i,P_targ_i,P_term_i,1,n_lindal,
                              SuperSatSelf_H2S,SuperSatSelf_NH3,SuperSatSelf_PH3,
                              SuperSatSelf_H2O,supersatNH3,supersatH2S,
                              AutoStep_constant,fp,dP_init,dP_fine,P_fine_start,P_fine_stop,use_dz,frain,select_ackerman)
        layerKeys=['P','T','z','XH2','XHe','XH2S','XNH3','XH2O','XCH4','XPH3',
                   'DNH4SH','DH2S','DNH3','DH2O','DCH4','DPH3','DSOL',
                   'g','mu','DSOL_NH3','P_real']
	layer={}
	for key in layerKeys:layer[key]=[]			

	for j in range(0, nlevels+1):
		for i in range(0,len(layerKeys)):
			layer[layerKeys[i]].append(TCM.getFloatValues(i,j))
        me,n=numpy.shape(P)

        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        # Remove fictional clouds from DeBoer TCM using built-in filter 'clouds'
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        #[DNH4SH,DH2S,DNH3,DH2O,DCH4,DPH3,DSOL]=filter_clouds(clouds,DNH4SH_i,DH2S_i,DNH3_i,DH2O_i,DCH4_i,DPH3_i,DSOL_i)

        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        #  Scale Polar profile according to oblateness factor
        # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
        clouds=numpy.zeros([me])
        #Equatorial
        tcme[0:me,0:23]=numpy.r_[layer['P'][0:me],layer['T'][0:me],layer['z'][0:me],layer['XH2'][0:me],layer['XHe'][0:me],layer['XH2S'][0:me],layer['XNH3'][0:me],layer['XH2O'][0:me],layer['XCH4'][0:me],layer['XPH3'][0:me],
                 clouds[0:me],layer['DNH4SH'][0:me],layer['DH2S'][0:me],layer['DNH3'][0:me],layer['DH2O'][0:me],layer['DCH4'][0:me],layer['DPH3'][0:me],layer['DSOL'][0:me],
                 layer['g'][0:me],layer['mu'][0:me],clouds,clouds,layer['P_real'][0:me]]
                 #ref_w_o[0:me],ref_w[0:me],Preal[0:me]]
        #Polar
        tcmp[0:me,0:23]=numpy.r_[tcme[0:me,0:2],oblateness_factor*tcme[0:me,2],tcme[0:me,3:23]]
        return me,tcme,tcmp,DSOL_NH3
 # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
