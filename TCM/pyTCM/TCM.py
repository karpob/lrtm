def runTCM(inputPar):
        """
        input pars: dz, 
        XHe_i,
        XH2S_i,
        XNH3_i,
        XH2O_i,
        XCH4_i,
        XPH3_i,
        XCO,
        P_temp,
        T_temp,
        g0_i,
        R0_i,
        P0_i,
        T_targ_i,
        P_targ_i,
        P_term_i,
        use_lindal_in,
        n_lindal_pts_i,
        SuperSatSelf1_i,
        SuperSatSelf2_i,
        SuperSatSelf3_i,
        SuperSatSelf4_i,
        supersatNH3_i,
        supersatH2S_i,
        AutoStep_constant,
        Hydrogen_Curve_Fit_Select,  
        dP_init,
        dP_fine,
        P_fine_start,
        P_fine_stop,
        use_dz,
        frain,
        select_ackerman
        """     
        from LAYERS import init_atm
        from LAYERS import new_layer       
        eflag=0 	  
        layer,others = init_atm(inputPar)
	import pylab 

        #while error in temperature is greater than Tlimit, or if you've gone on long enough
        #  without converging on a solution ...  
        T_err=100.
        TLIMIT=0.0001
        MAXTRIES=100
        MAXLAYERS=1e6
        T_targ=max(others['TfL'])
        jcntr=0
        while (abs(T_err) > TLIMIT and jcntr < MAXTRIES):
      
                eflag = 0; #reset flag
                jj=0;     #reset jj
                jcntr+=1;  #increment jcntr
           		
                print "Iteration: %d      \n"%jcntr
            	j=1	
                while(j<MAXLAYERS and eflag!=99): #go through the layers
                        eflag,dP,layer=new_layer(j,eflag, inputPar,layer,others)
                        P = layer['P'][j]
                        if (eflag == 98 ):#  #check target temperature 
                                # eflag       P      #
                                eflag = 97                         #  97      < P_targ 
                                
                                """
                                REPLACE FOLLOWING LINE
                                """          
                                T = layer['T'][j]                    #  99      <=P_term  
                                T_err = 100.*(T - T_targ)/T_targ  #  98      <=P_targ  
                                print "Terr %f %f %d \n"%(T_err,TLIMIT,j)
			        if (abs(T_err)>TLIMIT):
			                for key in layer.keys():
			                        try:print len(layer[key]),key
			                        except:pass
                                        
                                        inputPar['T_temp'] = layer['T'][0] - (T - T_targ)*2.0;
                                        inputPar['P_temp'] = layer['P'][0];
                                        print inputPar['T_temp'],inputPar['P_temp']
                                        layer,others = init_atm(inputPar)
                                        j=0
                                        eflag = 99;
                                        
		        j+=1			
      

      
        #Parr=numpy.asarray(layer['P']))
        #Zarr=numpy.asarray(layer['z']))
        #z_offset=min(Zarr[Parr>=P0])
        
        return layer,others		  



