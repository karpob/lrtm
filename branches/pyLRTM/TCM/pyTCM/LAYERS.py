
def init_atm(inputPars):
        """
        #****************************************************************************************************/
        #***************************************************************************************************/
        Initialize the bottom of the atmosphere data structure (in details in model.h) with initial values
        int init_atm()
                Inputs:
                      --> double XHe: Deep Abundance (Mole fraction) of Helium 
                      --> double XH2S: Deep Abundance (Mole fraction) of Hydrogen Sulfide
                      --> double XNH3: Deep Abundance (Mole fraction) of Ammonia
                      --> double XH2O: Deep Abundance (Mole fraction) of Water
                      --> double XCH4: Deep Abundance (Mole fraction) of Methane
                      --> double XPH3: Deep Abundance (Mole fraction) of Phosphine
                      --> double P_temp: Deep Pressure (in bars) of the bottom layer
                      --> double T_temp: Temperature of the bottom layer (in K)
                      --> float g0: Acceleration due to gravity  at pressure level P0
                      --> float R0: Altitude (km) associated with P0
                      --> float P0: Pressure value associated with P0 (bars)
                      --> char use_lindal: select whether or not to use initial Temp/pressure from input file
                      --> float T_targ: Temperature Target value (K)
                      --> float P_targ: Target pressure level associated with T_targ_i (bars)
                      --> float P_term: Top of atmosphere (TOA) value in bars
                      --> int n_lindal_pts: number of values specificed in reference TP profile input file
                      --> float SuperSatSelf1: Amount of supersaturation associated with H2S
                      --> float SuperSatSelf2: Amount of supersaturation associated with NH3
                      --> float SuperSatSelf3: Amount of supersaturation associated with PH3
                      --> float SuperSatSelf4: Amount of supersaturation associated with water
                      --> float supersatNH3: Amount of supersaturation with respect to NH3 associated with NH4SH cloud
                      --> float supersatH2S: Amount of supersaturation with respect to H2S associated with NH4SH cloud
        
                Output:
                      <-- Returns a layer data structure, and others data structure
        /****************************************************************************************************/
        /****************************************************************************************************/
        """
        from numpy import zeros
        from gravity import gravity
        from modelParams import AMU_H2,AMU_NH3,AMU_H2O,AMU_He,AMU_H2S,AMU_CH4,AMU_PH3
        IN=inputPars
        layer={}
        ZERO=1e-30       
        # ***** place abundances in layer data structure *****#	  
        layer['XHe'] = [IN['XHe'],]
        layer['XH2S'] =[IN['XH2S'],]
        layer['XNH3'] = [IN['XNH3'],]
        layer['XH2O'] = [IN['XH2O'],]
        layer['XCH4']= [IN['XCH4'],]
        layer['XPH3'] = [IN['XPH3'],]
        layer['XH2']  = [1.0-IN['XH2S']-IN['XNH3']-IN['XH2O']-IN['XCH4']-IN['XPH3']-IN['XHe']-IN['XCO'],]
        layer['clouds']= ['Nope',]
        layer['DH2S']= [ZERO,]
        layer['DNH3'] = [ZERO,]
        layer['DH2O'] = [ZERO,]
        layer['DCH4'] = [ZERO,]
        layer['DPH3'] = [ZERO,]
        layer['DNH4SH'] = [ZERO,]
        layer['DSOL'] = [ZERO,]
        layer['z'] = [0.0,]
        #*****************************************************#
        #******* Fix deep Temperature and Pressure ********#
        layer['T']=[IN['T_temp'],]
        layer['P']=[IN['P_temp'],]


        #**** Calculate average molecular weight of the atmosphere often called Mu ****#
        layer['mu'] = [AMU_H2*layer['XH2'][0] + AMU_He*IN['XHe'] + AMU_H2S*IN['XH2S'] + AMU_NH3*IN['XNH3'] + AMU_H2O*IN['XH2O'] + AMU_CH4*IN['XCH4'] + AMU_PH3*IN['XPH3'],]
        layer['P0']=IN['P0']
        layer['R0']=IN['R0']
        layer['g0']=IN['g0']
        
        layer['g'] =[gravity(layer,0),] #calculate gravity and put in datastructure

	  
        #******* Initialize Read in Lindal?, target Temp, target pressure, end pressure level ******#
        use_lindal=IN['use_lindal']
        
        #****** Read in temperature pressure profile ******#	  
        if (use_lindal == 'Y'):
                pfp=open('TP.TCM','r')
                lines=pfp.readlines()
                PfL=[]
                TfL=[]
                for line in lines:
                        tP,tT=line.split('    ')
                        tT,junk=tT.split('\n')
                        PfL.append(float(tP))
                        TfL.append(float(tT))
                                    
			
                pfp.close()
                  		
                P_targ=max(PfL)  #set P_target to first element
                T_targ=max(TfL)  #set T_taget to first element
                P_term=min(PfL) #set final pressure to final element
        # Initialize variables with correct values for super saturation
        else:
                P_targ=IN['P_targ']
                T_targ=IN['T_targ']
                P_term=IN['P_term']
                PfL=[]
                TfL=[]        
        supersatNH3=IN['supersatNH3']
        supersatH2S=IN['supersatH2S']
        SuperSatSelf=zeros(4)
        SuperSatSelf[0]=IN['SuperSatSelf1']
        SuperSatSelf[1]=IN['SuperSatSelf2']
        SuperSatSelf[2]=IN['SuperSatSelf3']
        SuperSatSelf[3]=IN['SuperSatSelf4']
        f=open('inputsol.dat','r')
        lines=f.readlines()
        nLists=1
        x1=[]
        a1=[]
        b1=[]
        c1=[]
        
        x2=[]
        a2=[]
        b2=[]
        c2=[]
        
        x3=[]
        a3=[]
        b3=[]
        c3=[]
        
        for line in lines:
                spl=line.split('   ')
                print spl
                x=float(spl[1])
                a=float(spl[2])
                b=float(spl[3])
                try:c=float(spl[4])
                except:c=float(spl[4][0:len(spl[4])-4])
                if(spl[1]=='0.00'):nLists+=1
                        
                if(nLists==1):
                        x1.append(x)
                        a1.append(a)
                        b1.append(b)
                        c1.append(c)
                if(nLists==2):
                        x2.append(x)
                        a2.append(a)
                        b2.append(b)
                        c2.append(c)
                if(nLists==3):
                        x3.append(x)
                        a3.append(a)
                        b3.append(b)
                        c3.append(c)
        f.close()                
        others={}
        others['P_targ']=P_targ
        others['T_targ']=T_targ
        others['P_term']=P_term
        others['supersatNH3']=supersatNH3
        others['supersatH2S']=supersatH2S
        others['SuperSatSelf']=SuperSatSelf
        others['PfL']=PfL
        others['TfL']=TfL
        others['x1']=x1
        others['x2']=x2
        others['x3']=x3
        others['a1']=a1
        others['a2']=a2
        others['a3']=a3
        others['b1']=b1
        others['b2']=b2
        others['b3']=b3
        others['c1']=c1
        others['c2']=c2
        others['c3']=c3
        
	return layer,others



def new_layer(j,
              eflag, 
              inputParams,
              layer,
              others):
              
               
        """
        // void new_layer() 
        //                 This function each calculates how much of a constituent
        //                 is, or isn't condensed for a given atmospheric layer j 
        //                 values are stored in the layer data structure (in model.h)
        //         Inputs:
        //              -->int j : layer index assiociated with layer data structure (see model.h)
        //              -->float dz : altitude step, if one chooses stepping in altitude
        //              -->int *eflag: flag which will cause function to return to main
        //              -->float dP_init: initial pressure step (if one chooses to step in pressure)
        //              -->float dP_fine: pressure step after P_start (if one chooses to step in pressure)
        //              -->float P_fine_start: pressure level to begin stepping with dP_fine
        //              -->float P_fine_stop: pressure level to stop stepping with dP_fine (step with reference input file after this)
        //              -->float frain: Ackerman and Marley value for f_rain loss mechanism
        //              -->float select_ackerman: select whether or not and where to apply Ackerman & Marley cloud loss mechanism 
        //
        //        Output:
        //              <-- values in layer data structure (see model.h)
        /*******************************************************************************/
        """
        from numpy import zeros
        from get_dP_using_dP import get_dP_using_dP
        from get_dT import get_dT
        LX=zeros([6])
        L2X=zeros([6])
#      if (j==1) hereonout=0;
#      /*  Get new P,dP and T,dT values  */

        if (inputParams['use_dz']):eflag,dP,layer = get_dP_using_dz(inputParams['dz'],j,eflag, inputParams,layer,others)         
        else: eflag,dP,layer=get_dP_using_dP(j,eflag, inputParams,layer,others)
        
      
        # Calculated Partial pressures of previous step
      
        dT,layer = get_dT(j,layer['T'][j-1],layer['P'][j],dP,LX,L2X,inputParams['fp'],layer,others);
        return eflag,dP,layer
        """
      
      #dry adiabat
      
      P  = layer[j].P;
      T  = layer[j].T;
      layer[j].first=0;
      PH2  = layer[j-1].XH2*P;
      PHe  = layer[j-1].XHe*P;
      PH2S = layer[j-1].XH2S*P;
      PNH3 = layer[j-1].XNH3*P;
      PPH3 = layer[j-1].XPH3*P;
      PH2O = layer[j-1].XH2O*P;
      PCH4 = layer[j-1].XCH4*P;
      PNH3p= PNH3 - supersatNH3*P; #Partial Pressure of H2S subtracting supersaturation for NH4SH cloud
      PH2Sp= PH2S - supersatH2S*P; #Partial Pressure of H2S subtracting supersaturation for NH4SH cloud
                                   # **Note! supersaturation is defined as an additional mole fraction
                                   #   whereas SuperSatSelf is defined as fraction of supersaturation (ie. 1=100%)
      H = R*T/(layer[j-1].mu*layer[j-1].g); //Scale Height
      
      
      //looking for density used Pideal      
      
      alr = 1e5*dT/(dP*H/layer[j].P);
      dry_adiabatic_lapse_rate=alr;


      /* For the solution clouds see:  Weidenschilling and Lewis 1973, Icarus 20:465.
            Lewis 1969, Icarus 10:365. Briggs and Sackett 1989, Icarus 80:77.
            Atreya and Romani 1985, Recent Advances in Planetary Meteorology.
            Atreya 1986, Atmospheres and Ionospheres of the Outer Planets. */


/********************** If the user selects the solution cloud***********************************
*************************************************************************************************/
/************************************************************************************************
  And so it begins....the endless maze of nested if statements.
  The first set of cases are checking to see whether or not a cloud may form, and calculating the necessary
  values for latent heat and coefficients specified in DeBoer's Thesis
  The next set will be denoted by "Updating Mixing Ratios"
***************************************************************************************************/



/********************** If the user selects the solution cloud***********************************
*************************************************************************************************/
      if (sol_cloud)
      {
            C_sol_NH3 = solution_cloud(T,PNH3,PH2O,&SPNH3,&SPH2O);

            /*This from Briggs and Sackett (from Cook--see also D-230 of CRC)*/
            if (C_sol_NH3 == -1.0)
                  freeze = 273.1;
            else
                  freeze = 273.1 - 124.167*C_sol_NH3 - 189.963*SQ(C_sol_NH3) + 2084.370*CU(C_sol_NH3);

            if (T < freeze)
            {
                  sol_cloud = 0;  /*Water starts to freeze*/
            }
            else if (C_sol_NH3 == -1.0)
                  C_sol_zero=0;
            else if (C_sol_NH3 == 0.0)
                  C_sol_zero=1;
            else
            {
                  CH2O = 1;
                  if (C_sol_NH3!=0.0)
                  {
                        CNH3 = 1;
                   /* Comment out if don't want H2S in solution cloud */
                        C_sol_H2S = h2s_dissolve(j,&SPH2S);
                        if (C_sol_H2S != 0.0)
                        {
                              CH2S = 1;
                        }
                  }
                  LH2O = C_sol_NH3*(4949.75+2022.11*(SQ(C_sol_NH3)-2.0*C_sol_NH3));
                  LH2O+= (1.0-C_sol_NH3)*(5540.48+2022.11*SQ(C_sol_NH3));
                  LH2O*= R/(1.0-C_sol_NH3);
                  LX[3]  = LH2O*layer[j-1].XH2O;
                  L2X[3] = LX[3]*LH2O/(R*T*T);
            }
      }
/***************************************************************************************************/

/*       If the user specifies no solution cloud, or the value for the solution cloud is zero.
          Calculate the saturation vapor pressures for each type of cloud that could form 
***************************************************************************************************/ 
      if(!sol_cloud || C_sol_zero)
      {
            layer[j].DSOL = ZERO;

            /*check for condensation of H2O, H2S and NH3 if out of solution cloud*/
            if (T < TRIPLEPT_H2O)
                  strcpy(phase_H2O,"H2O_over_ice");
            else
                  strcpy(phase_H2O,"H2O_over_water");
            SPH2O = sat_pressure(phase_H2O,T);

            if (T < TRIPLEPT_H2S)
                  strcpy(phase_H2S,"H2S_over_H2S_ice");
            else
                  strcpy(phase_H2S,"H2S_over_liquid_H2S");
            SPH2S = sat_pressure(phase_H2S,T);

            if (T < TRIPLEPT_NH3)
                  strcpy(phase_NH3,"NH3_over_NH3_ice");
            else
                  strcpy(phase_NH3,"NH3_over_liquid_NH3");
            SPNH3 = sat_pressure(phase_NH3,T);

            if(PH2O - SPH2O*SuperSatSelf[3] > SPH2O)
            {
                  CH2O = 1;
                  LH2O = latent_heat(phase_H2O,T);
                  LX[3] = LH2O*layer[j-1].XH2O;
                  L2X[3]= LX[3]*LH2O/(R*T*T);
            }

            if(PH2S - SPH2S*SuperSatSelf[0] > SPH2S)
            {
                  CH2S = 1;
                  LH2S = latent_heat(phase_H2S,T);
                  LX[1] = LH2S*layer[j-1].XH2S;
                  L2X[1]= LX[1]*LH2S/(R*T*T);
            }

            if(PNH3 - SPNH3*SuperSatSelf[1] > SPNH3)
            {
                  CNH3 = 1;
                  LNH3 = latent_heat(phase_NH3,T);
                  LX[2] = LNH3*layer[j-1].XNH3;
                  L2X[2]= LX[2]*LNH3/(R*T*T);
            }
      }

      /*check for condensation of components:  CH4, NH4SH, PH3*/

      if (T < TRIPLEPT_CH4) 
            strcpy(phase_CH4,"CH4_over_CH4_ice");
      else 
            strcpy(phase_CH4,"CH4_over_liquid_CH4");
      SPCH4 = sat_pressure(phase_CH4,T);

      strcpy(phase_PH3,"PH3_over_PH3_ice");
      SPPH3 = sat_pressure(phase_PH3,T);

      KNH4SH = sat_pressure("NH4SH",T);

      if(PH2Sp*PNH3p > KNH4SH) //Does NH4SH cloud form? If so, calculate latent heat.
      {
            if ( !NH4SH_cloud_base )
            {
                  NH4SH_cloud_base=1;
                  diff_P_base = PNH3p - PH2Sp;
            }
            CNH4SH = 1;
            LNH4SH = 9.312E11; /* Atreya book.  1.6E12 Briggs and Sackett */
            LX[0] = 2.0*LNH4SH*PH2Sp*PNH3p/(P*(PH2Sp+PNH3p));
            L2X[0]= LX[0]*5417.0/(T*T);
      }

      if(PCH4 > SPCH4) //Does a CH4 cloud form?
      {
            CCH4 = 1;
            LCH4 = latent_heat(phase_CH4,T);
            LX[4] = LCH4*layer[j-1].XCH4;
            L2X[4]= LX[4]*LCH4/(R*T*T);
      }

      if(PPH3 - SPPH3*SuperSatSelf[2] > SPPH3) //Does a Phosphine cloud form? If so calculate latent heat.
      {
            CPH3 = 1;
            LPH3 = latent_heat(phase_PH3,T);
            LX[5] = LPH3*layer[j-1].XPH3;// LX part of the numerator of DeBoer's Thesis eqn (3.19)
            L2X[5]= LX[5]*LPH3/(R*T*T);//L2X is part of the denominator of DeBoer's Thesis eqn (3.19)
      }

      C = CNH4SH + CH2S + CNH3 + CH2O + CCH4 + CPH3;
      /*  New temperature:  wet adiabat if any of the above clouds condense */
      if (C)
      {
            dT = get_dT(j,layer[j-1].T,layer[j].P,dP,LX,L2X,hereonout,Hydrogen_Curve_Fit_Select);
            T = layer[j].T;
            H = R*T/(layer[j-1].mu*layer[j-1].g);
            alr = 1e5*dT/(dP*H/P);
            wet_adiabatic_lapse_rate=alr;
      }

      if (use_lindal=='Y' && layer[j].P==P_targ && !hereonout) hereonout=1;
/*******************************************************  Updating Mixing Ratios  **************************************************
   Once the latent heat values and coefficients are computed, the cloud density values and saturated abundances of the precursor
  (remaining gas) are calculated. To find the end of this maze of if statements, look for "End phosphine condensate case"
/********************************************************************************************************************************/

/************************************************************************************************************************************
                                                Ammonium Hydrosulfide cloud case
*************************************************************************************************************************************/
      if(CNH4SH)
      {
            dXNH4SH = PH2Sp*PNH3p/(P*(PH2Sp+PNH3p))*(10834.0*dT/T/T - 2.0*dP/P);
            layer[j].DNH4SH = 1e6*AMU_NH4SH*P*P*dXNH4SH/(R*T*dP);
            if(layer[j].DNH4SH > COUNT_CLOUD)
                  layer[j].clouds+= 20L;
            /*  Taken from Hofstadter dissertation (p23)  */
            tdppa = sqrt( SQ(diff_P_base) + 4.0*KNH4SH );
            layer[j].XH2S = 0.5*( (double) tdppa - (double) diff_P_base )/(double) P;
            layer[j].XNH3 = 0.5*( (double) tdppa + (double) diff_P_base )/(double) P;
            if (layer[j].XH2S>layer[j-1].XH2S) layer[j].XH2S=layer[j-1].XH2S;
            if (layer[j].XNH3>layer[j-1].XNH3) layer[j].XNH3=layer[j-1].XNH3;
      }
      else
      {
            dXNH4SH = 0.0;
            layer[j].DNH4SH = ZERO;
            layer[j].XH2S = layer[j-1].XH2S;
            layer[j].XNH3 = layer[j-1].XNH3;
      }
/*************************************************   End Ammonium hydrosulfide      ********************************************/


/************************************************************************************************************************************
                                                Hydrogen Sulfide cloud case
*************************************************************************************************************************************/
      if(CH2S)
      {
            //layer[j].XH2S=(double)SPH2S/(double)P;
            if (sol_cloud)
            {
                  dXH2S = layer[j].XH2S - layer[j-1].XH2S;
                  layer[j].DH2S = ZERO;
                  layer[j].XH2S = (double) SPH2S/(double) P;
            }
            else
            {
                  dXH2S = (LH2S*layer[j-1].XH2S)*dT/(R*T*T) - layer[j-1].XH2S*dP/P;
                  layer[j].DH2S = 1e6*AMU_H2S*P*P*dXH2S/(R*T*dP);
                  layer[j].XH2S = (double) SPH2S/(double) P + ((double)SuperSatSelf[0]*(double)SPH2S)/(double)P;
                  if (layer[j].DH2S > COUNT_CLOUD)
                  {
                        if (T > TRIPLEPT_H2S) layer[j].clouds+=1000L;
                        else layer[j].clouds+=2000L;
                  }
            }
      }
      else
      {
            dXH2S = 0.0;
            layer[j].DH2S = ZERO;
      }

/*************************************************   End Hydrogen Sulfide      ********************************************/


/************************************************************************************************************************************
                            Ammonia condensate case (sub cases of solution cloud, and ammonia ice cloud)
*************************************************************************************************************************************/

      if(CNH3)
      {
            //layer[j].XNH3=(double)SPNH3/(double)P;
            if (sol_cloud)
            {
                  dXNH3 = layer[j].XNH3 - layer[j-1].XNH3;
                  layer[j].DNH3 = ZERO;
                  layer[j].XNH3=(double)SPNH3/(double)P;
            }
            else
            {
                  dXNH3 =(LNH3*layer[j-1].XNH3)*dT/(R*T*T) - layer[j-1].XNH3*dP/P;
                  if(select_ackerman==2 || select_ackerman==3)
                  {
                  dXNH3= -dXNH3;
                  Teff=124; //in Kelvin
                  q_c_nh3_ice=cloud_loss_ackerman_marley(j,Teff, T, P,H, wet_adiabatic_lapse_rate, dry_adiabatic_lapse_rate,layer[j].z, \
                                                layer[j-1].z, layer[j-1].q_c_nh3_ice,layer[j].XNH3,layer[j-1].XNH3, XH2, XHe, XH2S, XNH3, XH2O,\
                                                XCH4, XPH3, dXNH3, frain,Hydrogen_Curve_Fit_Select);
                  layer[j].DNH3 = 1e6*AMU_NH3*P*P*q_c_nh3_ice/(R*T*-dP);
                  layer[j].q_c_nh3_ice=q_c_nh3_ice;
                  layer[j].XNH3 = (double) SPNH3/(double) P + ((double)SPNH3*(double)SuperSatSelf[1])/(double)P;
                  layer[j].first_nh3=1;
                  }
                  else
                  {
                  layer[j].DNH3=1e6*AMU_NH3*P*P*dXNH3/(R*T*dP);
                  layer[j].XNH3 = (double) SPNH3/(double) P + ((double)SPNH3*(double)SuperSatSelf[1])/(double)P;
                  }
                  
                  if (layer[j].DNH3 > COUNT_CLOUD)
                  {
                        if (T > TRIPLEPT_NH3) layer[j].clouds+=100L;
                        else layer[j].clouds+=200L;
                  }
            }
      }
      else
      {
            dXNH3 = 0.0;
            layer[j].DNH3 = ZERO;
      }

/*************************************************   End Ammonia condensate     ********************************************/


/**************************************************************************************************************************
                            Water condensate clould (sub case of solution cloud, and ice cloud)
***************************************************************************************************************************/

      if(CH2O)
      {
            //layer[j].XH2O=(double)SPH2O/(double)P;
            if (sol_cloud)
            {
                  dXH2O = layer[j].XH2O - layer[j-1].XH2O;
                  //Insert Ackerman & Marley Procedure Here
                  
                  if(select_ackerman==1 || select_ackerman==3) //If frain is applied to solution cloud, or both solution cloud and ammonia ice
                  {
                      Teff=124; //Kelvin
                                        
                      q_c=cloud_loss_ackerman_marley(j,Teff, T, P,H, wet_adiabatic_lapse_rate, dry_adiabatic_lapse_rate,layer[j].z, \
                                                layer[j-1].z, layer[j-1].q_c,layer[j].XH2O,layer[j-1].XH2O, XH2, XHe, XH2S, XNH3, XH2O,\
                                                XCH4, XPH3, (-1)*dXH2O*(1-C_sol_NH3), frain,Hydrogen_Curve_Fit_Select);
                      
                      q_c_nh3=cloud_loss_ackerman_marley(j,Teff, T, P,H, wet_adiabatic_lapse_rate, dry_adiabatic_lapse_rate,layer[j].z, \
                                                layer[j-1].z, layer[j-1].q_c_nh3,layer[j].XNH3,layer[j-1].XNH3, XH2, XHe, XH2S, XNH3, XH2O,\
                                                XCH4, XPH3,(-1)*dXNH3*C_sol_NH3, frain,Hydrogen_Curve_Fit_Select);
                      
                      layer[j].DSOL =1e6*((q_c*AMU_H2O+q_c_nh3)*P*P)/(R*T*-dP); //Calulate cloud density g/cm^3 of solution cloud
                      layer[j].DSOL_NH3=1e6*(q_c_nh3*AMU_NH3*P*P)/(R*T*-dP);  // Calculate cloud density g/cm^3 of solution cloud that is actually NH3
                      
                      layer[j].q_c_nh3=q_c_nh3;
                      layer[j].q_c=q_c;

                  if(layer[j-1].first!=1)
                  { 
                        wet_adiabatic_lapse_rate=alr;
                        layer[j].first=1;
                  }      
                  else
                    {
                      LX[3]  = LH2O*layer[j-1].q_c;
                      L2X[3] = LX[3]*LH2O/(R*T*T);
                      dT = get_dT(j,layer[j-1].T,layer[j].P,dP,LX,L2X,hereonout,Hydrogen_Curve_Fit_Select);
                      T = layer[j].T;
                      H = R*T/(layer[j-1].mu*layer[j-1].g);
                      alr = 1e5*dT/(dP*H/P);
                      wet_adiabatic_lapse_rate=alr;
                    }
                    layer[j].XH2O = (double) SPH2O/(double) P + ((double)SuperSatSelf[3]*(double)SPH2O)/(double)P; //adding supersaturation of water as a fraction of saturation
                   
                   layer[j].first=1;//Set the first flag to indicate that this is NOT the first level we've seen a solution cloud.
                  
                  }
                  else
                   {
                    layer[j].DSOL = 1e6*( (1.0-C_sol_NH3)*AMU_H2O*dXH2O + C_sol_NH3*AMU_NH3*dXNH3)*P*P/(R*T*dP);
                    layer[j].DSOL_NH3=1e6*(C_sol_NH3*AMU_NH3*dXNH3*P*P)/(R*T*dP);
                    layer[j].XH2O = (double) SPH2O/(double) P + ((double)SuperSatSelf[3]*(double)SPH2O)/(double)P; //adding supersaturation of water as a fraction of saturation
               
                  }
                  layer[j].DH2O = ZERO;
                  //if we're above the cloud threshould put this "cloud mask address" in layer[j].clouds
                  if (layer[j].DSOL > COUNT_CLOUD)
                        layer[j].clouds+=1L;
            }
            else
            {
                  dXH2O = (LH2O*layer[j-1].XH2O)*dT/(R*T*T) - layer[j-1].XH2O*dP/P;
                  layer[j].DH2O = 1e6*AMU_H2O*P*P*dXH2O/(R*T*dP);   //NOTE!!!! no supersaturation included for ice cloud, if you want it, add it with caution
                  layer[j].XH2O= (double) SPH2O/(double) P;         //     (ie. I wouldn't recommend sticking the same value for solution cloud here
                  if (layer[j].DH2O > COUNT_CLOUD)                  //      as it wouldn't be physically represenative. Different phases will have different
                        layer[j].clouds+=2L;                        //      CCN (Cloud Condensation Nuclei))!!!!!
            }
      }
      else    //If there isn't a water condensate...Don't Make one via numerical error! Keep mole fraction for next level.
      {
            dXH2O = 0.0;
            layer[j].XH2O = layer[j-1].XH2O;
            layer[j].DH2O = ZERO;
      }

/*************************************************      End water condensate case          ********************************************/

/************************************************************************************************************************************
                                                       Methane condensate case
*************************************************************************************************************************************/
      if(CCH4) //If a methane cloud condensate forms. No supersaturation here!
      {
            dXCH4 = (LCH4*layer[j-1].XCH4)*dT/(R*T*T) - layer[j-1].XCH4*dP/P;
            layer[j].XCH4 = (double) SPCH4/(double) P;
            layer[j].DCH4 = 1e6*AMU_CH4*P*P*dXCH4/(R*T*dP);
            if(layer[j].DCH4 > COUNT_CLOUD)
            {
                  if (T > TRIPLEPT_CH4) layer[j].clouds+=10000L;
                  else layer[j].clouds+=20000L;
            }
      }
      else //If it doesn't don't make one vial numerical error! Keep mole fraction for next level.
      {
            dXCH4 = 0.0;
            layer[j].XCH4 = layer[j-1].XCH4;
            layer[j].DCH4 = ZERO;
      }

/***************************************************** End Methane condensate case **************************************************/

/************************************************************************************************************************************
                                                       Phosphine condensate case
*************************************************************************************************************************************/

      if(CPH3)
      {
            dXPH3 = (LPH3*layer[j-1].XPH3)*dT/(R*T*T) - layer[j-1].XPH3*dP/P;
            layer[j].XPH3 = (double) SPPH3/(double) P + ((double)SuperSatSelf[2]*(double)SPPH3)/(double)P; //adding supersaturation of phosphine as a fraction of saturation
            layer[j].DPH3 = 1e6*AMU_PH3*P*P*dXPH3/(R*T*dP);
            if(layer[j].DPH3 > COUNT_CLOUD)
                  layer[j].clouds+=100000L;
      }
      else
      {
            dXPH3 = 0.0;
            layer[j].XPH3 = layer[j-1].XPH3;
            layer[j].DPH3 = ZERO;
      }
/***************************************************** End phosphine condensate case **************************************************/

      //This last piece updates the value for Helium (doesn't condense, so use the previous step)
      layer[j].XHe=layer[j-1].XHe;
      // Check to see if the condensible species are above the "Gone" threshold.
      if(layer[j].XNH3 < GONE) layer[j].XNH3 = ZERO;
      if(layer[j].XH2S < GONE) layer[j].XH2S = ZERO;
      if(layer[j].XCH4 < GONE) layer[j].XCH4 = ZERO;
      if(layer[j].XH2O < GONE) layer[j].XH2O = ZERO;
      if(layer[j].XPH3 < GONE) layer[j].XPH3 = ZERO;

      XHe  = layer[j].XHe;
      XH2S = layer[j].XH2S;
      XNH3 = layer[j].XNH3;
      XH2O = layer[j].XH2O;
      XCH4 = layer[j].XCH4;
      XPH3 = layer[j].XPH3;
      
      //calculated the abundance of H2 as the remainder of the other constituents.
      layer[j].XH2  = (1.0-XH2S-XNH3-XH2O-XCH4-XPH3-XHe-XCO);
      if(Hydrogen_Curve_Fit_Select==666.0)
       {
         
      	vals=get_P_from_python(layer[j].T,layer[j].XH2*layer[j].P,layer[j].XHe*layer[j].P,layer[j].XCH4*layer[j].P,layer[j].XH2O*layer[j].P);
        P_real=float(vals[0])+layer[j].XNH3*P+layer[j].XH2S*P;
        
        
      	layer[j].P_real=P_real;
       }
      XH2  = layer[j].XH2;
      //Calculate the molecular weight (mu in Deboer's thesis chapter 3, not to be confused with cos(theta)) of the "Air"
      layer[j].mu = AMU_H2*XH2 + AMU_He*XHe + AMU_H2S*XH2S + AMU_NH3*XNH3 + AMU_H2O*XH2O + AMU_CH4*XCH4 + AMU_PH3*XPH3;
      layer[j].g = gravity(j);//calculated the acceleration due to gravity cm/sec^2
      return; //Go back to main function in TCM.C returning nothing, but having updated values in layers data structure (see model.h)
}

"""
