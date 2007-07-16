/*This version of layers deals with H2O, NH3 and H2S as a solution system*/

/****************************layers.c**************************************/
#include "model.h"

/*extern struct ATM_LAYER layer[MAXLAYERS];*/
extern struct ATM_LAYER *layer;
//fio extern FILE *lfp;
extern float P_targ, T_targ, P_term, g0, R0, P0, zP0, supersatNH3, supersatH2S;
extern float *TfL, *PfL, lawf[];
extern char use_lindal, line[];
extern int sol_cloud, NH4SH_cloud_base, AutoStep, CrossP0;
extern double XCO,AutoStep_constant;
static int n_lindal_pts, hereonout;
extern float Hydrogen_Curve_Fit_Select;
float gravity(int j);
float specific_heat(int j, float T, float P);
float sat_pressure(char component[], float T);
float solution_cloud(float T, float PNH3, float PH2O, float *SPNH3, float *SPH2O);
float h2s_dissolve(int j, float *SPH2S);
float latent_heat(char component[], float T);
float get_dP(int j, int *eflag, float dz);
float get_dT(int j, float T, float P, float dP, float *LX, float *L2X);
float SuperSatSelf[5];

int init_atm(int n,double XHe,double XH2S,double XNH3,double XH2O,double XCH4,double XPH3,double P_temp,double T_temp,float g0_i,float R0_i, float P0_i,char use_lindal_i, float T_targ_i, float P_targ_i, float P_term_i,int n_lindal_pts_i,float SuperSatSelf1_i,float SuperSatSelf2_i, float SuperSatSelf3_i, float SuperSatSelf4_i, float supersatNH3_i, float supersatH2S_i)    /*set deepest layer*/
{

	  int cntr, i;
      double XH2,temp;
      char c, tp_file[15], tmp[151];
	  FILE *pfp;
	  g0=g0_i; 
	  R0=R0_i; 
	  P0=P0_i;
	  // ***** place abundances in layer data structure *****//	  
	  layer[n].XHe  = XHe;
      layer[n].XH2S = XH2S;
      layer[n].XNH3 = XNH3;
      layer[n].XH2O = XH2O;
      layer[n].XCH4 = XCH4;
      layer[n].XPH3 = XPH3;
      layer[n].XH2  = (1.0-XH2S-XNH3-XH2O-XCH4-XPH3-XHe-XCO);
      layer[n].clouds = 0L;
      layer[n].DH2S = ZERO;
      layer[n].DNH3 = ZERO;
      layer[n].DH2O = ZERO;
      layer[n].DCH4 = ZERO;
      layer[n].DPH3 = ZERO;
      layer[n].DNH4SH = ZERO;
      layer[n].DSOL = ZERO;
      XH2 = layer[n].XH2;
      layer[n].z = 0.0;
	  //*****************************************************//
	  //******* Fix deep Temperature and Pressure ********//
	  layer[n].T=T_temp;
	  layer[n].P=P_temp;
      //*******************************************************************************//
	  
	  //**** Place abundances in data structure ****//
      layer[n].mu   = AMU_H2*XH2 + AMU_He*XHe + AMU_H2S*XH2S + AMU_NH3*XNH3 + AMU_H2O*XH2O + AMU_CH4*XCH4 + AMU_PH3*XPH3;
      layer[n].g = gravity(n); //calculate gravity and put in datastructure
      //*******************************************//
	  
	  //******* Initialize Read in Lindal?, target Temp, target pressure, end pressure level ******//
	  use_lindal=use_lindal_i;
	  T_targ=T_targ_i;
	  P_targ=P_targ_i;
	  P_term=P_term_i;
	  n_lindal_pts=n_lindal_pts_i;
	  //********************************************************************************//
	  
	  //***************** Write values to layer.out file *******************************//  
   //fio   fprintf(lfp,"XH2S=%g XNH3=%g XPH3=%g XH2O=%g XCH4=%g XHE=%g\n",layer[n].XH2S,layer[n].XNH3,layer[n].XPH3,layer[n].XH2O,layer[n].XCH4,layer[n].XHe);
      
	//fio  fprintf(lfp,"T(deep)=%g P(deep)=%g T(targ)=%g P(targ)=%g P(term)=%g\n",layer[n].T,layer[n].P,T_targ,P_targ,P_term);
	  
   //fio   fprintf(lfp,"Supersat:  NH3=%g\tH2S=%g\n",supersatNH3,supersatH2S);
   //fio   fprintf(lfp,"Use Lindal? %c\t\t", use_lindal);
	  //********************************************************************************//
	  
      if (use_lindal == 'Y')
      {
            pfp=fopen("TP.TCM","r");
            TfL = (float *) malloc(n_lindal_pts*sizeof(float));
            PfL = (float *) malloc(n_lindal_pts*sizeof(float));
			
			//** print filename, and number of points in the layers.out file ** //
    //fio        fprintf(lfp,"file used:  %s   points used:  %d\n","TP.SAT",n_lindal_pts);
            //******************************************************************//
			
			if (TfL == NULL || PfL == NULL)
            {
                  printf("Error:  Out of memory...needed %d bytes\n",2*sizeof(float)*n_lindal_pts);
                  return(2);
            }
			
        //    printf("Reading T-P profile from '%s'...","TP.SAT"); //print name of TP file to console
			
			//****** Read in temperature pressure profile ******//
            for(cntr=0;cntr < n_lindal_pts; ++cntr)
                  fscanf(pfp,"%f %f\n",PfL+cntr,TfL+cntr);
			//	  printf("%d",cntr);
            printf("Done.\n"); //print Done to the console
			//**************************************************//
      
		    P_targ=PfL[0]; //set P_target to first element
            T_targ=TfL[0]; // set T_taget to first element
            P_term=PfL[n_lindal_pts-1]; //set final pressure to final element
      }
	  
	  supersatNH3=supersatNH3_i;
	  supersatH2S=supersatH2S_i;
      SuperSatSelf[0]=SuperSatSelf1_i;
	  SuperSatSelf[1]=SuperSatSelf2_i;
	  SuperSatSelf[2]=SuperSatSelf3_i;
	  SuperSatSelf[3]=SuperSatSelf4_i;
	  
	    if ((SuperSatSelf[0]+SuperSatSelf[1]+SuperSatSelf[2]+SuperSatSelf[3])>0.0)
      {
      	printf("Self supersaturations:\t"); //print values to console
         for (i=0;i<4;++i)
         	if (SuperSatSelf[i]>0.0)
            	printf("%g  ",SuperSatSelf[i]);
      }
      else
      {
      	printf("No self superaturation.\n"); //print no values to console
      }
	  
	  /* ...flag value */
	  sol_cloud=1;
	  /********************/
	  
      printf("\nDo %s consider solution (windex) cloud formation.\n",sol_cloud?"":"not"); //print sol_cloud to console
      
//fio	  fputc('\n',lfp); //write carriage return to layers.out file
      
	  fclose(pfp);  // close planet file or TP profile 
	  
      CrossP0=0;
	  
      return (1);
}

void new_layer(int j, float dz, int *eflag)
{
      int CNH4SH=0, CH2S=0, CNH3=0, CH2O=0, CCH4=0, CPH3=0, C=0, C_sol_zero=0;
      float LX[6]={0.0,0.0,0.0,0.0,0.0,0.0}, L2X[6]={0.0,0.0,0.0,0.0,0.0,0.0};
      float P, T, dT, dP, PNH3p, PH2Sp, freeze;
      float PH2, PHe, PH2S, PNH3, PH2O, PCH4, PPH3;
      double XH2, XHe, XH2S, XNH3, XH2O, XCH4, XPH3;
      float SPH2O=1E6, SPCH4=1E6, SPH2S=1E6, SPNH3=1E6, SPPH3=1E6, KNH4SH;
      float LH2S, LNH3, LNH4SH, LH2O, LCH4, LPH3;
      float dXH2S, dXNH3, dXNH4SH, dXH2O, dXCH4, dXPH3;
      float H, alr, tdppa, C_sol_NH3, C_sol_H2S;
      static float diff_P_base;
      char phase_H2S[21], phase_NH3[21], phase_PH3[21];
      char phase_H2O[21], phase_CH4[21];
      FILE *alrfp;

      layer[j].clouds=0L;
      if (j==1) hereonout=0;
      /*  Get new P,dP and T,dT values  */
      dP = get_dP(j,eflag,dz);
      dT = get_dT(j,layer[j-1].T,layer[j].P,dP,LX,L2X); /*dry adiabat*/
      P  = layer[j].P;
      T  = layer[j].T;
 //fio     fprintf(lfp,"New pressure:  %g + %g = %g\n",layer[j-1].P,dP,P);
      PH2  = layer[j-1].XH2*P;
      PHe  = layer[j-1].XHe*P;
      PH2S = layer[j-1].XH2S*P;
      PNH3 = layer[j-1].XNH3*P;
      PPH3 = layer[j-1].XPH3*P;
      PH2O = layer[j-1].XH2O*P;
      PCH4 = layer[j-1].XCH4*P;
      PNH3p= PNH3 - supersatNH3*P;
      PH2Sp= PH2S - supersatH2S*P;
      H = R*T/(layer[j-1].mu*layer[j-1].g);
      alr = 1e5*dT/(dP*H/P);
 //fio     fprintf(lfp,"dry adiabatic lapse rate %g K/km\n",alr);
 //fio     fprintf(lfp,"dry new T:  %g + %g = %g\n",layer[j-1].T,dT,T);

      /* For the solution clouds see:  Weidenschilling and Lewis 1973, Icarus 20:465.
            Lewis 1969, Icarus 10:365. Briggs and Sackett 1989, Icarus 80:77.
            Atreya and Romani 1985, Recent Advances in Planetary Meteorology.
            Atreya 1986, Atmospheres and Ionospheres of the Outer Planets. */
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
//fio                  fprintf(lfp,"Reached solution cloud top:  %g concentrate.\n",C_sol_NH3);
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
   //fio                           fprintf(lfp,"Solution cloud H2S content:  %g\n",C_sol_H2S);
                        }
/****************************************************/
                  }
                  LH2O = C_sol_NH3*(4949.75+2022.11*(SQ(C_sol_NH3)-2.0*C_sol_NH3));
                  LH2O+= (1.0-C_sol_NH3)*(5540.48+2022.11*SQ(C_sol_NH3));
                  LH2O*= R/(1.0-C_sol_NH3);
                  LX[3]  = LH2O*layer[j-1].XH2O;
                  L2X[3] = LX[3]*LH2O/(R*T*T);
    //fio              fprintf(lfp,"Solution cloud:  %g concentrate.  LSOL=%g\n",C_sol_NH3,LH2O);
    //fio              fprintf(lfp,"SPH2O=%g\tPH2O=%g\tSPNH3=%g\tPNH3=%g\n",SPH2O,PH2O,SPNH3,PNH3);
            }
      }

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

            if(PH2O - SuperSatSelf[3] > SPH2O)
            {
                  CH2O = 1;
                  LH2O = latent_heat(phase_H2O,T);
     //fio             fprintf(lfp,"%s:  ",phase_H2O);
      //fio            fprintf(lfp,"PH2O=%g  SPH2O=%g  LH2O=%g\n",PH2O,SPH2O,LH2O);
                  LX[3] = LH2O*layer[j-1].XH2O;
                  L2X[3]= LX[3]*LH2O/(R*T*T);
            }

            if(PH2S - SuperSatSelf[0] > SPH2S)
            {
                  CH2S = 1;
                  LH2S = latent_heat(phase_H2S,T);
       //fio           fprintf(lfp,"%s:  ",phase_H2S);
       //fio           fprintf(lfp,"PH2S=%g  SPH2S=%g  LH2S=%g\n",PH2S,SPH2S,LH2S);
                  LX[1] = LH2S*layer[j-1].XH2S;
                  L2X[1]= LX[1]*LH2S/(R*T*T);
            }

            if(PNH3 - SuperSatSelf[1] > SPNH3)
            {
                  CNH3 = 1;
                  LNH3 = latent_heat(phase_NH3,T);
        //fio          fprintf(lfp,"%s:  ",phase_NH3);
        //fio          fprintf(lfp,"PNH3=%g  SPNH3=%g  LNH3=%g\n",PNH3,SPNH3,LNH3);
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

      if(PH2Sp*PNH3p > KNH4SH)
      {
            if ( !NH4SH_cloud_base )
            {
                  NH4SH_cloud_base=1;
                  diff_P_base = PNH3p - PH2Sp;
            }
            CNH4SH = 1;
            LNH4SH = 9.312E11; /* Atreya book.  1.6E12 Briggs and Sackett */
     //fio       fprintf(lfp,"NH4SH cloud:  ");
     //fio       fprintf(lfp,"PNH3PH2S=%g  KNH4SH=%g  LNH4SH=%g",PH2Sp*PNH3p,KNH4SH,LNH4SH);
     //fio       fprintf(lfp,"  XNH3=%g  XNH3p=%g  XH2S=%g  XH2Sp=%g\n",layer[j-1].XNH3,PNH3p/P,layer[j-1].P,PH2Sp);
            LX[0] = 2.0*LNH4SH*PH2Sp*PNH3p/(P*(PH2Sp+PNH3p));
            L2X[0]= LX[0]*5417.0/(T*T);
      }

      if(PCH4 > SPCH4)
      {
            CCH4 = 1;
            LCH4 = latent_heat(phase_CH4,T);
     //fio       fprintf(lfp,"%s:  ",phase_CH4);
     //fio       fprintf(lfp,"PCH4=%g  SPCH4=%g  LCH4=%g\n",PCH4,SPCH4,LCH4);
            LX[4] = LCH4*layer[j-1].XCH4;
            L2X[4]= LX[4]*LCH4/(R*T*T);
      }

      if(PPH3 - SuperSatSelf[2] > SPPH3)
      {
            CPH3 = 1;
            LPH3 = latent_heat(phase_PH3,T);
     //fio       fprintf(lfp,"%s:  ",phase_PH3);
      //fio      fprintf(lfp,"PPH3=%g  SPPH3=%g  LPH3=%g\n",PPH3,SPPH3,LPH3);
            LX[5] = LPH3*layer[j-1].XPH3;
            L2X[5]= LX[5]*LPH3/(R*T*T);
      }

      C = CNH4SH + CH2S + CNH3 + CH2O + CCH4 + CPH3;
      /*  New temperature:  wet adiabat if C != 0 */
      if (C)
      {
            dT = get_dT(j,layer[j-1].T,layer[j].P,dP,LX,L2X);
            T = layer[j].T;
            H = R*T/(layer[j-1].mu*layer[j-1].g);
            alr = 1e5*dT/(dP*H/P);
    //fio        fprintf(lfp,"wet adiabatic lapse rate %g K/km\n",alr);
    //fio        fprintf(lfp,"wet new T:  %g + %g = %g\n",layer[j-1].T,dT,T);
      }

      if (use_lindal=='Y' && layer[j].P==P_targ && !hereonout) hereonout=1;
      /*  Update mixing ratios  */
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
    //fio  fprintf(lfp,"dXNH4SH = %g  ",dXNH4SH);

      if(CH2S)
      {
            layer[j].XH2S = (double) SPH2S/(double) P + SuperSatSelf[0];
            if (sol_cloud)
            {
                  dXH2S = layer[j].XH2S - layer[j-1].XH2S;
                  layer[j].DH2S = ZERO;
            }
            else
            {
                  dXH2S = (LH2S*layer[j-1].XH2S)*dT/(R*T*T) - layer[j-1].XH2S*dP/P;
                  layer[j].DH2S = 1e6*AMU_H2S*P*P*dXH2S/(R*T*dP);
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
      //fio fprintf(lfp,"dXH2S = %g  ",dXH2S);

      if(CNH3)
      {
            layer[j].XNH3 = (double) SPNH3/(double) P + SuperSatSelf[1];
            if (sol_cloud)
            {
                  dXNH3 = layer[j].XNH3 - layer[j-1].XNH3;
                  layer[j].DNH3 = ZERO;
            }
            else
            {
                  dXNH3 = (LNH3*layer[j-1].XNH3)*dT/(R*T*T) - layer[j-1].XNH3*dP/P;
                  layer[j].DNH3 = 1e6*AMU_NH3*P*P*dXNH3/(R*T*dP);
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
     //fio fprintf(lfp,"dXNH3 = %g  ",dXNH3);

      if(CH2O)
      {
            layer[j].XH2O = (double) SPH2O/(double) P + SuperSatSelf[3];
            if (sol_cloud)
            {
                  dXH2O = layer[j].XH2O - layer[j-1].XH2O;
                  layer[j].DSOL = 1e6*( (1.0-C_sol_NH3)*AMU_H2O*dXH2O + C_sol_NH3*AMU_NH3*dXNH3)*P*P/(R*T*dP);
                  layer[j].DH2O = ZERO;
                  if (layer[j].DSOL > COUNT_CLOUD)
                        layer[j].clouds+=1L;
            }
            else
            {
                  dXH2O = (LH2O*layer[j-1].XH2O)*dT/(R*T*T) - layer[j-1].XH2O*dP/P;
                  layer[j].DH2O = 1e6*AMU_H2O*P*P*dXH2O/(R*T*dP);
                  if (layer[j].DH2O > COUNT_CLOUD)
                        layer[j].clouds+=2L;
            }
      }
      else
      {
            dXH2O = 0.0;
            layer[j].XH2O = layer[j-1].XH2O;
            layer[j].DH2O = ZERO;
      }
    //fio  fprintf(lfp,"dXH2O = %g  ",dXH2O);

      if(CCH4)
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
      else
      {
            dXCH4 = 0.0;
            layer[j].XCH4 = layer[j-1].XCH4;
            layer[j].DCH4 = ZERO;
      }
    //fio  fprintf(lfp,"dXCH4 = %g  ",dXCH4);

      if(CPH3)
      {
            dXPH3 = (LPH3*layer[j-1].XPH3)*dT/(R*T*T) - layer[j-1].XPH3*dP/P;
            layer[j].XPH3 = (double) SPPH3/(double) P + SuperSatSelf[2];
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
     //fio fprintf(lfp,"dXPH3 = %g\n",dXPH3);

      layer[j].XHe=layer[j-1].XHe;

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
      layer[j].XH2  = (1.0-XH2S-XNH3-XH2O-XCH4-XPH3-XHe-XCO);
      XH2  = layer[j].XH2;

      layer[j].mu = AMU_H2*XH2 + AMU_He*XHe + AMU_H2S*XH2S + AMU_NH3*XNH3 + AMU_H2O*XH2O + AMU_CH4*XCH4 + AMU_PH3*XPH3;
      layer[j].g = gravity(j);
   //fio   fprintf(lfp,"cloud structure %06ld",layer[j].clouds);
      return;
}

float get_dP(int j, int *eflag, float dz)
{
      float P, T, H, dP;

      if (AutoStep)
      { 
            /*dz = log10(layer[j-1].P + 1.6); */
           /* dz = log10(layer[j-1].P + 5.0); */
		   printf("%f",AutoStep_constant);
		   dz=log10(layer[j-1].P + AutoStep_constant );  //what dave actually uses for Priscilla's stuff
            if (dz > 2.0) dz=2.0;
      }

  //fio    fprintf(lfp,"\tdz = %f km\n",dz);
      H = R*layer[j-1].T/(layer[j-1].mu*layer[j-1].g);
      dP = -1.0*(layer[j-1].P/H)*(1.0E5*dz);
      P  = layer[j-1].P + dP;

      if (P <= P_term)
      {
            *eflag = 99;
            P = P_term;
            dP = P - layer[j-1].P;
      }
      else if ( P <= P_targ && *eflag != 97 )
      {
            *eflag = 98;
            P = P_targ;
            dP = P - layer[j-1].P;
      }

      layer[j].P = P;
      layer[j].z = layer[j-1].z + dz;
      if (P<P0 && !CrossP0) 
      {
           zP0 = layer[j].z;
           CrossP0 = 1;
      }

      return (dP);
}


float get_dT(int j, float T, float P, float dP, float *LX, float *L2X)
{
      int i;
      float m, b, Cp, dT_num, dT_den, dT; 

      if (hereonout)   /*linear interpolation from Lindal's points*/
      {
            for (i=0; P <= PfL[i]; ++i)  ;
            m = (TfL[i] - TfL[i-1])/(PfL[i] - PfL[i-1]);
            b = TfL[i] - m*PfL[i];
            layer[j].T = m*P + b;
            dT = layer[j].T - layer[j-1].T;
      }
      else
      {
            Cp = specific_heat(j, T, P);
            dT_num = (R*T + LX[0]+LX[1]+LX[2]+LX[3]+LX[4]+LX[5])*dP;
            dT_den = P*(Cp + L2X[0]+L2X[1]+L2X[2]+L2X[3]+L2X[4]+L2X[5]);
            dT = dT_num/dT_den;
            layer[j].T = layer[j-1].T + dT;
      }
      return (dT);
}

float gravity(int j)
{
      float gamma, b;

      gamma = 2.0*log(layer[j].P/P0)*(R*layer[j].T)/(R0*layer[j].mu);
      b = g0 + gamma;
      return  0.5*( b + sqrt( SQ(b) - SQ(gamma) ) );
}
