/**************************** TCM.C **************************************/
/*  This model requires the data files:                                  */
/*                                       inputsoln.dat                   */
/*************************************************************************/
#include "model.h"
//struct ATM_LAYER layer[MAXLAYERS];
struct ATM_LAYER *layer;
//fio FILE *lfp;
char use_lindal;
char line[151];
int NH4SH_cloud_base=0;
float P_targ, T_targ, P_term, g0, R0, P0, supersatNH3, supersatH2S,use_lindal_in;
float *TfL, *PfL, lawf[3],zP0;
double *funoutput;
double XCO;
int AutoStep,CrossP0;
double AutoStep_constant;
int use_dz,jj;
/*set 0 to turn off, 1 to turn on*/
#define USE_SOL_CLOUD  1      /*consider its formation?*/
double P_temp,T_temp;
float z_offset;
double **tcm;
long *cldFlags;
void new_layer(int j, float dz, int *eflag, float dP_init, float dP_fine, float P_fine_start, float P_fine_stop, float frain,float select_ackerman);
void new_layer_original(int j, float dz, int *eflag, float dP_init, float dP_fine, float P_fine_start, float P_fine_stop);
float specific_heat(int j, float T, float P,float Cp_in);
int init_atm(int n,double XHe,double XH2S,double XNH3,double XH2O,double XCH4,double XPH3,double P_temp,double T_temp,float g0_i,float R0_i, float P0_i,char use_lindal_i, float T_targ_i, float P_targ_i, float P_term_i,int n_lindal_pts_i,float SuperSatSelf1_i,float SuperSatSelf2_i, float SuperSatSelf3_i, float SuperSatSelf4_i,float supersatNH3_i,float supersatH2S_i);
int intoTheVoid( float dz,
                  double XHe_i,
                  double XH2S_i,
                  double XNH3_i,
                  double XH2O_i,
                  double XCH4_i,
                  double XPH3_i,
                  double XCO,
                  float P_temp,
                  float T_temp,
                  float g0_i,
                  float R0_i,
                  float P0_i,
                  float T_targ_i,
                  float P_targ_i,
                  float P_term_i,
                  float use_lindal_in,
                  int n_lindal_pts_i,
                  float SuperSatSelf1_i,                                             
                  float SuperSatSelf2_i,
                  float SuperSatSelf3_i,
                  float SuperSatSelf4_i,
                  float supersatNH3_i,
                  float supersatH2S_i,
                  int AutoStep_constant,
                  float Hydrogen_Curve_Fit_Select,
                  float dP_init,
                  float dP_fine,
                  float P_fine_start,
                  float P_fine_stop,
                  int use_dz,
                  float frain,
                  float select_ackerman);
double getFloatValues(int i);
long getCloudFlags(int i);
void init_soln_cloud(int mode);
void dt_tag(FILE *fp);
int sol_cloud=USE_SOL_CLOUD;

// Function declarations.
// -----------------------------------------------------------------
float Hydrogen_Curve_Fit_Select;


/*********************************************************************
 void mexFunction(): Main function which returns values to Matlab (tm)
              
            Input:
                -->int nlhs: number of variables on the left hand side (output)
                -->int nrhs: number of variables on the right hand side (input)
                -->mxArray prhs[]: variables from matlab to the mexFunction.
            Output:
                <-- mxArray plhs[]: variables returned to Maltab (tm)   

*********************************************************************/
int intoTheVoid( float dz, 
      		  double XHe_i,
      		  double XH2S_i,
      		  double XNH3_i,
                  double XH2O_i,
                  double XCH4_i,
                  double XPH3_i,
                  double XCO,
                  float P_temp,
                  float T_temp,
                  float g0_i,
                  float R0_i,
                  float P0_i,
                  float T_targ_i,
                  float P_targ_i,
                  float P_term_i,
                  float use_lindal_in,
                  int n_lindal_pts_i,
                  float SuperSatSelf1_i,
                  float SuperSatSelf2_i,
                  float SuperSatSelf3_i,
                  float SuperSatSelf4_i,
                  float supersatNH3_i,
                  float supersatH2S_i,
                  int AutoStep_constant,
                  float Hydrogen_Curve_Fit_Select,
                  float dP_init,
                  float dP_fine,
                  float P_fine_start,
                  float P_fine_stop,
                  int use_dz,
                  float frain,
                  float select_ackerman)
{
      int bottom, top, j, eflag=0, jcntr=0, found,cross_my_P0=0;
      float P, T, T_err=100.0, f=1.0, T_bright, T_ob, W,Pr;
      float refr, refr2, dawf;
      char outfile[15],use_lindal_i='Y';
      //FILE *ofp1;
	  

      printf("Starting DeBoer/Steffes/Karpowicz Thermo-chemical Model\n");		
      
       //Allocate memory for layer data structure.
      layer = (struct ATM_LAYER *) calloc(MAXLAYERS, sizeof(struct ATM_LAYER));
      
      if (layer==NULL) {printf("Insufficient memory.\n");return 0;}
      
      //initialize the bottom layer with values from Matlab

      found = init_atm(0, XHe_i,XH2S_i,XNH3_i,XH2O_i,XCH4_i,XPH3_i,P_temp,T_temp,g0_i, R0_i, \
                      P0_i,use_lindal_i, T_targ_i, P_targ_i, P_term_i,n_lindal_pts_i,\
                      SuperSatSelf1_i,SuperSatSelf2_i,SuperSatSelf3_i,SuperSatSelf4_i,\
                      supersatNH3_i,supersatH2S_i);
	  
      if (USE_SOL_CLOUD)
            init_soln_cloud(1);
      printf("giggity.");
      if (dz==0.0) {printf("Using AutoStep.\n"); AutoStep=1;}
      else AutoStep=0;
	  
      //while error in temperature is greater than Tlimit, or if you've gone on long enough
      //  without converging on a solution ...  
      
      while (fabs(T_err) > TLIMIT && jcntr < MAXTRIES)
      {
            eflag = 0; //reset flag
            jj=0;     //reset jj
            ++jcntr;  //increment jcntr
           		
            printf("iteration: %d      \n",jcntr);
            		
            for(j=1;eflag!=99 && j<MAXLAYERS;++j) //go through the layers          
            {                                               
              
              new_layer(j,dz,&eflag, dP_init, dP_fine, P_fine_start, P_fine_stop, frain,select_ackerman);
              P = layer[j].P;
              //printf("error flag Pressure, level,T_targ %d %f %d %f\n",eflag,P,j,T_targ);
                  
                  if (eflag == 98)  /* check target temperature */
                  {                                         /* eflag       P      */
                        eflag = 97;                         /*  97      < P_targ  */
                        T = layer[j].T;                     /*  99      <=P_term  */
                        T_err = 100.0*(T - T_targ)/T_targ;  /*  98      <=P_targ  */
                        //printf("Terr %f %f \n",T_err,TLIMIT);
                      if (fabs(T_err)>TLIMIT)
                        {
                              layer[0].T = layer[0].T - (T - T_targ)*2.0;
                              P = layer[0].P;
                              eflag = 99;
                              sol_cloud = USE_SOL_CLOUD;
                              NH4SH_cloud_base = 0;
                        }
                     
                  }
            }
			
      }
  
      if (P > P_term)
      {
             printf("Abnormal termination. P= %f bigger than Pterm=%f                       \n",P,P_term);
             printf("jj val %d",jj);
             return 0;
      }

      top = j-1;

      if (USE_SOL_CLOUD)
            init_soln_cloud(0); //free solution cloud memory allocation
  
      if (use_lindal == 'Y' || use_lindal == 'y')
      {
            free(PfL); //free lindal profile memory allocation
            free(TfL);
      }
      printf("Sending Data your way via Matlab mex! \n \n");
      
      //ofp1=fopen("tcm.out","w");
      //funoutput = (double *) malloc(top*sizeof(double));
      /***********************************************
       Output data to Matlab (tm)
      **********************************************/
      /*1 2  3   4   5    6    7    8     9      10    11   12   13   14   15   16   17  18 19 20       21        */
     /*P T XH2 XHe XH2S XNH3 XH2O XCH4 XPH3  clouds DNH4SH DH2S DNH3 DH2O DCH4 DPH3 DSOL g mu refr_w/o refr_ w/  */
     /* P_out = mxGetPr(plhs[0]); //1
     	T_out=mxGetPr(plhs[1]); //2
     	XH2_out=mxGetPr(plhs[2]); //3
    	 XHe_out=mxGetPr(plhs[3]); //4
     	XH2S_out=mxGetPr(plhs[4]); //5
    	 XNH3_out=mxGetPr(plhs[5]); //6
     	XH2O_out=mxGetPr(plhs[6]); //7
     	XCH4_out=mxGetPr(plhs[7]); //8 
     	XPH3_out=mxGetPr(plhs[8]); //9
     	clouds_out=mxGetPr(plhs[9]); //10
     	DNH4SH_out=mxGetPr(plhs[10]); //11
     	DH2S_out=mxGetPr(plhs[11]); //12
     	DNH3_out=mxGetPr(plhs[12]); //13
     	DH2O_out=mxGetPr(plhs[13]); //14
     	DCH4_out=mxGetPr(plhs[14]); //15
    	 DPH3_out=mxGetPr(plhs[15]); //16
     	DSOL_out=mxGetPr(plhs[16]); //17
     	g_out=mxGetPr(plhs[17]); //18
     	mu_out=mxGetPr(plhs[18]); //19
     	refr_w_o_out=mxGetPr(plhs[19]); //20
     	refr_w_out=mxGetPr(plhs[20]); //21
     	z_out=mxGetPr(plhs[21]); //22
     	DSOL_NH3_out=mxGetPr(plhs[22]);//23
     	P_real=mxGetPr(plhs[23]);
     */
     for (j =0 ;j<=top;++j)
      {     
         if(layer[j].P<P0 && !cross_my_P0)
            {
                  z_offset=layer[j].z;
                  cross_my_P0=1;
            }
      }

/*******************Print sutff to output file, just for fun, and fill the matlab arrays. ***********/
//      for(j = 0;j<=top;++j)
//      {
//            
//            T = layer[j].T;
//            P = layer[j].P;
//            Pr= layer[j].P_real;
//            P_real[j]=Pr;
//            P_out[j]=P;
//            T_out[j]=T;
//            XH2_out[j]=layer[j].XH2;
//            XHe_out[j]=layer[j].XHe;
//            XH2S_out[j]=layer[j].XH2S;
//            XNH3_out[j]=layer[j].XNH3;
//            XH2O_out[j]=layer[j].XH2O;
//            XCH4_out[j]=layer[j].XCH4;
//            XPH3_out[j]=layer[j].XPH3;
//            clouds_out[j]=layer[j].clouds;
//            DNH4SH_out[j]=layer[j].DNH4SH;
//            DH2S_out[j]=layer[j].DH2S;
//            DNH3_out[j]=layer[j].DNH3;
//            DH2O_out[j]=layer[j].DH2O;
//            DCH4_out[j]=layer[j].DCH4;
//            DPH3_out[j]=layer[j].DPH3;
//            DSOL_out[j]=layer[j].DSOL;
//            g_out[j]=layer[j].g;
//            mu_out[j]=layer[j].mu;
//            DSOL_NH3_out[j]=layer[j].DSOL_NH3;
//            refr = layer[j].XH2*REFR_H2  +layer[j].XHe*REFR_He  +layer[j].XH2S*REFR_H2S+layer[j].XPH3*REFR_PH3;
//            refr+= layer[j].XNH3*REFR_NH3+layer[j].XCH4*REFR_CH4+layer[j].XH2O*REFR_H2O(T);
//            refr*= ( P*(293.0/T) );
//			
//            refr2= (layer[j].XH2*REFR_H2  +layer[j].XHe*REFR_He)*P*(293.0/T);
//			refr_w_o_out[j]=refr;
//			refr_w_out[j]=refr2;
//			z_out[j]=layer[j].z - z_offset;
 //            /*1 2      3   4   5    6    7    8     9      10    11   12   13   14   15   16   17  18 19 20       21        */
 //            /*P T dz   XH2 XHe XH2S XNH3 XH2O XCH4 XPH3  clouds DNH4SH DH2S DNH3 DH2O DCH4 DPH3 DSOL g mu refr_w/o refr_ w/  */
 //            fprintf(ofp1,"%7.3f\t%6.3f\t%7.3f\t%10.3e\t%10.3e\t",P,T,z_out[j],layer[j].XH2,layer[j].XHe);
 //           
 //            fprintf(ofp1,"%10.3e\t%10.3e\t%10.3e\t%10.3e\t%10.3e\t",layer[j].XH2S,layer[j].XNH3,layer[j].XH2O,layer[j].XCH4,layer[j].XPH3);
 //           
 //            fprintf(ofp1,"%06ld\t%10.3e\t%10.3e\t%10.3e\t%10.3e\t%10.3e\t%10.3e\t%10.3e\t",layer[j].clouds,layer[j].DNH4SH,layer[j].DH2S,\
 //                   layer[j].DNH3,layer[j].DH2O,layer[j].DCH4,layer[j].DPH3,layer[j].DSOL);
 //           
 //            fprintf(ofp1,"%g\t%g\t%g\t%g\n",layer[j].g,layer[j].mu,refr,refr2);
 //     }
 //     
 //     fclose(ofp1);//close file pointer
 
/***********************************************************************/
      printf("End of Thermo-chemical model, return to Matlab.\n\n");
     
      return (top);		  
}


long getCloudFlags(int i)
{
	return layer[i].clouds;
}

double getFloatValues(int i,int j)
{

		if(i==0) return layer[j].P;
		if(i==1) return layer[j].T;
		if(i==2) return layer[j].XH2;
		if(i==3) return layer[j].XHe;
		if(i==4) return layer[j].XH2S;
		if(i==5) return layer[j].XNH3;
		if(i==6)return layer[j].XH2O;
		if(i==7)return layer[j].XCH4;
		if(i==8)return layer[j].XPH3;
		if(i==9)return layer[j].DNH4SH;
		if(i==10)return layer[j].DH2S;
		if(i==11)return layer[j].DNH3;
		if(i==12)return layer[j].DH2O;
		if(i==13)return layer[j].DCH4;
		if(i==14)return layer[j].DPH3;
		if(i==16)return layer[j].DSOL;
		if(i==17)return layer[j].g;
		if(i==18)return layer[j].mu;
		if(i==19)return layer[j].z-z_offset;
		if(i==20)return layer[j].DSOL_NH3;
		if(i==21)return layer[j].P_real;
}



/**********Misc functions, Not really used anymore *****/
int y_n_answer(void)
{
	char lp;

SNF:
      lp = getc(stdin);;
      if (lp == 'n' || lp == 'N')
      {
            printf("No\n");
		return(0);
      }
      else if (lp == 'y' || lp == 'Y')
      {
            printf("Yes\n");
		return(1);
      }
      else 
            goto SNF;
}
/*int getline(FILE *fp, char *s, int lim)
{
      char c;
      int i;
	  
      for (i=0;i<lim-2 && (c=getc(fp))!=EOF && c!='\n';++i)
           s[i] = c;
      s[i] = '\0';

      return (i);
}*/
/**************That's All Folks!********************************/
