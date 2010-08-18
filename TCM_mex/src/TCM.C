/**************************** TCM.C **************************************/
/*  This model requires the data files:                                  */
/*                                       inputsoln.dat                   */
/*************************************************************************/
#include "model.h"
#include "mex.h"
/*struct ATM_LAYER layer[MAXLAYERS];*/
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
void new_layer(int j, float dz, int *eflag, float dP_init, float dP_fine, float P_fine_start, float P_fine_stop, float frain,float select_ackerman);
void new_layer_original(int j, float dz, int *eflag, float dP_init, float dP_fine, float P_fine_start, float P_fine_stop);
float specific_heat(int j, float T, float P);
int init_atm(int n,double XHe,double XH2S,double XNH3,double XH2O,double XCH4,double XPH3,double P_temp,double T_temp,float g0_i,float R0_i, float P0_i,char use_lindal_i, float T_targ_i, float P_targ_i, float P_term_i,int n_lindal_pts_i,float SuperSatSelf1_i,float SuperSatSelf2_i, float SuperSatSelf3_i, float SuperSatSelf4_i,float supersatNH3_i,float supersatH2S_i);
void init_soln_cloud(int mode);
void dt_tag(FILE *fp);
int sol_cloud=USE_SOL_CLOUD;
//********mex stuff********//
extern void _main();
// Function declarations.
// -----------------------------------------------------------------
char getMatlabCharacter(const mxArray* ptr);
double  getMatlabScalar    (const mxArray* ptr);
int getMatlabInt (const mxArray* ptr);
double& createMatlabScalar (mxArray*& ptr);
const int numInputArgs  = 33;
const int numOutputArgs = 24;
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
void mexFunction(int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{
      int bottom, top, j, eflag=0, jcntr=0, found,cross_my_P0=0;
      float P, T, T_err=100.0, f=1.0, T_bright, T_ob, W,Pr;
      float refr, refr2, dawf,z_offset;
      char outfile[15],use_lindal_i='Y';
      double *P_out,*T_out,*XH2_out,*XHe_out,*XH2S_out,*XNH3_out,*DH2S_out,*XH2O_out,\
             *XCH4_out,*XPH3_out,*clouds_out,*DNH4SH_out,*DNH3_out,*DH2O_out,*DCH4_out,\
             *DPH3_out,*DSOL_out,*g_out,*mu_out,*refr_w_o_out,*refr_w_out,*z_out,*DSOL_NH3_out,\
             *P_real;
	  
      FILE *ofp1;
	  
      //*************Error messages for mex (matlab)*****************//
      if (nrhs != numInputArgs) mexErrMsgTxt("Incorrect number of input arguments");
      if (nlhs != numOutputArgs) mexErrMsgTxt("Incorrect number of output arguments");
      //***************************************************//

      //*************** Get values from matlab **********//
      float dz  = getMatlabScalar(prhs[0]);
      double XHe_i=getMatlabScalar(prhs[1]);
      double XH2S_i=getMatlabScalar(prhs[2]);
      double XNH3_i=getMatlabScalar(prhs[3]);
      double XH2O_i=getMatlabScalar(prhs[4]);
      double XCH4_i=getMatlabScalar(prhs[5]);
      double XPH3_i=getMatlabScalar(prhs[6]);
      double XCO=getMatlabScalar(prhs[7]);
      float P_temp=getMatlabScalar(prhs[8]);
      float T_temp=getMatlabScalar(prhs[9]);
      float g0_i=getMatlabScalar(prhs[10]);
      float R0_i=getMatlabScalar(prhs[11]);
      float P0_i=getMatlabScalar(prhs[12]);
      float T_targ_i=getMatlabScalar(prhs[13]);
      float P_targ_i=getMatlabScalar(prhs[14]);
      float P_term_i=getMatlabScalar(prhs[15]);
      float use_lindal_in=getMatlabScalar(prhs[16]);
      int n_lindal_pts_i=getMatlabInt(prhs[17]);
      float SuperSatSelf1_i=getMatlabScalar(prhs[18]);
      float SuperSatSelf2_i=getMatlabScalar(prhs[19]);
      float SuperSatSelf3_i=getMatlabScalar(prhs[20]);
      float SuperSatSelf4_i=getMatlabScalar(prhs[21]);
      float supersatNH3_i=getMatlabScalar(prhs[22]);
      float supersatH2S_i=getMatlabScalar(prhs[23]);
      AutoStep_constant=getMatlabScalar(prhs[24]);
      Hydrogen_Curve_Fit_Select=getMatlabScalar(prhs[25]);
      float dP_init=getMatlabScalar(prhs[26]);
      float dP_fine=getMatlabScalar(prhs[27]);
      float P_fine_start=getMatlabScalar(prhs[28]);
      float P_fine_stop=getMatlabScalar(prhs[29]);
      use_dz=getMatlabScalar(prhs[30]);
      float frain=getMatlabScalar(prhs[31]);
      float select_ackerman=getMatlabScalar(prhs[32]);
      //************************************************//
      
      printf("Starting DeBoer/Steffes/Karpowicz Thermo-chemical Model\n");		
      
       //Allocate memory for layer data structure.
      layer = (struct ATM_LAYER *) calloc(MAXLAYERS, sizeof(struct ATM_LAYER));
      
      if (layer==NULL) {printf("Insufficient memory.\n");return;}
      
      //initialize the bottom layer with values from Matlab

      found = init_atm(0, XHe_i,XH2S_i,XNH3_i,XH2O_i,XCH4_i,XPH3_i,P_temp,T_temp,g0_i, R0_i, \
                      P0_i,use_lindal_i, T_targ_i, P_targ_i, P_term_i,n_lindal_pts_i,\
                      SuperSatSelf1_i,SuperSatSelf2_i,SuperSatSelf3_i,SuperSatSelf4_i,\
                      supersatNH3_i,supersatH2S_i);
	  
      if (USE_SOL_CLOUD)
            init_soln_cloud(1);
      
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
             return ;
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
      
      ofp1=fopen("tcm.out","w");
      funoutput = (double *) malloc(top*sizeof(double));
      /***********************************************
       Output data to Matlab (tm)
      **********************************************/
      plhs[0] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[1] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[2] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[3] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[4] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[5] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[6] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[7] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[8] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[9] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[10] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[11] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[12] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[13] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[14] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[15] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[16] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[17] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[18] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[19] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[20] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[21] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[22] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      plhs[23] = mxCreateDoubleMatrix(top+1, 1, mxREAL);
      /*1 2  3   4   5    6    7    8     9      10    11   12   13   14   15   16   17  18 19 20       21        */
     /*P T XH2 XHe XH2S XNH3 XH2O XCH4 XPH3  clouds DNH4SH DH2S DNH3 DH2O DCH4 DPH3 DSOL g mu refr_w/o refr_ w/  */
     P_out = mxGetPr(plhs[0]); //1
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
     for (j =0 ;j<=top;++j)
      {     
         if(layer[j].P<P0 && !cross_my_P0)
            {
                  z_offset=layer[j].z;
                  cross_my_P0=1;
            }
      }

/*******************Print sutff to output file, just for fun, and fill the matlab arrays. ***********/
      for(j = 0;j<=top;++j)
      {
            
            T = layer[j].T;
            P = layer[j].P;
            Pr= layer[j].P_real;
            P_real[j]=Pr;
            P_out[j]=P;
            T_out[j]=T;
            XH2_out[j]=layer[j].XH2;
            XHe_out[j]=layer[j].XHe;
            XH2S_out[j]=layer[j].XH2S;
            XNH3_out[j]=layer[j].XNH3;
            XH2O_out[j]=layer[j].XH2O;
            XCH4_out[j]=layer[j].XCH4;
            XPH3_out[j]=layer[j].XPH3;
            clouds_out[j]=layer[j].clouds;
            DNH4SH_out[j]=layer[j].DNH4SH;
            DH2S_out[j]=layer[j].DH2S;
            DNH3_out[j]=layer[j].DNH3;
            DH2O_out[j]=layer[j].DH2O;
            DCH4_out[j]=layer[j].DCH4;
            DPH3_out[j]=layer[j].DPH3;
            DSOL_out[j]=layer[j].DSOL;
            g_out[j]=layer[j].g;
            mu_out[j]=layer[j].mu;
            DSOL_NH3_out[j]=layer[j].DSOL_NH3;
            refr = layer[j].XH2*REFR_H2  +layer[j].XHe*REFR_He  +layer[j].XH2S*REFR_H2S+layer[j].XPH3*REFR_PH3;
            refr+= layer[j].XNH3*REFR_NH3+layer[j].XCH4*REFR_CH4+layer[j].XH2O*REFR_H2O(T);
            refr*= ( P*(293.0/T) );
			
            refr2= (layer[j].XH2*REFR_H2  +layer[j].XHe*REFR_He)*P*(293.0/T);
			refr_w_o_out[j]=refr;
			refr_w_out[j]=refr2;
			z_out[j]=layer[j].z - z_offset;
             /*1 2      3   4   5    6    7    8     9      10    11   12   13   14   15   16   17  18 19 20       21        */
             /*P T dz   XH2 XHe XH2S XNH3 XH2O XCH4 XPH3  clouds DNH4SH DH2S DNH3 DH2O DCH4 DPH3 DSOL g mu refr_w/o refr_ w/  */
             fprintf(ofp1,"%7.3f\t%6.3f\t%7.3f\t%10.3e\t%10.3e\t",P,T,z_out[j],layer[j].XH2,layer[j].XHe);
            
             fprintf(ofp1,"%10.3e\t%10.3e\t%10.3e\t%10.3e\t%10.3e\t",layer[j].XH2S,layer[j].XNH3,layer[j].XH2O,layer[j].XCH4,layer[j].XPH3);
            
             fprintf(ofp1,"%06ld\t%10.3e\t%10.3e\t%10.3e\t%10.3e\t%10.3e\t%10.3e\t%10.3e\t",layer[j].clouds,layer[j].DNH4SH,layer[j].DH2S,\
                    layer[j].DNH3,layer[j].DH2O,layer[j].DCH4,layer[j].DPH3,layer[j].DSOL);
            
             fprintf(ofp1,"%g\t%g\t%g\t%g\n",layer[j].g,layer[j].mu,refr,refr2);
      }
      
      fclose(ofp1);//close file pointer
/***********************************************************************/
      printf("End of Thermo-chemical model, return to Matlab.\n\n");
	  
}

/*************************** End of Main Matlab Function*****************/




double getMatlabScalar (const mxArray* ptr) {

  // Make sure the input argument is a scalar in double-precision.
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != 1)
    mexErrMsgTxt("The input argument must be a double-precision scalar");

  return *mxGetPr(ptr);
}

int getMatlabInt (const mxArray* ptr) {

  return *mxGetPr(ptr);
}

double& createMatlabScalar (mxArray*& ptr) { 
  ptr = mxCreateDoubleMatrix(1,1,mxREAL);
  return *mxGetPr(ptr);
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
int getline(FILE *fp, char *s, int lim)
{
      char c;
      int i;
	  
      for (i=0;i<lim-2 && (c=getc(fp))!=EOF && c!='\n';++i)
           s[i] = c;
      s[i] = '\0';

      return (i);
}
/**************That's All Folks!********************************/
