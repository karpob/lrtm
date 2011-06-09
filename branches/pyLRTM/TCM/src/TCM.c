/**************************** TCM.C **************************************/
/*  This model requires the data files:                                  */
/*                                       inputsoln.dat 
/*                                       and                                */
/*************************************************************************/
#include "model.h"

/*Globals*/
struct ATM_LAYER *layer;
int hereonout;
char use_lindal;
int NH4SH_cloud_base=0;
float P_targ, T_targ, P_term, g0, R0, P0, supersatNH3, supersatH2S,use_lindal_in;
float *TfL, *PfL, lawf[3],zP0;
double XCO;
int AutoStep,CrossP0;
double AutoStep_constant;
int use_dz,jj;
/*set 0 to turn off, 1 to turn on*/
#define USE_SOL_CLOUD  1      /*consider its formation?*/
double P_temp,T_temp;
float z_offset;
//float Hydrogen_Curve_Fit_Select;
int sol_cloud=USE_SOL_CLOUD;

/*       Prototypes      */
float get_dP_using_dz(int j, int *eflag, float dz);
float get_dP_using_dP(int j, int *eflag, float dP_init, float dP_fine, float P_fine_start, float P_fine_stop);
void new_layer(int j, float dz, int *eflag,float dP_init, float dP_fine, float P_fine_start, \
               float P_fine_stop,float frain, float select_ackerman,float Hydrogen_Curve_Fit_Select);
int new_layer_original(int j, float dz, int eflag, float dP_init, float dP_fine, float P_fine_start, float P_fine_stop);
float specific_heat(int j, float T, float P);
int init_atm(int n,double XHe,double XH2S,double XNH3,double XH2O,double XCH4,double XPH3,double P_temp,double T_temp,float g0_i,float R0_i, float P0_i,char use_lindal_i, float T_targ_i, float P_targ_i, float P_term_i,int n_lindal_pts_i,float SuperSatSelf1_i,float SuperSatSelf2_i, float SuperSatSelf3_i, float SuperSatSelf4_i,float supersatNH3_i,float supersatH2S_i);
int intoTheVoid( float dz, double XHe_i,double XH2S_i,double XNH3_i,double XH2O_i,double XCH4_i,double XPH3_i,double XCO, float P_temp,float T_temp,float g0_i,float R0_i,float P0_i,float T_targ_i,float P_targ_i,float P_term_i,float use_lindal_in,int n_lindal_pts_i,float SuperSatSelf1_i,float SuperSatSelf2_i,float SuperSatSelf3_i,float SuperSatSelf4_i,float supersatNH3_i,float supersatH2S_i,int AutoStep_constant,float Hydrogen_Curve_Fit_Select,float dP_init,float dP_fine,float P_fine_start,float P_fine_stop,int use_dz,float frain,float select_ackerman);
double getFloatValues(int i,int j);
long getCloudFlags(int i);
void init_soln_cloud(int mode);

//void dt_tag(FILE *fp);
/*
int main(void)
{
float dz=1;
double XHe_i=0.13464647667;
double XH2S_i=6.60368070295e-05;
double XNH3_i=0.000390217496083;
double XH2O_i=0.00547488011318;
double XCH4_i=0.00180100382808;
double XPH3_i=5.14572522308e-07;
double XCO=0.0;
float P_temp=1000.0;
float T_temp=1583.7824707;
float g0_i=2330.0;
float R0_i=7149200000.0;
float P0_i=1.0;
float T_targ_i=164.944290161;
float P_targ_i=1.0;
float P_term_i=0.141;
float use_lindal_i=1;
int n_lindal_pts_i=61;
float SuperSatSelf1_i=0.0;
float SuperSatSelf2_i=0.0;
float SuperSatSelf3_i=0.0;
float SuperSatSelf4_i=0.0;
float supersatNH3_i=0.0;
float supersatH2S_i=0.0;
int AutoStep_constant=8;
float Hydrogen_Curve_Fit_Select=-1.0;
float dP_init=10.0;
float dP_fine=0.25;
float P_fine_start=13.0;
float P_fine_stop=1.0;
int use_dz=1;
float frain=3.0;
float select_ackerman=0;
int i;
int nlay=0;
float Val;

nlay=intoTheVoid( dz, 
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
                   select_ackerman);
	for(i=1;i<nlay;++i)
	{
		Val=getFloatValues(16,i);
		printf("Val %f\n",Val);
		printf("all done.");
	}
} */
/***********************************************************************************
Main function which runs the bulk of TCM, returns the number of layers.
Use helper functions getCloudFlags and getFloatValues to grab values from the TCM.
**********************************************************************************/
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
      int  top, j, jcntr=0, found,cross_my_P0=0;
      double P, T, T_err=100.0;
      char use_lindal_i='Y';
      
      int eflag=0;
      int tmp;	  
      
      
      
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
             
              new_layer(j,dz,&eflag, dP_init, dP_fine, P_fine_start, P_fine_stop, frain,select_ackerman,Hydrogen_Curve_Fit_Select);
	     
	     
              P = layer[j].P;
            
	     
	      
	      
              //printf("error flag Pressure, level,T_targ %d %f %d %f\n",layer[j].eflag,P,j,T_targ,layer[j].T);
                  
                  if (eflag == 98 )  /* check target temperature */
                  {                                         /* eflag       P      */
                        eflag = 97;                         /*  97      < P_targ  */
                        T = layer[j].T;                     /*  99      <=P_term  */
                        T_err = 100.*(T - T_targ)/T_targ;  /*  98      <=P_targ  */
                        printf("Terr %f %f %d \n",T_err,TLIMIT,j);
			//printf("eflag %d",eflag);
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
      
      /*if (P > P_term)
      {
             printf("Abnormal termination. P= %f bigger than Pterm=%f                       \n",P,P_term);
             printf("jj val %d",jj);
             return 0;
      } */

      top = j-1;
      
      if (USE_SOL_CLOUD)
            init_soln_cloud(0); //free solution cloud memory allocation
  
      //if (use_lindal == 'Y' || use_lindal == 'y')
     // {
     //       free(PfL); //free lindal profile memory allocation
     //       free(TfL);
     // }
      
     for (j =0 ;j<=top;++j)
      {     
         if(layer[j].P<P0 && !cross_my_P0)
            {
                  z_offset=layer[j].z;
                  cross_my_P0=1;
            }
      }
      

      return (top);		  
}


/*********************************************************
	Returns Cloud Flags Bits packed with information
	about cloud type.
**********************************************************/

long getCloudFlags(int i)
{      
	return layer[i].clouds;
}
/*********************************************************
Helper function to grab values out of the TCM one by one.
It seems SWIG will only allow passing via this method?
Worth revisiting if we can pass a pointer or something...
**********************************************************/
double getFloatValues(int i,int j)
{

		if(i==0) return layer[j].P;
		if(i==1) return layer[j].T;
		if(i==2)return layer[j].z-z_offset;
		if(i==3) return layer[j].XH2;
		if(i==4) return layer[j].XHe;
		if(i==5) return layer[j].XH2S;
		if(i==6) return layer[j].XNH3;
		if(i==7)return layer[j].XH2O;
		if(i==8)return layer[j].XCH4;
		if(i==9)return layer[j].XPH3;
		if(i==10)return layer[j].DNH4SH;
		if(i==11)return layer[j].DH2S;
		if(i==12)return layer[j].DNH3;
		if(i==13)return layer[j].DH2O;
		if(i==14)return layer[j].DCH4;
		if(i==15)return layer[j].DPH3;
		if(i==16)return layer[j].DSOL;
		if(i==17)return layer[j].g;
		if(i==18)return layer[j].mu;
		if(i==19)return layer[j].DSOL_NH3;
		if(i==20)return layer[j].P_real;
                
		printf("exceed i dimension");
		return 0;
}


/**************That's All Folks!********************************/
