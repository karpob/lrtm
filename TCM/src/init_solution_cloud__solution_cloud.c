/* Initialize solution cloud (read in Romani's data file), and calculate concentration of NH3 that condenses in solution cloud*/
#include "model.h"
void init_soln_cloud(int mode);
float solution_cloud(float TC, float PNH3, float PH2O, float *SPNH3, float *SPH2O);

/*  These declarations are for computing the H2O-NH3 solution cloud  */
#define ROMANI        1           /* 1=use Romani, 0=use Briggs/Sackett */
static int KA=0, KLT=0, KHT=0;
static float S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12, S13, S14, S15;
static float *YA,   *CA,   *RTA,   *GTA;
static float *YWLT, *CWLT, *RTWLT, *GTWLT;
static float *YWHT, *CWHT, *RTWHT, *GTWHT;
static float *AA1, *AA2, *AA3, *DXA, *DYA1, *DYA2, *DYA3;
static float *ALT1,*ALT2,*ALT3,*DXLT,*DYLT1,*DYLT2,*DYLT3;
static float *AHT1,*AHT2,*AHT3,*DXHT,*DYHT1,*DYHT2,*DYHT3;
void STP(int N, float *X, float *Y1, float *Y2, float *Y3, float *A1, float *A2, float *A3, float *DX, float *DY1, float *DY2, float *DY3);
void SPN(float XX, float *ST1, float *ST2, float *ST3, int *K, int N, float *X, float *Y1, float *Y2, float *Y3, float *A1, float *A2, float *A3, float *DX, float *DY1, float *DY2, float *DY3);
void PPRES(float C1, float T1, float *VPH, float *VPN);
#define  CONVT  1.0E-06          /*Convert dynes/cm^2 to bars*/
#define  NA     20
#define  NWLT   18
#define  NWHT   20

/****************************************************************************************************
      init_soln_cloud()
                       This reads in coefficients for some spline fits
                       given by Paul Romani.  This is DeBoer's adaptation in C
                       from the original FORTRAN code. Uses the coefficients
                       in spline fits for calculation of parameters in solution_cloud().
            input: 
                   --> int mode: If 1 is selected, this will read in and plug in to 
                                 spline fits.
                                 If anything else, it will remove the coefficients from memory.
            output:
                      None, really. However, several cryptic global variables are calculated which
                       go into solution_cloud()

****************************************************************************************************/

void init_soln_cloud(int mode)
{
      int j;
	  double CCA, CCW;
      FILE *id;

      /*  File INPUSTSOLN.DAT - Input file for my equations and method  */
      if(mode)
      {
            printf("Computing solution cloud coefficients.\n");
            /* Allocate memory and initialize arrays */
            YA = (float *) calloc(NA, sizeof(float));
			CA = (float *) calloc(NA, sizeof(float));
            RTA = (float *) calloc(NA, sizeof(float));
            GTA = (float *) calloc(NA, sizeof(float));
            YWLT = (float *) calloc(NWLT, sizeof(float));
            CWLT = (float *) calloc(NWLT, sizeof(float));
            RTWLT = (float *) calloc(NWLT, sizeof(float));
            GTWLT = (float *) calloc(NWLT, sizeof(float));
            YWHT = (float *) calloc(NWHT, sizeof(float));
            CWHT = (float *) calloc(NWHT, sizeof(float));
            RTWHT = (float *) calloc(NWHT, sizeof(float));
            GTWHT = (float *) calloc(NWHT, sizeof(float));
            AA1 = (float *) calloc(NA, sizeof(float));
            AA2 = (float *) calloc(NA, sizeof(float));
            AA3 = (float *) calloc(NA, sizeof(float));
            ALT1 = (float *) calloc(NWLT, sizeof(float));
            ALT2 = (float *) calloc(NWLT, sizeof(float));
            ALT3 = (float *) calloc(NWLT, sizeof(float));
            AHT1 = (float *) calloc(NWHT, sizeof(float));
            AHT2 = (float *) calloc(NWHT, sizeof(float));
            AHT3 = (float *) calloc(NWHT, sizeof(float));
            DXA  = (float *) calloc(NA, sizeof(float));
            DXLT = (float *) calloc(NWLT, sizeof(float));
            DXHT = (float *) calloc(NWHT, sizeof(float));
            DYA1 = (float *) calloc(NA, sizeof(float));
            DYA2 = (float *) calloc(NA, sizeof(float));
            DYA3 = (float *) calloc(NA, sizeof(float));
            DYLT1 = (float *) calloc(NWLT, sizeof(float));
            DYLT2 = (float *) calloc(NWLT, sizeof(float));
            DYLT3 = (float *) calloc(NWLT, sizeof(float));
            DYHT1 = (float *) calloc(NWHT, sizeof(float));
            DYHT2 = (float *) calloc(NWHT, sizeof(float));
            DYHT3 = (float *) calloc(NWHT, sizeof(float));
					    
            if (DYHT3 == NULL)
            {
                  printf("Error:  Insufficient memory for solution cloud.\n");
                  exit(1);
            }
			
            id=fopen("inputsol.dat","r");
			
            for(j=0;j<NA;++j)
              fscanf(id,"%f %f %f %f", YA+j, CA+j, RTA+j, GTA+j);
            for(j=0;j<NWLT;++j)
              fscanf(id,"%f %f %f %f", YWLT+j, CWLT+j, RTWLT+j, GTWLT+j);
            for(j=0;j<NWHT;++j)
              fscanf(id,"%f %f %f %f", YWHT+j, CWHT+j, RTWHT+j, GTWHT+j);
            fclose(id);

/*  Get interpolation points to subroutines  */
            STP(NA,  YA,  CA,  RTA,  GTA,  AA1, AA2, AA3, DXA, DYA1, DYA2, DYA3);
            STP(NWLT,YWLT,CWLT,RTWLT,GTWLT,ALT1,ALT2,ALT3,DXLT,DYLT1,DYLT2,DYLT3);
            STP(NWHT,YWHT,CWHT,RTWHT,GTWHT,AHT1,AHT2,AHT3,DXHT,DYHT1,DYHT2,DYHT3);

/*******************************************************************/
/*  Calculate the coefficients of the vapor pressure equations to  */
/*  correct at low/high concentrations. See subrotuine PPRES for   */
/*  details.                                                       */
/*******************************************************************/
            CCA = 0.00;
            CCW = 1.00;
            SPN(CCA,&S7, &S8, &S9, &KA, NA,  YA,  CA,  RTA,  GTA,  AA1, AA2, AA3, DXA, DYA1, DYA2, DYA3);
            SPN(CCW,&S10,&S11,&S12,&KLT,NWLT,YWLT,CWLT,RTWLT,GTWLT,ALT1,ALT2,ALT3,DXLT,DYLT1,DYLT2,DYLT3);
            SPN(CCW,&S13,&S14,&S15,&KHT,NWHT,YWHT,CWHT,RTWHT,GTWHT,AHT1,AHT2,AHT3,DXHT,DYHT1,DYHT2,DYHT3);
      }
      else if(!mode)
      {
            printf("\nFreeing solution cloud coefficients.\n\n");
            /* Allocate memory and initialize arrays */
            free(YA);
            free(CA);
            free(RTA);
            free(GTA);
            free(YWLT);
            free(CWLT);
            free(RTWLT);
            free(GTWLT);
            free(YWHT);
            free(CWHT);
            free(RTWHT);
            free(GTWHT);
            free(AA1);
            free(AA2);
            free(AA3);
            free(ALT1);
            free(ALT2);
            free(ALT3);
            free(AHT1);
            free(AHT2);
            free(AHT3);
            free(DXA);
            free(DXLT);
            free(DXHT);
            free(DYA1);
            free(DYA2);
            free(DYA3);
            free(DYLT1);
            free(DYLT2);
            free(DYLT3);
            free(DYHT1);
            free(DYHT2);
            free(DYHT3);
      }

}

#define LOWER_CONC    0.0
#define UPPER_CONC    0.5
#define CONC_STEP     0.0001
/***************************************************************************************************************
   solution_cloud(): This routine calculates the mole fraction of ammonia that goes into solution for the water-solution cloud.

         Input: 
              -->TC: Temperature
              -->PNH3: Partial Pressure of NH3 in bars.
              -->PH2O: Partial Pressure of H2O in bars.
              -->SPNH3: Saturation Pressure of Ammonia
              -->SPH2O: Saturation Pressure of Water

         Output:
              <--C_sol: mole fraction of ammonia that dissolves in the solution.
**************************************************************************************************************/


float solution_cloud(float TC, float PNH3, float PH2O, float *SPNH3, float *SPH2O)
{
      int have_cloud=0;
      float C_sol, PNH3C,PH2S;

      if (ROMANI)  /* Is there a prayer of a solution cloud? */
      {
            SPN(UPPER_CONC, &S1, &S2, &S3, &KA, NA,  YA,  CA,  RTA,  GTA,  AA1, AA2, AA3, DXA, DYA1, DYA2, DYA3);
            if (TC < 300.0)
                  SPN(UPPER_CONC,&S4,&S5,&S6,&KLT,NWLT,YWLT,CWLT,RTWLT,GTWLT,ALT1,ALT2,ALT3,DXLT,DYLT1,DYLT2,DYLT3);
            else
                  SPN(UPPER_CONC,&S4,&S5,&S6,&KHT,NWHT,YWHT,CWHT,RTWHT,GTWHT,AHT1,AHT2,AHT3,DXHT,DYHT1,DYHT2,DYHT3);
            PPRES(UPPER_CONC, TC, SPH2O, SPNH3);
            *SPH2O *= CONVT;
            *SPNH3 *= CONVT;
            if (*SPH2O > PH2O)  /* No.  Then return. */
                  return (-1.0);
      }

      *SPNH3 = 0.0;
      PNH3C = PNH3;
      for (C_sol=LOWER_CONC; *SPNH3<PNH3C && C_sol<UPPER_CONC; C_sol+=CONC_STEP)
      {
            PNH3C = PNH3*(1.0-C_sol);
            if (ROMANI)
            {
                  /*  Romani version  */
                  SPN(C_sol, &S1, &S2, &S3, &KA, NA,  YA,  CA,  RTA,  GTA,  AA1, AA2, AA3, DXA, DYA1, DYA2, DYA3);
                  if (TC < 300.0)
                        SPN(C_sol,&S4,&S5,&S6,&KLT,NWLT,YWLT,CWLT,RTWLT,GTWLT,ALT1,ALT2,ALT3,DXLT,DYLT1,DYLT2,DYLT3);
                  else
                        SPN(C_sol,&S4,&S5,&S6,&KHT,NWHT,YWHT,CWHT,RTWHT,GTWHT,AHT1,AHT2,AHT3,DXHT,DYHT1,DYHT2,DYHT3);
                  PPRES(C_sol, TC, SPH2O, SPNH3);
                  *SPH2O *= CONVT;
                  *SPNH3 *= CONVT;
            }
            else if (!ROMANI)
            {
                  /*  Briggs and Sackett */
                  *SPH2O = log(1.0-C_sol) + 29.0423 + 4.0134*SQ(C_sol) - (5540.48 + 2022.11*SQ(C_sol))/TC;
                  *SPH2O = exp(*SPH2O) * CONVT;
                  *SPNH3 = log(C_sol) + 30.0048 + 4.0134*(SQ(C_sol) - 2.0*C_sol) - (4949.75 + 2022.11*(SQ(C_sol) - 2.0*C_sol))/TC;
                  *SPNH3 = exp(*SPNH3) * CONVT;
            }

            if (PH2O > *SPH2O && PNH3C > *SPNH3) /* We have a solution cloud */
                  have_cloud = 1;                /*    but what C?           */
            if (have_cloud  &&  PNH3C < *SPNH3)  /* have_cloud               */
            {
                  *SPNH3 = PNH3C;
		  //SPH2S = h2s_sol_cloud(TC, C_sol,PH2S);  /*DON'T Ignore H2S in cloud */
                  return(C_sol);
            }
      }
      return(-1.0);
}

/*********************************************************************/
/*                                                                   */
/*  NAASA  6.1.002 SPLINE   FTN   03-16-80     MARILYN WOLFSON       */
/*  INTERPOLATION BY PIECEWISE CUBIC SPLINES                         */
/*                                                                   */
/*  COMPUTES COEFFICIENTS FOR SPLINE INTERPOLATION                   */
/*  INPUT.. N = NUMBER OF DATA POINTS                                */
/*          X(I),Y1(I),Y2(I),Y3(I) = DATA TO BE INTERPOLATED         */
/*          X(1).LT.X(2).LT. ... .LT.X(N)                            */
/*  OUTPUT.. SAVED INTERNALLY FOR USE BY SPLINE                      */
/*********************************************************************/
void STP(int N, float *X, float *Y1, float *Y2, float *Y3, float *A1, float *A2, float *A3, float *DX, float *DY1, float *DY2, float *DY3)
{
      int i, N1=N-1, ib;
      float T[NA], PIV;

      /*  A(I) = S''(X(I))/6  */
      for(i=0;i<N1;++i)
      {
            DX[i] = X[i+1] - X[i];
            DY1[i] = (Y1[i+1] - Y1[i]) / DX[i];
            DY2[i] = (Y2[i+1] - Y2[i]) / DX[i];
            DY3[i] = (Y3[i+1] - Y3[i]) / DX[i];
      }
      T[0] = 0.E0;
      for(i=1;i<N1;++i)
      {
            PIV = 2.E0 * (DX[i-1] + DX[i]) - DX[i-1] * T[i-1];
            T[i] = DX[i] / PIV;
            A1[i] = (DY1[i] - DY1[i-1] - DX[i-1]*A1[i-1]) / PIV;
            A2[i] = (DY2[i] - DY2[i-1] - DX[i-1]*A2[i-1]) / PIV;
            A3[i] = (DY3[i] - DY3[i-1] - DX[i-1]*A3[i-1]) / PIV;
      }
      for(ib=1;ib<N1;++ib)
      {
            i = N1 - ib;
            A1[i] = A1[i] - T[i] * A1[i+1];
            A2[i] = A2[i] - T[i] * A2[i+1];
            A3[i] = A3[i] - T[i] * A3[i+1];
      }
}

/************************************************************/
/*  INPUT.. XX = ANY VALUE                                  */
/*  OUTPUT.. ST1,ST2,ST3 = SPLINE EVALUATED AT XX           */
/************************************************************/
void SPN(float XX, float *ST1, float *ST2, float *ST3, int *K,int N, float *X, float *Y1, float *Y2, float *Y3,float *A1, float *A2, float *A3, float *DX,float *DY1, float *DY2, float *DY3)
{
      int N1=N-1;
      float w, v;

FORTY:
      if (XX < X[*K])   goto SIXTY;
FIFTY:
      if (XX > X[*K+1]) goto SEVENTY;

      w = (XX - X[*K]) / DX[*K];
      v = 1.E0-w;
      *ST1 = w*Y1[*K+1]+v*Y1[*K]+SQ(DX[*K])*((CU(w)-w)*A1[*K+1]+(CU(v)-v)*A1[*K]);
      *ST2 = w*Y2[*K+1]+v*Y2[*K]+SQ(DX[*K])*((CU(w)-w)*A2[*K+1]+(CU(v)-v)*A2[*K]);
      *ST3 = w*Y3[*K+1]+v*Y3[*K]+SQ(DX[*K])*((CU(w)-w)*A3[*K+1]+(CU(v)-v)*A3[*K]);
      return;

SIXTY:
      *K = *K - 1;
      if (*K != -1) goto FORTY;
      *K = 0;
      *ST1 = Y1[*K] + (DY1[*K] - DX[*K]*A1[*K+1]) * (XX - X[*K]);
      *ST2 = Y2[*K] + (DY2[*K] - DX[*K]*A2[*K+1]) * (XX - X[*K]);
      *ST3 = Y3[*K] + (DY3[*K] - DX[*K]*A3[*K+1]) * (XX - X[*K]);
      return;

SEVENTY:
      *K = *K + 1;
      if (*K != N1) goto FIFTY;
      *K = N - 2;
      *ST1 = Y1[N1] + (DY1[*K] + DX[*K]*A1[*K]) * (XX - X[N1]);
      *ST2 = Y2[N1] + (DY2[*K] + DX[*K]*A2[*K]) * (XX - X[N1]);
      *ST3 = Y3[N1] + (DY3[*K] + DX[*K]*A3[*K]) * (XX - X[N1]);
      return;
}

/************************************************************************/
/*  This subroutine computes the partial equilibrium saturation vapor   */
/*  pressures of NH3 and H2O over aqueous NH3 solution at concentration */
/*  of C1 in mole fraction of NH3.                                      */
/************************************************************************/
void PPRES(float C1, float T1, float *VPH, float *VPN)
{
      float GOLT1, GOLPN, GOLPH, GOLP0N, P0N, GOLP0H, P0H, X;

      GOLT1 = log(T1);
      GOLPN = S1 + S2/T1 + S3*GOLT1;
      GOLPH = S4 + S5/T1 + S6*GOLT1;
      *VPN = exp(GOLPN);
      *VPH = exp(GOLPH);

/*  When C1 equals zero, the vapor pressure of NH3 should be zero.
    But, alas, the vapor pressure equation will return a non zero value
    due to the way the coefficients for the equation are calculated. To
    correct this, subtact off the fictous non-zero vapor pressure in a
    linear manner until at the lowest concentration for which there are
    non-interpolated/extrapolated coefficients. */
      if (C1 <= .05)
      {
            GOLP0N = S7 + S8/T1 + S9*GOLT1;
            P0N = exp(GOLP0N);
            X = (0.0500 - C1)/0.0500;
            P0N = X*P0N;
            *VPN -= P0N;
      }
/*  Same problem as above except for H2O at high concentration.*/
      if (C1 >= .95)
      {
            if (T1 <= 300.0 )
                  GOLP0H = S10 + S11/T1 + S12 * GOLT1;
            else
                  GOLP0H = S13 + S14/T1 + S15 * GOLT1;
            P0H = exp(GOLP0H);
            X = (C1 - 0.9500)/0.0500;
            P0H = X*P0H;
            *VPH -= P0H;
      }
      return;
}
