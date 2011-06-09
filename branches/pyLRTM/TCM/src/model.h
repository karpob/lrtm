#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#define nbgetc()        bioskey(0) & 0xFF
#define MAXLAYERS       5000
#define R               8.3143E7    /* Universal gas constant [erg/K/mol] */
#define AMU             1.66056E-24 /* Atomic Mass Unit in grams */
#define AMU_H2          2.01594     /*updated karpowicz*/
#define AMU_He          4.0026      /*updated bk*/  /* constituent masses [amu] */
#define AMU_H2S         34.076
#define AMU_NH3         17.030
#define AMU_H2O         18.015486068204702 /*bk*/ 
#define AMU_CH4         16.0428  /*bk*/
#define AMU_PH3         33.997
#define AMU_NH4SH       51.110
#define H2S_SOLAR       3.18E-5
#define NH3_SOLAR       1.47E-4     /* solar abundances */
#define H2O_SOLAR       1.17e-3     /* dePater, Romani, Atreya */
#define CH4_SOLAR       7.06e-4
#define TRIPLEPT_H2S    187.61
#define TRIPLEPT_NH3    195.5       /* triple points */
#define TRIPLEPT_CH4    90.7
#define TRIPLEPT_PH3    0.0
#define TRIPLEPT_H2O    273.16
#define REFR_H2         124.43                 /* refractivity coefficients */
#define REFR_He         35.832                 /* R = REFRACT_X*PX*(293/T)  */
#define REFR_H2S        2247.0                 /*      n = (R/1E6) + 1      */
#define REFR_CH4        413.0
#define REFR_PH3          0.0
#define REFR_NH3        2700.0  /*Spilker--VERY average since near resonances*/
#define REFR_H2O(T)     ( 245.0 + 1.28E6/(T) ) /*Janssen p218--or is the h2o cloud really almost an ocean,see liquid water p298 UFM*/
#define GONE            1e-30       /* value when a constituent is considered gone */
#define ZERO            1e-30       /* zero cloud density (for plotting) */
#define COUNT_CLOUD     1e-3       /* count cloud opacity if greater than this number */
#define MAXTRIES        100          /* Maximum number of iterations to fit deep temperature */
#define TLIMIT          0.05         /* allowed % error in T */
#define TAULIMIT        20          /* max tau used for T */
#define SQ(x)           ( (x)*(x) )
#define CU(x)           ( (x)*(x)*(x) )



struct ATM_LAYER {
      float alpha;            /* layer absorptivity  */
      float tau;              /* layer optical depth */
      float DW;               /* disk averaged weighting function */
      float T;                /* layer temperature   */
      float P;                /* layer pressure  (ideal pressure)    */
      float P_real;           /* Real pressure */
      float z;                /* altitude */
      double XH2;
      double XHe;
      double XH2S;            /* mixing ratios */
      double XNH3;
      double XH2O;
      double XCH4;
      double XPH3;
      long clouds;            /* cloud structure: 111111: PH3,CH4,H2S,NH3,NH4SH,H2O: 1=liq,2=ice */
      float DSOL;
      float DNH4SH;
      float DH2S;
      float DNH3;             /* cloud densities     */
      float DH2O;
      float DCH4;
      float DPH3;
      float g;                /* gravity             */
      float mu;               /* average mass        */
      double DSOL_NH3;        /*in g/cm^3 density of cloud that is composed of NH3 */
      int first;
      float q_c;
      float q_c_nh3;
      float q_c_nh3_ice;
      int first_nh3;
      

                  };
