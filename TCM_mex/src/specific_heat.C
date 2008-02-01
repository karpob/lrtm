/****************************sp_heat.c**************************************/
#include "model.h"
//#define  NUMBER_H2  16
//#define  fp         0.0  /* Para fraction.  0.0=equilibrium,  -1.0=intermediate
//                                            0.25 = normal     */

/*extern struct ATM_LAYER layer[MAXLAYERS];*/
extern struct ATM_LAYER *layer;
extern float Hydrogen_Curve_Fit_Select;

float specific_heat(int j, float T, float P)
{
      int i;
      float XH2, XHe, XH2S, XNH3, XH2O, XCH4, XPH3;
      float Cp_H2, Cp_He, Cp_H2S, Cp_NH3, Cp_H2O, Cp_CH4, Cp_PH3, Cp;

/**************versions of CpH2*******************************************/
/*      fit for equilibrium   ***/
      	if (T<120.0 && Hydrogen_Curve_Fit_Select==0.0)
        	Cp_H2 = 8.8477 - 1.0421*T + 6.1631e-2*SQ(T) - 1.6961e-3*pow(T,3.0) +
                    2.6360e-5*pow(T,4.0) - 2.5341e-7*pow(T,5.0) + 1.5665e-9*pow(T,6.0) -
                    6.2470e-12*pow(T,7.0) + 1.5544e-14*pow(T,8.0) - 2.19549e-17*pow(T,9.0) +
                    1.344336e-20*pow(T,10.0);
      	else if (T>=120.0 && T<=130.0 && Hydrogen_Curve_Fit_Select==0.0)
            Cp_H2 = -6.875E-3*(T-120.0) + 3.31494;
      	else if (T>130.0 && T<=358.0 && Hydrogen_Curve_Fit_Select==0.0)
            Cp_H2 = 3.0795 + 0.001174*T;
      	else if (T > 358.0 && Hydrogen_Curve_Fit_Select==0.0)
            Cp_H2 = 3.5;
       

/***      fit for intermediate or frozen equilibrium   ***/
  
      	if (T<=314.0 && Hydrogen_Curve_Fit_Select==-1.0)
            Cp_H2 = 2.6114 - 9.3144e-3*T + 2.0612e-4*SQ(T) - 1.217e-6*pow(T,3.0) +
                    3.135e-9*pow(T,4.0) - 3.0537e-12*pow(T,5.0);
      	else if (Hydrogen_Curve_Fit_Select==-1.0)
            Cp_H2 = 3.5;

/***      fit for normal hydrogen      ***/
	if (T<=310.0 && Hydrogen_Curve_Fit_Select==0.25)
           	Cp_H2 = 2.5785 - 5.6608e-3*T + 9.4152e-5*SQ(T) - 1.8306e-7*pow(T,3.0) -
                   	6.239e-10*pow(T,4.0) + 1.696e-12*pow(T,5.0);
        else if (Hydrogen_Curve_Fit_Select==0.25)
            	Cp_H2 = 3.5;
/*************************************************************************/
      Cp_He = 2.503;   /*                  */
      Cp_H2S = 4.013;
      Cp_NH3 = 4.459;  /*Briggs and Sackett*/
      Cp_H2O = 4.0; 
      Cp_CH4 = 4.5;    /*                  */
	Cp_PH3 = 0.0;

      XH2  = layer[j-1].XH2;
      XHe  = layer[j-1].XHe;
      XH2S = layer[j-1].XH2S;
      XNH3 = layer[j-1].XNH3;
      XH2O = layer[j-1].XH2O;
      XCH4 = layer[j-1].XCH4;
	XPH3 = layer[j-1].XPH3;
      Cp = (XH2*Cp_H2 + XHe*Cp_He + XH2S*Cp_H2S + XNH3*Cp_NH3 + XPH3*Cp_PH3)*R;

      return ( Cp );
}
