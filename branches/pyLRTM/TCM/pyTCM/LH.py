#include "model.h"


def coefficients(s):
        """
        /************************************************************************
        void coefficients(): This function returns the coefficients for 
                        saturation pressure and latent heat.
                        See DeBoer's thesis Table 3.1.
          Input:
             -->char s[]: String with desired components NH3_over_NH3_ice etc.
             -->float T: Temperature (not really used, there for fun?)
          Output:   
             <--float a[]: Array of coefficients returned 
        **********************************************************************/
        """
        from numpy import zeros
        a=zeros([6])
        if( s=="NH3_over_NH3_ice" ):  #/*Briggs and Sackett*/
                a[1] = -4122.0;
                a[2] = 27.8632;
                a[3] = -1.8163;
                a[4] = 0.0;
                a[5] = 0.0;
      
        elif( s=="NH3_over_liquid_NH3" ):  #/*Briggs and Sackett*/

                a[1] = -4409.3512;
                a[2] = 63.0487;
                a[3] = -8.4598;
                a[4] = 5.51E-3;
                a[5] = 6.80E-6;
      
        elif( s=="H2S_over_H2S_ice" ):#  /*Allen, Giauque/Blue*/

                a[1] = -2920.6;
                a[2] = 14.156;
                a[3] = 0.0;
                a[4] = 0.0;
                a[5] = 0.0;
      
        elif( s=="H2S_over_liquid_H2S" ):#  /*Allen, Giauque/Blue*/
                a[1] = -2434.62;
                a[2] = 11.4718;
                a[3] = 0.0;
                a[4] = 0.0;
                a[5] = 0.0;
      
        elif ( s=="H2O_over_ice" ):  #/*Briggs and Sackett*/
                a[1] = -5631.1206;
                a[2] = -22.1791;
                a[3] = 8.2312;
                a[4] = -3.861449e-2;
                a[5] = 2.77494e-5;
      
        elif ( s=="H2O_over_water" ):#  /*Briggs and Sackett*/
                a[1] = -2313.0338;
                a[2] = -177.848;
                a[3] = 38.053682;
                a[4] = -0.13844344;
                a[5] = 7.4465367e-5;
      
        elif( s=="CH4_over_CH4_ice" ):#  /*dePater and Massie*/
                a[1] = -1168.1;
                a[2] = 10.710;
                a[3] = 0.0;
                a[4] = 0.0;
                a[5] = 0.0;
 
        elif ( s=="CH4_over_liquid_CH4" ):#  /*dePater and Massie*/
                a[1] = -1032.5;
                a[2] = 9.216;
                a[3] = 0.0;
                a[4] = 0.0;
                a[5] = 0.0;
      
        elif( s=="NH4SH" ):  #/*Lewis*/
                a[1] = -10834.0;
                a[2] = 34.151;
                a[3] = 0.0;
                a[4] = 0.0;
                a[5] = 0.0;
      
        elif ( s=="PH3_over_PH3_ice" ):  #/*Orton/Kaminski*/
                a[1] = -1830.0;
                a[2] = 9.8225;
                a[3] = 0.0;
                a[4] = 0.0;
                a[5] = 0.0;
      
        else:
                a[1] = 0.0;
                a[2] = 0.0;
                a[3] = 0.0;
                a[4] = 0.0;
                a[5] = 0.0;
                print "Invalid phase\n"
      
        return a



def sat_pressure(s, T):
        from math import exp,log
        #a=zeros([6])
        a=coefficients(s)
        SP = a[1]/T + a[2] + a[3]*log(T) + a[4]*T + a[5]*T*T;
        if(s=="H2S_over_H2S_ice" ): return (exp(float(SP)) );
	else: return exp(float(SP)) ;




def latent_heat(s, T):
        """
        /*************************************************************************
        float latent_heat(): calculates the Latent Heat of a given
                       condensible species.

         Input:
                -->char s[]: String of the desired condensation mechanism 
                        H2S_over_H2S_ice, etc.
                -->float T: Temperature (K)
         Output:
               <--L: Latent Heat.
        *************************************************************************/
        """
        from modelParams import R
        a=coefficients(s)
        L = R*(-a[1] + a[3]*T + a[4]*T*T + 2.0*a[5]*T*T*T);
        return  L 


