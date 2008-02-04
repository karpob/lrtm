/****************************sat_lat.c************************************/
#include "model.h"

void coefficients(char s[], float a[], float T)
{
      if( !strcmp(s,"NH3_over_NH3_ice") )  /*Briggs and Sackett*/
      {
            a[1] = -4122.0;
            a[2] = 27.8632;
            a[3] = -1.8163;
            a[4] = 0.0;
            a[5] = 0.0;
      }
      else if( !strcmp(s,"NH3_over_liquid_NH3") )  /*Briggs and Sackett*/
      {
            a[1] = -4409.3512;
            a[2] = 63.0487;
            a[3] = -8.4598;
            a[4] = 5.51E-3;
            a[5] = 6.80E-6;
      }
      else if( !strcmp(s,"H2S_over_H2S_ice") )  /*Allen, Giauque/Blue*/
      {
            a[1] = -2920.6;
            a[2] = 14.156;
            a[3] = 0.0;
            a[4] = 0.0;
            a[5] = 0.0;
      }
      else if( !strcmp(s,"H2S_over_liquid_H2S") )  /*Allen, Giauque/Blue*/
      {
            a[1] = -2434.62;
            a[2] = 11.4718;
            a[3] = 0.0;
            a[4] = 0.0;
            a[5] = 0.0;
      }
      else if ( !strcmp(s,"H2O_over_ice") )  /*Briggs and Sackett*/
      {
            a[1] = -5631.1206;
            a[2] = -22.1791;
            a[3] = 8.2312;
            a[4] = -3.861449e-2;
            a[5] = 2.77494e-5;
      }
      else if ( !strcmp(s,"H2O_over_water") )  /*Briggs and Sackett*/
      {
            a[1] = -2313.0338;
            a[2] = -177.848;
            a[3] = 38.053682;
            a[4] = -0.13844344;
            a[5] = 7.4465367e-5;
      }
      else if( !strcmp(s,"CH4_over_CH4_ice") )  /*dePater and Massie*/
      {
            a[1] = -1168.1;
            a[2] = 10.710;
            a[3] = 0.0;
            a[4] = 0.0;
            a[5] = 0.0;
      }
      else if ( !strcmp(s,"CH4_over_liquid_CH4") )  /*dePater and Massie*/
      {
            a[1] = -1032.5;
            a[2] = 9.216;
            a[3] = 0.0;
            a[4] = 0.0;
            a[5] = 0.0;
      }
      else if( !strcmp(s,"NH4SH") )  /*Lewis*/
      {
            a[1] = -10834.0;
            a[2] = 34.151;
            a[3] = 0.0;
            a[4] = 0.0;
            a[5] = 0.0;
      }
      else if ( !strcmp(s,"PH3_over_PH3_ice") )  /*Orton/Kaminski*/
      {
            a[1] = -1830.0;
            a[2] = 9.8225;
            a[3] = 0.0;
            a[4] = 0.0;
            a[5] = 0.0;
      }
      else
      {
            a[1] = 0.0;
            a[2] = 0.0;
            a[3] = 0.0;
            a[4] = 0.0;
            a[5] = 0.0;
            printf("Invalid phase\n");
      }
      return;
}

#define SUPERSAT 1.00
float sat_pressure(char s[], float T)
{
      float SP, a[6];

      coefficients(s,a,T);
      SP = a[1]/T + a[2] + a[3]*log(T) + a[4]*T + a[5]*T*T;
      if( !strcmp(s,"H2S_over_H2S_ice") ) return ( SUPERSAT*exp(SP) );
	else return ( exp(SP) );
}
#undef SUPERSAT

float latent_heat(char s[], float T)
{
      float L, a[6];

      coefficients(s,a,T);
      L = R*(-a[1] + a[3]*T + a[4]*T*T + 2.0*a[5]*T*T*T);
      return ( L );
}

