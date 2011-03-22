#include "layers.h"
extern float *TfL, *PfL, lawf[];
/********************************************************************************************/
/********************************************************************************************/
// float get_dP_using_dz()
//                         This function calculates the pressure step for a given value of dz
//             Input:
//                 --> int j: value of the layer index of layer data structure (see model.h)
//                 --> int eflag: error flag which is returned to main 
//                 --> dz : altitude step in km
//             Output:
//                 <-- dP: the pressure step in bars
//                 <-- layer data structure see model.h
/********************************************************************************************/

float get_dP_using_dz(int j,int *eflag, float dz)
{
      float P, H, dP;

      if (AutoStep)
      { 
            /*dz = log10(layer[j-1].P + 1.6); */
           /* dz = log10(layer[j-1].P + 5.0); */
		   printf("%f",AutoStep_constant);
		   dz=log10(layer[j-1].P + AutoStep_constant );  //what dave actually uses for Priscilla's stuff
            if (dz > 2.0) dz=2.0;
      }

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
            *eflag= 98;
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
/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
