#include "layers.h"
float specific_heat(int j, float T, float P);
/********************************************************************************************/
/********************************************************************************************/
// float get_dT()  This function uses Pressure temperature and the latent heats associated with
//                 each layer and calculates the temperature step associated with a step in
//                 pressure.
//
//          Input:
//            --> int j: index associated with layer data structure (see model.h)
//            --> float T: Temperature of layer in K
//            --> float P: Pressure of the layer in bars
//            --> float dP: Pressure step in bars
//            --> *LX: Array of latent heats times mole fraction (see DeBoer eqn 3.19)
//            --> *L2X: Array of latent heats squared times mole fraction (see DeBoer eqn 3.19)
/********************************************************************************************/

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
            dT = dT_num/dT_den; // eqn 3.19 in DeBoer's thesis
            layer[j].T = layer[j-1].T + dT;
      }
      return (dT);
}
/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
