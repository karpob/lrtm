#include "layers.h"
extern float *TfL, *PfL, lawf[];
/********************************************************************************************/
/********************************************************************************************/
// float get_dP_using_dP()
//                        This function calculates the pressure step and altitude step 
//                        when the user specifies that he/she wants to use a step in pressure
//
//
//                Input: 
//                    --> int j: the index associated with layer (see model.h)
//                    --> int *eflag: error flag return to main
//                    --> dP_init: intial pressure step (bars)
//                    --> dP_fine: fine pressure step (bars)
//                    --> P_fine_start: bottom boundary for fine pressure step (above this fine)
//                    --> P_fine_stop: upper boundary for fine pressure step
//
//                Output: 
//                    <-- dP: the pressure step (bars)
//                    <-- layer data structure (see model.h)
/********************************************************************************************/

float get_dP_using_dP(int j, int *eflag, float dP_init, float dP_fine, \
                      float P_fine_start, float P_fine_stop)
{
      float dz, P, T, H, dP, new_P_fine, new_P_coarse;
      H = R*layer[j-1].T/(layer[j-1].mu*layer[j-1].g);
      new_P_fine=layer[j-1].P - dP_fine;
      new_P_coarse=layer[j-1].P - dP_init;
      
      
      if((new_P_fine > P_fine_start)||(new_P_coarse > P_fine_start))
      {
       dP= -1.0*dP_init;
       P= layer[j-1].P + dP;
      }
      else if((new_P_fine <= P_fine_start)&&(new_P_fine > P_targ))
      {
          dP= -1.0*dP_fine;
          P= layer[j-1].P + dP;
      }
      else if (P<= P_targ && *eflag != 97)
      {
           *eflag = 98;
           P = P_targ;
           dP = PfL[1]-P_targ;
           jj=j;
      }
      else if (P <=P_targ)
      {
           dP=PfL[j-jj+1]-layer[j-1].P;
           P=layer[j-1].P+dP;
          
      }
      if(P <= P_term)
      {
           *eflag = 99;
           P= P_term;
           dP=P - layer[j-1].P;
      }
      dz= -1.0*dP/((layer[j-1].P/H)*(1e5));
      layer[j].P = P;
      layer[j].z = layer[j-1].z + dz;
      if(P<P0 && !CrossP0)
      {
           zP0= layer[j].z;
           CrossP0=1;
      }
      return (dP);
}
/********************************************************************************************/
/********************************************************************************************/
/********************************************************************************************/
