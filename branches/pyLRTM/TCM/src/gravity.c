#include "layers.h"
/******************************************************************************************/
// float gravity()
//          Input: 
//            --> int j: index associated with layer data structure.
//          Output:
//            <-- returns float with the gravity associated with pressure level j
//                in (cm/sec^2)
//
/*****************************************************************************************/
float gravity(int j)
{
      float gamma, b;

      gamma = 2.0*log(layer[j].P/P0)*(R*layer[j].T)/(R0*layer[j].mu);
      b = g0 + gamma;
      return  0.5*( b + sqrt( SQ(b) - SQ(gamma) ) );
}
/********************************************************************************************/
