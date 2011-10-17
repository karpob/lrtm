#include "layers.h"
#******************************************************************************************/
# float gravity()
#          Input: 
#            --> int j: index associated with layer data structure.
#          Output:
#            <-- returns float with the gravity associated with pressure level j
#                in (cm/sec^2)
#
#*****************************************************************************************/
def gravity(layer,j):
      from numpy import log,sqrt
      from modelParams import R
      P0=layer['P0']
      R0=layer['R0']
      g0=layer['g0']
      gamma = 2.0*log(layer['P'][j]/P0)*(R*layer['T'][j])/(R0*layer['mu'][j])
      b = g0 + gamma
      return  0.5*( b + sqrt( b**2 - gamma**2 ) )


