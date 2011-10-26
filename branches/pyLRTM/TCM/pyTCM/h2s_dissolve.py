

def h2s_dissolve(j,layer):
        """
        /*   Returns CH2S.  Adapted from Romani's FORTRAN code.  */

        /***********************************************************************
        C
        C   Paul Romani  June 1980
        C
        C   This subroutine calculates how much H2S disolves into the already
        C   formed aqueous ammonia cloud. It does so by calculating, via an
        C   iterative method, an emperical expression that gives the partial
        C   pressure of H2S above an H2O-NH3-H2S solution as a function of
        C   temperature and the concentrations of NH3 and H2S in moles/liter.
        C   The source of the empircal expression is Leyko, Bull. Acad. Polon.
        C   Sci. Ser. Chim., 12, 275-276, 1964. The density of the aqueous
        C   ammonia cloud is from an equation given by Croft, Lunine, and Kargel
        C   in Icarus 1987 (1988?). The density is a function of the mass
        C   fraction of ammonia only, any temperature and pressure efects are
        C   ignored. If convergence of the concentration of H2S in the solution
        C   is not achieved (to 0.01% within 50 iterations), or the amount of
        C   H2S in the atmosphere becomes zero or less, no H2S is allowed to
        C   disolve into the solution.
        C
        ************************************************************************/
        /***********************************************************************
        h2s_dissolve(): Calculate the mole fraction of H2S that dissolves into
                   the solution cloud
                
         Input:
             -->int j: layer index
             -->SPH2S: Saturation vapor pressure of H2S

         Output:
             <--CH2S: the mole fraction of H2S that dissolves in the cloud.
        ***********************************************************************/
        """
        from math import log,exp
        from modelParams import R,AMU_H2,AMU_He,AMU_H2S,AMU_NH3,AMU_H2O ,AMU_CH4,AMU_PH3,AMU_NH4SH,TRIPLEPT_H2S,TRIPLEPT_NH3,TRIPLEPT_CH4,TRIPLEPT_PH3,TRIPLEPT_H2O,GONE,ZERO       
        bars_mmHg=750.062
        cm3lit = 1000.0;
        power = 1.0/( 1.13 + 1.8953 );
        oldC = 0.0;
        XH2S = layer['XH2S'][j];
        T = layer['T'][j];
        P = layer['P'][j];
        dXH2O = layer['XH2O'][j-1] - layer['XH2O'][j];
        dXNH3 = layer['XNH3'][j-1] - layer['XNH3'][j];
	SPH2S = XH2S*P;
        
        #/* Find the density of the solution cloud */
        X = ( dXNH3*AMU_NH3 ) / ( dXNH3*AMU_NH3 + dXH2O*AMU_H2O );
        DENSOL = 0.9991 + X*( -0.4336 + X*( 0.3303 + X*( 0.2833 + X*( -1.9716 + X*( 2.1396 - X*0.7294) ) ) ) );
        #/*  Calculte the partial pressure of H2S in the atmosphere in mm Hg.  */
        PH2S = ( bars_mmHg*P*XH2S );

        #/*  Find the concentration of the solution cloud in moles of NH3 / liter */
        CM3SOL = -( dXH2O*AMU_H2O + dXNH3*AMU_NH3 ) / DENSOL;
        VOLSOL = CM3SOL / cm3lit;
        CONSOL = -dXNH3 / VOLSOL;

        #/*  Iterate through the equation until a constant
        #    CH2S is found. If not do nothing.  */
        FUNTMP = exp( 22.221 - 5388.0/T );
        FUNNH3 = pow(CONSOL,1.8953);
        F = FUNNH3 / FUNTMP;
        for i in range(0,50):
                CH2S = pow( PH2S*F, power );
                try:eCH2S = abs( (CH2S - oldC) / CH2S ) * 100.0;
                except:eCH2S=0.0
                if ( eCH2S <= 0.001 ): 
                         
                        SPH2S-=( CH2S*VOLSOL )*P
                        return CH2S,SPH2S
                
                if ( PH2S <= 0.0 ):
                        return 0.0,SPH2S;
                oldC = CH2S;
                PH2S = ( XH2S - (CH2S*VOLSOL) ) * ( P/bars_mmHg );        
        return 0.0,SPH2S
