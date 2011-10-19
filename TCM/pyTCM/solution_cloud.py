

#define  NA     20
#define  NWLT   18
#define  NWHT   20




#define LOWER_CONC    0.0
#define UPPER_CONC    0.5
#define CONC_STEP     0.0001


def PPRES(C1, T1,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15): 
#/************************************************************************/
#/*  This subroutine computes the partial equilibrium saturation vapor   */
#/*  pressures of NH3 and H2O over aqueous NH3 solution at concentration */
#/*  of C1 in mole fraction of NH3.                                      */
#/************************************************************************/
        from math import log,exp
        
        try:GOLT1 = log(T1)
        except:print C1,T1
        GOLPN = S1 + S2/T1 + S3*GOLT1
        GOLPH = S4 + S5/T1 + S6*GOLT1
        
        VPN = exp(GOLPN)
        VPH = exp(GOLPH)
        """    
    When C1 equals zero, the vapor pressure of NH3 should be zero.
    But, alas, the vapor pressure equation will return a non zero value
    due to the way the coefficients for the equation are calculated. To
    correct this, subtact off the fictous non-zero vapor pressure in a
    linear manner until at the lowest concentration for which there are
    non-interpolated/extrapolated coefficients.
        """
        if (C1 <= .05):
                GOLP0N = S7 + S8/T1 + S9*GOLT1
                P0N = exp(GOLP0N)
                X = (0.0500 - C1)/0.0500
                P0N = X*P0N
                VPN -= P0N
      
    #/*  Same problem as above except for H2O at high concentration.*/
        if (C1 >= .95):
                if (T1 <= 300.0 ):
                        GOLP0H = S10 + S11/T1 + S12 * GOLT1
                else:
                        GOLP0H = S13 + S14/T1 + S15 * GOLT1;
                P0H = exp(GOLP0H);
                X = (C1 - 0.9500)/0.0500;
                P0H = X*P0H;
                VPH -= P0H;
      
        return VPH,VPN

def SQ(C):
        return C**2


def solution_cloud(TC, PNH3, PH2O,others):#, float *SPNH3, float *SPH2O):
        """
        /***************************************************************************************************************
        solution_cloud(): This routine calculates the mole fraction of ammonia that goes into solution for the water-solution cloud.

         Input: 
              -->TC: Temperature
              -->PNH3: Partial Pressure of NH3 in bars.
              -->PH2O: Partial Pressure of H2O in bars.
              

         Output:
              <--C_sol: mole fraction of ammonia that dissolves in the solution.
              <--SPNH3: Saturation Pressure of Ammonia
              <--SPH2O: Saturation Pressure of Water
        **************************************************************************************************************/
        """
        from scipy import interpolate
        from modelParams import CONVT
        from math import log,exp  
        ROMANI=True
        UPPER_CONC=0.5
        LOWER_CONC=0.0
        CONC_STEP=0.0001
        
        ammoniaA=interpolate.InterpolatedUnivariateSpline(others['x1'],others['a1'])
        ammoniaB=interpolate.InterpolatedUnivariateSpline(others['x1'],others['b1'])
        ammoniaC=interpolate.InterpolatedUnivariateSpline(others['x1'],others['c1'])
        
        waterLowTempA=interpolate.InterpolatedUnivariateSpline(others['x2'],others['a2'])
        waterLowTempB=interpolate.InterpolatedUnivariateSpline(others['x2'],others['b2'])
        waterLowTempC=interpolate.InterpolatedUnivariateSpline(others['x2'],others['c2'])
        
        waterHighTempA=interpolate.InterpolatedUnivariateSpline(others['x3'],others['a3'])
        waterHighTempB=interpolate.InterpolatedUnivariateSpline(others['x3'],others['b3'])
        waterHighTempC=interpolate.InterpolatedUnivariateSpline(others['x3'],others['c3'])
        
        S7=ammoniaA(0.0)
        S8=ammoniaB(0.0)
        S9=ammoniaC(0.0)
        
        S10=waterLowTempA(1.0)
        S11=waterLowTempB(1.0)
        S12=waterLowTempC(1.0)
        
        S13=waterHighTempA(1.0)
        S14=waterHighTempB(1.0)
        S15=waterHighTempC(1.0)
        if (ROMANI):  #/* Is there a prayer of a solution cloud? */
                
                S1=ammoniaA(UPPER_CONC)
                S2=ammoniaB(UPPER_CONC)
                S3=ammoniaC(UPPER_CONC)
                      
                if (TC < 300.0):
                        S4=waterLowTempA(UPPER_CONC)
                        S5=waterLowTempB(UPPER_CONC)
                        S6=waterLowTempC(UPPER_CONC)
                else:
                        S4=waterHighTempA(UPPER_CONC)
                        S5=waterHighTempB(UPPER_CONC)
                        S6=waterHighTempC(UPPER_CONC)
                        
                SPH2O,SPNH3=PPRES(UPPER_CONC, TC,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15)
                SPH2O *= CONVT
                SPNH3 *= CONVT
                if (SPH2O > PH2O):  #/* No.  Then return. */
                        return -1.0,SPH2O,SPNH3
      

        SPNH3 = 0.0
        PNH3C = PNH3
        C_sol=LOWER_CONC
        have_cloud=False
        #for (C_sol=LOWER_CONC; *SPNH3<PNH3C && C_sol<UPPER_CONC; C_sol+=CONC_STEP):
        while(C_sol<UPPER_CONC and SPNH3<PNH3C):        
                PNH3C = PNH3*(1.0-C_sol);
                if (ROMANI):#/*  Romani version  */
                        S1=ammoniaA(C_sol)
                        S2=ammoniaB(C_sol)
                        S3=ammoniaC(C_sol)

                        if (TC < 300.0): # low temperature fit for H2O (where the LT HT comes from in NWLT NWHT N-number of points W-water L-Low T-temperature)
                                S4=waterLowTempA(C_sol)
                                S5=waterLowTempB(C_sol)
                                S6=waterLowTempC(C_sol)
                                
                        else:
                                S4=waterHighTempA(C_sol)
                                S5=waterHighTempB(C_sol)
                                S6=waterHighTempC(C_sol)
                                
                        SPH2O,SPNH3=PPRES(C_sol, TC,S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15)
                        SPH2O *= CONVT
                        SPNH3 *= CONVT
            
                elif (not ROMANI):#/*  Briggs and Sackett */
                        SPH2O = log(1.0-C_sol) + 29.0423 + 4.0134*SQ(C_sol) - (5540.48 + 2022.11*SQ(C_sol))/TC
                        SPH2O = exp(SPH2O) * CONVT
                        SPNH3 = log(C_sol) + 30.0048 + 4.0134*(SQ(C_sol) - 2.0*C_sol) - (4949.75 + 2022.11*(SQ(C_sol) - 2.0*C_sol))/TC
                        SPNH3 = exp(SPNH3) * CONVT
            

                if (PH2O > SPH2O and PNH3C > SPNH3): #  We have a solution cloud 
                          have_cloud = True              #   but what C?           
                if (have_cloud  and  PNH3C < SPNH3):                
                        SPNH3 = PNH3C;
		        #SPH2S = h2s_sol_cloud(TC, C_sol,PH2S);  #/*DON'T Ignore H2S in cloud */
                        return C_sol,SPH2O,SPNH3
                    
                C_sol+=CONC_STEP      
        return -1.0,SPH2O,SPNH3




