
def get_dP_using_dz(dz,j,eflag, inputParams,layer,others):

        """
        float get_dP_using_dz()
                         This function calculates the pressure step for a given value of dz
             Input:
                 --> int j: value of the layer index of layer data structure (see model.h)
                 --> int eflag: error flag which is returned to main 
                 --> dz : altitude step in km
             Output:
                 <-- dP: the pressure step in bars
                 <-- layer data structure see model.h

        """
        from modelParams import R
        AutoStep=inputParams['AutoStep']
        AutoStep_constant=inputParams['AutoStep_constant']
        P_term=inputParams['P_term']
        P_targ=inputParams['P_targ']
      
        if (AutoStep):
            #/*dz = log10(layer[j-1].P + 1.6); */
            #/* dz = log10(layer[j-1].P + 5.0); */
            dz=log10(layer[j-1].P + AutoStep_constant )  #what dave actually uses for Priscilla's stuff
            if (dz > 2.0): dz=2.0
      
        H = R*layer['T'][j-1]/(layer['mu'][j-1]*layer['g'][j-1])
        dP = -1.0*(layer['P'][j-1]/H)*(1.0e5*dz);
        P  = layer['P'][j-1] + dP;
        if (P <= P_term):
                eflag = 99
                P = P_term
                dP = P - layer['P'][j-1]
        elif ( P <= P_targ and eflag != 97 ):
                eflag= 98
                P = P_targ
                dP = P - layer['P'][j-1]
        
        layer['P'].append(P)
        layer['z'].append(layer['z'][j-1] + dz)
        
        return eflag,dP,layer
        

