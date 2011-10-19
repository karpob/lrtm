
def get_dP_using_dP(j,eflag, inputParams,layer,others):
        """
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
        """                
        from modelParams import R
        from numpy import asarray,nonzero
        dP_init=inputParams['dP_init']
        dP_fine=inputParams['dP_fine']
        P_fine_start=inputParams['P_fine_start']
        P_fine_stop=inputParams['P_fine_stop'] 
        P_term=inputParams['P_term']
        P_targ=inputParams['P_targ']
        PfL=others['PfL']
        
        H = R*layer['T'][j-1]/(layer['mu'][j-1]*layer['g'][j-1])
        new_P_fine=layer['P'][j-1] - dP_fine
        new_P_coarse=layer['P'][j-1] - dP_init
        
        if((new_P_fine > P_fine_start) or (new_P_coarse > P_fine_start)):
                dP= -1.0*dP_init
                P= layer['P'][j-1] + dP
        elif((new_P_fine <= P_fine_start) and (new_P_fine > P_targ)):
                dP= -1.0*dP_fine
                P= layer['P'][j-1] + dP
        elif (new_P_fine<= P_targ and eflag!=97): #P=new_P_fine
                P = P_targ
                dP = PfL[1]-P_targ
                eflag=98
        elif(eflag==97):
                arrPfL=asarray(PfL)
                idx,=nonzero(arrPfL==layer['P'][j-1])
                P=float(arrPfL[idx+1])
                dP=float(arrPfL[idx+1]-layer['P'][j-1])
                #dP=float(arrPfL[idx])-layer['P'][j-1]
                #P=layer['P'][j-1]+dP
                
                #print j,eflag,P,dP,arrPfL[idx],idx,layer['P'][j-1]
                                        
        if(P <= P_term):
                eflag = 99
                P= P_term
                dP=P - layer['P'][j-1]
        
        dz= -1.0*dP/((layer['P'][j-1]/H)*(1e5))
        #print 'P',P_term
        layer['P'].append(P)
        layer['z'].append(layer['z'][j-1] + dz)
        
        return eflag,dP,layer


