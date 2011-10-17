
def get_dT( j,  T,  P,  dP, LX, L2X, Hydrogen_Curve_Fit_Select,layer,others):
        """
         float get_dT()  This function uses Pressure temperature and the latent heats associated with
                 each layer and calculates the temperature step associated with a step in
                 pressure.

          Input:
            --> int j: index associated with layer data structure (see model.h)
            --> float T: Temperature of layer in K
            --> float P: Pressure of the layer in bars
            --> float dP: Pressure step in bars
            --> *LX: Array of latent heats times mole fraction (see DeBoer eqn 3.19)
            --> *L2X: Array of latent heats squared times mole fraction (see DeBoer eqn 3.19)

        """
        from numpy import asarray,nonzero
        from modelParams import R
        from specific_heat import specific_heat
        
        PfL=others['PfL']
        TfL=others['TfL']
        if (P<=max(PfL)):   #/*linear interpolation from Lindal's points*/
                arrPfL=asarray(PfL)
                arrTfL=asarray(TfL)
                idx,=nonzero(PfL==P)
                dT=layer['T'][j-1]-float(arrTfL[idx])
                layer['T'][j]=layer['T'][j-1]+dT       
        else:
                Cp = specific_heat(j, T, P,Hydrogen_Curve_Fit_Select)
	        dT_num = (R*T + LX[0]+LX[1]+LX[2]+LX[3]+LX[4]+LX[5])*dP
                dT_den = P*(Cp + L2X[0]+L2X[1]+L2X[2]+L2X[3]+L2X[4]+L2X[5])
                dT = dT_num/dT_den                                         #eqn 3.19 in DeBoer's thesis
                layer[T][j] = layer['T'][j-1] + dT
        return dT,layer


