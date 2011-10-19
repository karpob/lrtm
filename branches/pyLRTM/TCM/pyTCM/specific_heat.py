

def specific_heat(j, T, P,Hydrogen_Curve_Fit_Select,layer):
        """

        specific_heat() this function calculates the specific heat for a 
                  given temperature and pressure.

              Input: 
                  -->int j : layer index.
                  -->float T: Temperature
                  -->float P: Pressure in bars  
              Output:
                  <--Cp: specific heat in erg/K/mol 
                     To select equilibrium normal or intermediate this is switched through
                     the global variable Hydrogen_Curve_Fit_select 
                                       0=equilibrium, -1=intermediate,0.25=normal


        """
        from modelParams import R
        from python_compressibility.calc_Cp.getCp import getCp
        #      fit for equilibrium   
      	if (T<120.0 and Hydrogen_Curve_Fit_Select==0.0):
      	        Cp_H2 = 8.8477 - 1.0421*T + 6.1631e-2*T**2 - 1.6961e-3*pow(T,3.0) + \
      	                2.6360e-5*pow(T,4.0) - 2.5341e-7*pow(T,5.0) + 1.5665e-9*pow(T,6.0) - \
                        6.2470e-12*pow(T,7.0) + 1.5544e-14*pow(T,8.0) - 2.19549e-17*pow(T,9.0) + \
                        1.344336e-20*pow(T,10.0)
                    
      	elif (T>=120.0 and T<=130.0 and Hydrogen_Curve_Fit_Select==0.0):
      	        Cp_H2 = -6.875E-3*(T-120.0) + 3.31494
      	elif (T>130.0 and T<=358.0 and Hydrogen_Curve_Fit_Select==0.0):
      	        Cp_H2 = 3.0795 + 0.001174*T
      	elif (T > 358.0 and Hydrogen_Curve_Fit_Select==0.0):
      	        Cp_H2 = 3.5
       

        # fit for intermediate or frozen equilibrium   
        if (T<=314.0 and Hydrogen_Curve_Fit_Select==-1.0):
                Cp_H2 = 2.6114 - 9.3144e-3*T + 2.0612e-4*T**2 - 1.217e-6*pow(T,3.0) +\
                        3.135e-9*pow(T,4.0) - 3.0537e-12*pow(T,5.0)
      	elif (Hydrogen_Curve_Fit_Select==-1.0):
      	        Cp_H2 = 3.5
            
        #      fit for normal hydrogen      ***/
	if (T<=310.0 and Hydrogen_Curve_Fit_Select==0.25):
           	Cp_H2 = 2.5785 - 5.6608e-3*T + 9.4152e-5*T**2 - 1.8306e-7*pow(T,3.0) -\
                   	6.239e-10*pow(T,4.0) + 1.696e-12*pow(T,5.0)
        elif (Hydrogen_Curve_Fit_Select==0.25):
            	Cp_H2 = 3.5
	
	if(Hydrogen_Curve_Fit_Select==666.0):
	        XH2  = layer['XH2'][j-1]
                XHe  = layer['XHe'][j-1]
                XH2S = layer['XH2S'][j-1]
                XNH3 = layer['XNH3'][j-1]
                XH2O = layer['XH2O'][j-1]
                XCH4 = layer['XCH4'][j-1]
                XPH3 = layer['XPH3'][j-1]

		Cp=getCp(T,XH2*P,XHe*P,XCH4*P,XH2O*P)
                
		
		Cp_H2S = 4.013
        	Cp_NH3 = 4.459  #/*Briggs and Sackett*/
		Cp=Cp+XH2S*Cp_H2S*R+XNH3*Cp_NH3*R
                

	else:
        	Cp_He = 2.503   
        	Cp_H2S = 4.013
        	Cp_NH3 = 4.459  #/*Briggs and Sackett*/
        	Cp_H2O = 4.0 
        	Cp_CH4 = 4.5    #/*                  */
		Cp_PH3 = 0.0
	        XH2  = layer['XH2'][j-1]
                XHe  = layer['XHe'][j-1]
                XH2S = layer['XH2S'][j-1]
                XNH3 = layer['XNH3'][j-1]
                XH2O = layer['XH2O'][j-1]
                XCH4 = layer['XCH4'][j-1]
                XPH3 = layer['XPH3'][j-1]
                Cp = (XH2*Cp_H2 + XHe*Cp_He + XH2S*Cp_H2S + XNH3*Cp_NH3 + XPH3*Cp_PH3)*R

        return Cp 

