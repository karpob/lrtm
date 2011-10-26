
#def specific_heat(int j, float T, float P,float Hydrogen_Curve_Fit_Select);
        
def cloud_loss_ackerman_marley( j,Teff, T,  P, H,  wet_adiabatic_lapse_rate,\
                                dry_adiabatic_lapse_rate, current_z,  previous_z, \
                                previous_q_c, current_q_v, previous_q_v, XH2, \
                                XHe, XH2S, XNH3, XH2O, XCH4, XPH3, \
                                delta_q_c, frain, Hydrogen_Curve_Fit_Select,layer,first):
        """
        /******************************************************************************************/
        // float cloud_loss_ackerman_marley()
        //                                   This function computes the mole fraction of condensate
        //                                   that remains in a cloud layer when considering vertical 
        //                                   Diffusion, and a paramter frain which represents potential
        //                                   sources of wet/dry deposition, following Ackerman & Marley, 2001.
        //                  Input:
        //                     --> int j: the index of layer data structure seem (model.h)
        //                     --> float Teff: Effective temperature associated with Flux (in Kelvin)
        //                     --> float P: pressure in bars.
        //                     --> float H: scale height in cm.
        //                     --> float wet_adiabatic_lapse_rate: moist adiabatic lapse rate (K/km)
        //                     --> float dry_adiabatic_lapse_rate: dry adiabatic lapse rate (K/km)
        //                     --> float current_z: Current altitude (km)
        //                     --> float previous_z: pressure at the previous step (j-1)
        //                     --> float previous q_c: previous value for mole fraction of condensate
        //                     --> float current_q_v: current value for constituent in the gas phase (mole fraction)
        //                     --> float previous_q_v: previous value (j-1) for constituent in the gas phase (mole fraction)
        //                     --> float XH2: mole fraction of hydrogen
        //                     --> float XHe: mole fraction of helium
        //                     --> float XH2S: mole fraction of hydrogen sulfide
        //                     --> float XNH3: mole fraction of ammonia
        //                     --> float XH2O: mole fractoin of water
        //                     --> float XCH4: mole fraction of methane
        //                     --> float XPH3: mole fraction of phosphine
        //                     --> delta_q_c: change in condensate (initial condition)
        //                     --> float frain: value for f_rain as stated in Ackerman & Marley, 2001.                
        /******************************************************************************************/
        """
                        

        from specific_heat import specific_heat        
        from modelParams import R,AMU_H2,AMU_He,AMU_H2S,AMU_NH3,AMU_H2O ,AMU_CH4,AMU_PH3,AMU_NH4SH,TRIPLEPT_H2S,TRIPLEPT_NH3,TRIPLEPT_CH4,TRIPLEPT_PH3,TRIPLEPT_H2O,GONE,ZERO                            
        boltz_sigma=5.6704e-8                 # W/m^2/K^4 (stefan boltzmann constant)
        Flux=boltz_sigma*(pow(Teff,4))         # W/m^2  (Jupiter's thermal emission)
        
        cp_temp=specific_heat(j,T,P,Hydrogen_Curve_Fit_Select,layer);  
        
        mu_temp=AMU_H2*XH2 + AMU_He*XHe + AMU_H2S*XH2S + AMU_NH3*XNH3 + AMU_H2O*XH2O + AMU_CH4*XCH4 + AMU_PH3*XPH3; #//  g/mol
                  
        if((wet_adiabatic_lapse_rate/dry_adiabatic_lapse_rate)>0.1): GAMMA_TEMP=(wet_adiabatic_lapse_rate/dry_adiabatic_lapse_rate);
        else: GAMMA_TEMP=0.1; #//Eqn 6 of Ackerman and Marley.

        mixing_L=H*GAMMA_TEMP;
        rho_temp=(mu_temp*P)/(R*T);
        #F=(GAMMA_TEMP-1.0)/gamma;
        F=1./3.
        #//eqn. 5 in Ackerman and Marley, 2001
        Eddy_Diffusion_Coef=(1./3.)*H*pow(mixing_L/H,(4.0/3.0))*pow((1e-3)*(R*Flux)/(mu_temp*rho_temp*cp_temp),1.0/3.0);
                                                                                                
                                                                                                #1e-3 is there to get consistent units.  
                                                                                                #Should be in cm^2/sec
        if(Eddy_Diffusion_Coef<1e5): Eddy_Diffusion_Coef=1e5; #Lower limit on Diffusion as suggested by Ackerman and Marley, 2001                                     
        my_dz=1e5*(layer['z'][j]-layer['z'][j-1]); # my_dz in cm
        wstar=Eddy_Diffusion_Coef/mixing_L;  # cm/sec
        #flag for being past first layer
        if(first!=1):
                q_c_old=delta_q_c
                q_v_old=previous_q_v
        else:
                q_c_old=previous_q_c; 
                q_v_old=previous_q_v;                    
         
        q_v=current_q_v;
        q_c=(Eddy_Diffusion_Coef*(q_v-q_v_old-q_c_old))/(Eddy_Diffusion_Coef+frain*wstar*my_dz) #//solution in discrete form of eqn. 4 in Ackerman and Marley
                  
        return q_c

