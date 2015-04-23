
import numpy
from nh3.HanleySteffes import HanleySteffes
from nh3.DevarajSteffesDuong import DevarajSteffesDuong
from h2s.DeBoerSteffes import DeBoerSteffes
from h2o.KarpowiczSteffes import KarpowiczSteffes
from h2o.KarpowiczSteffesModified import KarpowiczSteffesModified
from h2o.DeBoerCorrected import DeBoerCorrected
from h2o.Goodman import Goodman
from ph3.HoffmanSteffes import HoffmanSteffes
from clouds.get_complex_dielectric_constant_water import get_complex_dielectric_constant_water
from clouds.rayleigh_absorption import rayleigh_absorption
from h2.OrtonCIA import OrtonCIA

def findkappa(f,T,P,P_H2,P_He,P_NH3,P_H2O,P_CH4,P_PH3,P_H2S,XH2,XHe,XNH3,XH2O,DNH4SH,DH2S,DNH3,DH2O,DSOL,select_h2h2_model,select_ammonia_model,select_water_model,include_clouds):
    """
         function kappa=findkappa(f,T,P_H2,P_He,P_NH3,P_H2O,P_CH4,P_PH3)
     Finds the atmospheric absorption at each pressure level.
    
                            -->INPUT:
                                   ->f: Frequency in GHz
                                   ->T; Temperature Profile (K)
                                   ->P: Pressure Profile (bars)
                                   ->P_H2: Partial pressure of hydrogen/H2 (bars)
                                   ->P_He: Partial pressure of helium (bars)
                                   ->P_NH3: Partial pressure of ammonia (bars)
                                   ->P_H2S: Partial pressure of hydrogen sulfide (bars)
                                   ->XH2: absolute Mole fraction of hydrogen/H2 
                                   ->XHe: absolute Mole fraction of He
                                   ->XNH3: absolute Mole fraction ammonia
                                   ->XH2O: absolute mole fraction of water
                                   ->DNH4SH: cloud density of ammonium hydrosulfide cloud (g/cm^3)
                                   ->DNH3: cloud density of ammonia ice (g/cm^3)
                                   ->DH2O: cloud density of water ice (g/cm^3)
                                   ->DSOL: cloud density of water ammonia solution (g/cm^3)
                                   ->select_h2h2_model: select hydrogen continuum absorption model (1=Joiner Model, 2=Goodman, 3=Goodman by Joiner, 4=Borysow, 5=Borysow Orton modification
                                   ->select_ammonia_model: select ammonia absorption model (1=Spilker original Hoffman coding,2=Joiner Steffes 3=Berge Gulkis 4=Mohammed Steffes KA_band 5=Spilker Model
                                   ->select_water_model: select water absorption model (1=DeBoer Original, 2=DeBoer Corrected, 3=Goodman)
                                   ->include_clouds: select whether or not to include cloud absorption (1=Yes, 0=No)
                            <--OUTPUT:
                                   <-kappa:layer(s) absorption coefficient in cm^-1
     Output shoud be in optical depth per centimeter
    """ 

    # THIS VERSION USED FOR IMAGING DOES ALL PRESSURE LEVELS-TO BE USED MULTIPLE TIMES FOR
    # THE DIFFERENT LOOK ANGLES

    # Includes the opcacity due to Ammonia, Phosphine, Hydrogen Sulfide, Hydrogen, Water Vapor
    # inlcuding H2 collisions with CH4
    

    OpticaldepthstodB=434294.5
    stopindex=T.shape[0];
    stopindex=stopindex-1;
     # can't include bottom P pressure..no appropriate depth associated with it.
    # Call ammonia
    farr=numpy.array([[f]])
    
    alphanh3=numpy.zeros([stopindex])
    kappa_1=numpy.zeros([stopindex])
    
        
    for k in range(0,stopindex):
            if(select_ammonia_model==1):
                alphanh3_t=HanleySteffes(farr,T[k],P[k],XH2[k],XHe[k],XNH3[k])/OpticaldepthstodB
                alphanh3[k]=alphanh3_t[0]
            elif(select_ammonia_model==2):
                
                alphanh3_t=DevarajSteffesDuong(farr,T[k],P[k],XH2[k],XHe[k],XH2O[k],XNH3[k])/OpticaldepthstodB
                alphanh3[k]=alphanh3_t[0]
            else: print "sorry nothing other than Hanley-Steffes and Devaraj-Steffes-Duong Right now. %d is an invalid selection 1 or 2 are your options."%select_ammonia_model 
    


    # Call phosphine
    
    alphaph3=numpy.zeros([stopindex])
    for k in range(0,stopindex):
           alphaph3_t=HoffmanSteffes(farr,T[k],P_H2[k],P_He[k],P_PH3[k])/OpticaldepthstodB    # Vector in f
           alphaph3[k]=alphaph3_t[0]
    
    
    # Call Hydrogen Sulfide
    
    alphah2s=numpy.zeros([stopindex])
    for k in range(0,stopindex):
           alphah2s_t=DeBoerSteffes(farr,T[k],P_H2[k],P_He[k],P_H2S[k])/OpticaldepthstodB    # Vector in f
           alphah2s[k]=alphah2s_t[0]
    
    # Call Hydrogen
    
    alphah2=numpy.zeros([stopindex])
    for k in range(0,stopindex):
            
        alphah2_t=OrtonCIA(farr,T[k],P_H2[k],P_He[k],P_CH4[k])
        alphah2[k]=alphah2_t[0]
    
    alphah2o=numpy.zeros([stopindex])
    for k in range(0,stopindex):
        if(select_water_model==1):
            alphah2o_t=DeBoerCorrected(farr,T[k],P_H2[k],P_He[k],P_H2O[k])/OpticaldepthstodB
            alphah2o[k]=alphah2o_t[0]
        if(select_water_model==2):
            alphah2o_t=Goodman(farr,P[k],T[k],XH2O[k],XH2[k],XHe[k])/OpticaldepthstodB
            alphah2o[k]=alphah2o_t[0]
        if(select_water_model==3):
            density_h2_g_m3=(P_H2[k]*2.01594)/((8.314472e-5)*T[k]);
            density_he_g_m3=(P_He[k]*4.0026)/((8.314310e-5)*T[k]);
            density_h2o_g_m3=(P_H2O[k]*(8.314472/0.46141805))/((8.314472e-5)*T[k]);
            alphah2o_t=KarpowiczSteffes(farr,density_h2o_g_m3,density_h2_g_m3,density_he_g_m3,T[k])
            alphah2o[k]=alphah2o_t[0]/OpticaldepthstodB
        if(select_water_model==4):
            density_h2_g_m3=(P_H2[k]*2.01594)/((8.314472e-5)*T[k]);
            density_he_g_m3=(P_He[k]*4.0026)/((8.314310e-5)*T[k]);
            density_h2o_g_m3=(P_H2O[k]*(8.314472/0.46141805))/((8.314472e-5)*T[k]);
            alphah2o_t=KarpowiczSteffesModified(farr,density_h2o_g_m3,density_h2_g_m3,density_he_g_m3,T[k])
            
    complex_dielectric_SOL=numpy.zeros([stopindex],dtype=numpy.complex)
    complex_dielectric_SOL=numpy.zeros([stopindex],dtype=numpy.complex)
    complex_dielectric_H2S=numpy.zeros([stopindex],dtype=numpy.complex)
    complex_dielectric_NH3=numpy.zeros([stopindex],dtype=numpy.complex)
    complex_dielectric_H2O=numpy.zeros([stopindex],dtype=numpy.complex)
    complex_dielectric_NH4SH=numpy.zeros([stopindex],dtype=numpy.complex)
    
    alpha_NH4SH=numpy.zeros([stopindex])
    alpha_H2S=numpy.zeros([stopindex])
    alpha_NH3=numpy.zeros([stopindex])
    alpha_H2O=numpy.zeros([stopindex])
    alpha_SOL=numpy.zeros([stopindex])
        
    if(include_clouds==1):
            
            for k in range(0,stopindex):
                complex_dielectric_SOL[k]=get_complex_dielectric_constant_water(f,T[k]); #Use equations in Ulaby, Fung, Moore for water
                complex_dielectric_H2S[k]=0.0 #Don't care no H2S clouds here.
                complex_dielectric_NH3[k]=(((1.48)**2)-((8.73e-4)**2)) + 1j*(2*1.48*8.73e-4); #Use Howett, Carlson, Irwin,Calcutt ref index nu=1300 cm^-1
                complex_dielectric_H2O[k]=3.15 + 1j*1e-3;                                   #Worst Case at 33 GHz Matsuoka,Fujita,Mae 1996
                complex_dielectric_NH4SH[k]=(((2.72)**2)-((7.83e-4)**2)) + 1j*(2*2.72*7.83e-4);#Use Howett, Carlson, Irwin,Calcutt ref index nu=1300 cm^-1
                
                alpha_NH4SH[k]=rayleigh_absorption(f,DNH4SH[k],1.,complex_dielectric_NH4SH[k]);
                alpha_H2S[k]=rayleigh_absorption(f,DH2S[k],1.,complex_dielectric_H2S[k]);
                alpha_NH3[k]=rayleigh_absorption(f,DNH3[k],1.,complex_dielectric_NH3[k]);
                alpha_H2O[k]=rayleigh_absorption(f,DH2O[k],1.,complex_dielectric_H2O[k]);
                alpha_SOL[k]=rayleigh_absorption(f,DSOL[k],1.,complex_dielectric_SOL[k]);
            for k in range(0,stopindex):
                kappa_1[k]=alphanh3[k]+alphaph3[k]+alphah2o[k]+alphah2s[k]+alpha_NH4SH[k]+alpha_H2S[k]+alpha_NH3[k]+alpha_H2O[k]+alpha_SOL[k]+alphah2[k]
            
    else:
        for k in range(0,stopindex):
            kappa_1[k]=alphanh3[k]+alphaph3[k]+alphah2o[k]+alphah2s[k]+alphah2[k]
            
    

    #save alfs
        
    kappa=numpy.transpose(kappa_1);    # (in 1/cm) if limb-have passed same pressure twice-masterindex keeps track
        
    return kappa
