from numpy import interp
import numpy as np
import os
from pandas.io.parsers import read_table

def DevarajSteffesDuong(f,T,P,H2mr,Hemr,H2Omr,NH3mr):
    """
    THIS FUNCTION COMPUTES THE OPACITY OF AMMONIA BROADENED BY HYDROGEN, 
    HELIUM, AND WATER VAPOR UNDER JOVIAN CONDITIONS BETWEEN 1-150 GHz
    NAME:
      NH3H2HeH2OModel
    EXPLANATION:
        This function can be used to compute the opacity of pure ammonia and
        that of ammonia in a hydrogen/helium/water vapor atmosphere between
        1-150 GHz under jovian conditions.
    
    CALLING SEQUENCE:
       alphanh3=NH3H2HeH2O(f,T,P,H2mr,Hemr,H2Omr,NH3mr)

    INPUTS:
       f        - Array of frequencies (GHz)
       T        - Temperature (K)
       P        - Pressure (bars)
       H2mr     - Hydrogen mixing ratio (as mole fraction)
       Hemr     - Helium mixing ratio (as mole fraction)
       H2Omr    - Water Vapor mixing ratio (as mole fraction)
       NH3mr    - Ammonia mixing ratio (as mole fraction)

    OUTPUTS:
      alphanh3  - Array of ammonia opacity (dB/km) computed at the input frequencies

    METHOD:
       A modified Ben Reuven (Ben Reuven, 1966) lineshape is used for computing
       ammonia opacity due to inversion lines, and a Gross lineshape (Gross, 1955)
       is used for computing ammonia opacity due to the rotational lines and
       the v2 roto-vibrational lines.

    History:
       written by Kiruthika Devaraj at Stanford,  April, 2014
       Ported to Python by Bryan Karpowicz in Boulder, March 2015. 
    """

    # The data files containing the frequency, line intensity and lower state 
    # energy for the ammonia transitions as given in the latest JPL spectral 
    # line catalog provided by Shanshan Yu and Brian Droiun  (personal 
    # communication, 2010), and the self and foreign gas broadening 
    # parameters for the various transitions as given by Devaraj et al. Icarus,
    # 2010 are loaded in the following steps. 

    # Inversion lines: fo is frequency in GHz, Io is line intensity in
    # cm^-1/(molecule./cm^2), Eo is lower state energy in cm^-1, gammaNH3o and
    # H2HeBroad are self and foreign gas broadening parameters.
    # old matlab -->[fo, Io, Eo, gammaNH3o, H2HeBroad] = textread('ammonia_inversion.dat','#f #f #f #f #f','headerlines',1);
    #pandas read_table to read in variable whitespace "stuff"
    a=read_table(os.path.join(os.getcwd(),'nh3LineParameters','ammonia_inversion.dat'),sep=r"\s*",header=0)
    fo,Io,Eo,gammaNH3o,H2HeBroad=np.asarray(a['fo']),np.asarray(a['Io']),np.asarray(a['Eo']),np.asarray(a['gammaNH3o']),np.asarray(a['H2HeBroad'])
    mm=fo.shape[0]
    fo=fo.reshape(mm,1)
    Io=Io.reshape(mm,1)
    Eo=Eo.reshape(mm,1)
    # Rotational lines: fo_rot is frequency in GHz, Io_rot is line intensity in
    # cm^-1/(molecule./cm^2), Eo_rot is lower state energy in cm^-1, gNH3_rot,
    # gH2_rot, gHe_rot are broadening parameters for rotational lines.
    #old matlab-->[fo_rot, Io_rot, Eo_rot, gNH3_rot, gH2_rot gHe_rot] = textread('ammonia_rotational.dat','#f #f #f #f #f #f','headerlines',1);
    a=read_table(os.path.join(os.getcwd(),'nh3LineParameters','ammonia_rotational.dat'),sep=r"\s*", header=0)
    fo_rot, Io_rot, Eo_rot, gNH3_rot, gH2_rot,gHe_rot=np.asarray(a['fo_rot']),np.asarray(a['Io_rot']), np.asarray(a['Eo_rot']), np.asarray(a['gNH3_rot']), np.asarray(a['gH2_rot']),np.asarray(a['gHe_rot'])    
    mm=fo_rot.shape[0]
    fo_rot=fo_rot.reshape(mm,1)
    Io_rot=Io_rot.reshape(mm,1)
    Eo_rot=Eo_rot.reshape(mm,1)
    gNH3_rot=gNH3_rot.reshape(mm,1)
    gH2_rot=gH2_rot.reshape(mm,1)
    gHe_rot=gHe_rot.reshape(mm,1)
    
    # v2 roto-vibrational lines: fo_v2 is frequency in GHz, Io_v2 is line intensity in
    # cm^-1/(molecule./cm^2), Eo_v2 is lower state energy in cm^-1,
    #old matlab --> [fo_v2, Io_v2, Eo_v2] = textread('ammonia_rotovibrational.dat','#f #f #f','headerlines',1);
    a=read_table(os.path.join('nh3LineParameters','ammonia_rotovibrational.dat'),sep=r"\s*", header=0)
    fo_v2,Io_v2,Eo_v2=np.asarray(a['fo_v2']),np.asarray(a['Io_v2']),np.asarray(a['Eo_v2'])
    mm=fo_v2.shape[0]
    fo_v2=fo_v2.reshape(mm,1)
    Io_v2=Io_v2.reshape(mm,1)
    Eo_v2=Eo_v2.reshape(mm,1)
    
    #clear out a in case it's big and hogs up RAM
    a=[]

    GHztoinv_cm=1.0/29.9792458           # for converting GHz to inverse cm
    OpticaldepthstodB=434294.5         # convert from cm^-1 to dB/km
    torrperatm=760                     # convert from atm to torr
    bartoatm=0.987                     # convert from bat to atm
    GHztoMHz=1000.0                    # convert from GHz to MHz
    hc=19.858252418e-24                #planks (J.s) light (cm/s)
    k=1.38e-23                    #boltzmann's in J/K or N.m/K
    No=6.02297e23                      #Avogadros Number [mole^-1]
    R=8.31432e7                        #Rydberg's [erg/mole-K]
    To=300.                         #Ref temp for P/P Catalogue
    dynesperbar=1.0e6                  #dyne=bar/1e6;
    coef=dynesperbar*No/R              #See Appendix D: Using the Poyter-Pickett Catalogs

    # Compute the partial pressure of H2, He, and NH2
    PH2=P*H2mr
    PHe=P*Hemr
    PH2O=P*H2Omr
    PNH3=P*NH3mr

    # Compute the temperature factor
    Tdiv=To/T

    #  Coefficient for symmetric top molecule
    eta=3./2.
    ## Computing the opacity due to the inversion lines

    # Pressure Dependent Switch for the parameters of the inversion transitions
    if P>20.0:
        gnu_H2=1.6361
        gnu_He=0.4555
        gnu_NH3=0.7298
        gnu_H2O=4.8901
        GAMMA_H2=0.8
        GAMMA_He=0.5
        GAMMA_NH3=1
        GAMMA_H2O=0.5
        zeta_H2=1.1313
        zeta_He=0.1
        zeta_NH3=0.5152
        zeta_H2O=2.7310
        Z_H2=0.6234
        Z_He=0.5
        Z_NH3=2./3.
        Z_H2O=1.
        d=0.2
        Con=1.3746
    elif P<=10:
        gnu_H2=1.7465
        gnu_He=0.9779
        gnu_NH3=0.7298
        gnu_H2O=4.8901
        GAMMA_H2=0.8202
        GAMMA_He=1.
        GAMMA_NH3=1.
        GAMMA_H2O=0.5
        zeta_H2=1.2163
        zeta_He=0.0291
        zeta_NH3=0.5152
        zeta_H2O=2.7310
        Z_H2=0.8873
        Z_He=0.8994
        Z_NH3=2./3.
        Z_H2O=1.0
        d=-0.0627
        Con=0.9862
    else:
        Pinterp=np.asarray([10,20])
        gnu_H2=interp(P,Pinterp,np.flipud(np.asarray([1.6361,1.7465])))
        gnu_He=interp(P,Pinterp,np.flipud(np.asarray([0.4555,0.9779])))
        gnu_NH3=0.7298
        gnu_H2O=4.8901
        GAMMA_H2=interp(P,Pinterp,np.flipud(np.asarray([0.8,0.8202])))
        GAMMA_He=interp(P,Pinterp,np.flipud(np.asarray([0.5,1.0])))
        GAMMA_NH3=1.0
        GAMMA_H2O=0.5
        zeta_H2=interp(P,Pinterp,np.flipud(np.asarray([1.1313,1.2163])))
        zeta_He=interp(P,Pinterp,np.flipud(np.asarray([0.1,0.0291])))
        zeta_NH3=0.5152
        zeta_H2O=2.7310
        Z_H2=interp(P,Pinterp,np.flipud(np.asarray([0.6234,0.8873])))
        Z_He=interp(P,Pinterp,np.flipud(np.asarray([0.5,0.8994])))
        Z_NH3=2./3
        Z_H2O=1.
        d=interp(P,Pinterp,np.flipud(np.asarray([0.2,-0.0627])))
        Con=interp(P,Pinterp,np.flipud(np.asarray([1.3746,0.9862])))
    # Individual broadening parameters
    gH2=gnu_H2*PH2
    gHe=gnu_He*PHe
    gH2O=gnu_H2O*PH2O
    gNH3=gnu_NH3*PNH3*gammaNH3o
    # Broadening parameter
    gamma=((gH2)*((Tdiv)**(GAMMA_H2))+(gHe)*((Tdiv)**(GAMMA_He))+(gH2O)*((Tdiv)**(GAMMA_H2O))+gNH3*(295.0/T)**(GAMMA_NH3))
    gamma=gamma.reshape(gamma.shape[0],1)
    # Shift parameter
    delt=d*gamma

    # Individual coupling parameters
    zH2=zeta_H2*PH2
    zHe=zeta_He*PHe
    zH2O=zeta_H2O*PH2O
    zNH3=zeta_NH3*PNH3*gammaNH3o
    # Coupling parameter
    zeta=(zH2)*((Tdiv)**(Z_H2))+(zHe)*((Tdiv)**(Z_He))+(zH2O)*((Tdiv)**(Z_H2O))+zNH3*(295.0/T)**(Z_NH3)
    zeta=zeta.reshape(zeta.shape[0],1)
    zetasize=fo.shape[0]
    pst=delt;                                   # answer in GHz
    #Coupling element, pressure shift and dnu or gamma are in GHz, need to convert brlineshape to inverse cm which is done below
    try:n=f.shape[1]  #returns the number of columns in f
    except:n=1
    m=fo.shape[0] #returns the number of rows in fo
    # f1 f2 f3 f4 ....fn  n times where n is the number of frequency steps
    # f1 f2 f3 f4 ....fn                in the observation range
    # ...
    # f1 f2 f3 f4 ....fn
    # m times where m is the number of spectral lines

    nones=np.ones([1,n])
    mones=np.ones([m,1])
     
    f_matrix=mones*f
    fo_matrix=fo.reshape(m,1)*nones

    # The 10^6 allows use of P(bar) for P(dynes/cm^2)
    
    expo=-(1./T-1./To)*Eo*hc/k
    ST=Io*np.exp(expo)    # S(T) =S(To)converted for temperature
    alpha_noshape=Con*coef*(PNH3/To)*((To/T)**(eta+2.0))*ST  
    
    #Alpha Max Found

    #Ben Reuven lineshape calculated by the brlineshape function gives the answer in GHz
    #Here we change from GHz to inverse cm.
    dnu_matrix=gamma*nones
    ce_matrix=zeta*nones
    pst_matrix=pst*nones
    Aa=(2./np.pi)*(np.power((f_matrix/fo_matrix),2))  
             
    Bb=(dnu_matrix-ce_matrix)*(np.power(f_matrix,2))
    Cc=dnu_matrix+ce_matrix
    Dd=(np.power((fo_matrix+pst_matrix),2)) + (np.power(dnu_matrix,2))-(np.power(ce_matrix,2))
    Ee=np.power(f_matrix,2)
    Jj=np.power((fo_matrix+pst_matrix),2)
    Gg=np.power(dnu_matrix,2)
    Hh=np.power(ce_matrix,2)
    Ii=4.*(np.power(f_matrix,2))*(np.power(dnu_matrix,2))
    Ff=Aa*(Bb+Cc*Dd)/((np.power((Ee-Jj-Gg+Hh),2))+Ii)

    Fbr=(1./GHztoinv_cm)*Ff
    alpha_noshape_matrix=alpha_noshape*nones
    br_alpha_matrix=alpha_noshape_matrix*Fbr

    ## Computing the opacity due to rotational lines
    # Computing the absorption contributed by the rotational lines
    ST_rot=Io_rot*(np.exp((1.0/To-1.0/T)*Eo_rot*hc/k))
    # Factor GAMMA:
    eta_H2=0.8730
    eta_He=2.0/3.0
    eta_NH3=1.0;
    # Factor nu:
    gnu_H2=0.2984
    gnu_He=0.75
    gnu_NH3=3.1789
    # Factor Con
    Con=2.4268

    # Individual broadening

    gH2=gnu_H2*PH2*gH2_rot
    gHe=gnu_He*PHe*gHe_rot
    gNH3=gnu_NH3*PNH3*gNH3_rot
    
    # Total broadening
    gamma_rot=((gH2)*((Tdiv)**(eta_H2))+(gHe)*((Tdiv)**(eta_He))+gNH3*(Tdiv)**(eta_NH3))
    

    try:n=f.shape[1]  #returns the number of columns in f
    except:n=1
    m=fo_rot.shape[0]
    nones=np.ones([1,n])
    mones=np.ones([m,1])
    f_matrix=mones*f
    fo_matrix=fo_rot*nones
    dnu_matrix=gamma_rot*nones

    # Gross Lineshape
    Aa=4.0/np.pi*(np.power(f_matrix,2))*dnu_matrix
    Bb=np.power(np.power(fo_matrix,2)-np.power(f_matrix,2),2)
    Cc=4.0*np.power(f_matrix,2)*np.power(dnu_matrix,2)
    
    F_rot=Aa/(Bb+Cc)
    Fbr_rot=(1.0/GHztoinv_cm)*F_rot;
    alpha_rot=Con*coef*(PNH3/To)*((To/T)**(eta+2.0))*ST_rot*nones*Fbr_rot

    ## Computing the opacity due to v2 roto-vibrational lines

    # Computing the absorption contributed by the v2 rotovibrational lines
    
    
    ST_v2=Io_v2*(np.exp((1.0/To-1.0/T)*Eo_v2*hc/k))
    
    # Broadening parameters for the v2 vibrational lines
    
    gH2_v2=1.4*np.ones([fo_v2.shape[0],1])
    gHe_v2=0.68*np.ones([fo_v2.shape[0],1])
    gNH3_v2=9.5*np.ones([fo_v2.shape[0],1])
    
    # Factor GAMMA:
    eta_H2=0.73
    eta_He=0.5716
    eta_NH3=1
    # Factor Con
    Con=1.1206

    # Individual broadening parameters
    gH2=PH2*gH2_v2
    gHe=PHe*gHe_v2
    gNH3=PNH3*gNH3_v2
    #Total broadening
    gamma_v2=((gH2)*((Tdiv)**(eta_H2))+(gHe)*((Tdiv)**(eta_He))+gNH3*(Tdiv)**(eta_NH3));
    
    #gamma_v2=gamma_v2.reshape(gamma_v2.shape[0],1)
    
    try:n=f.shape[1]  #returns the number of columns in f
    except:n=1
    m=fo_v2.shape[0]
    nones=np.ones([1,n])
    mones=np.ones([m,1])
    f_matrix=mones*f
    fo_matrix=fo_v2*nones
    dnu_matrix=gamma_v2*nones

    # Gross Lineshape
    Aa=4.0/np.pi*np.power(f_matrix,2)*dnu_matrix
    Bb=np.power(np.power(fo_matrix,2)-np.power(f_matrix,2),2)
    Cc=4.0*np.power(f_matrix,2)*np.power(dnu_matrix,2)
    
    F_v2=Aa/(Bb+Cc)
    Fbr_v2=(1.0/GHztoinv_cm)*F_v2
    alpha_v2=Con*coef*(PNH3/To)*((To/T)**(eta+2.0))*ST_v2*nones*Fbr_v2
    
    ## Computing the total opacity

    alpha_opdep=np.sum(br_alpha_matrix,axis=0)+np.sum(alpha_v2,axis=0)+np.sum(alpha_rot,axis=0)
    alphanh3=alpha_opdep*434294.5
    return alphanh3
if __name__ == '__main__':
    from scipy.io import loadmat
    import pylab
    o=loadmat(os.path.join('testData','testDataDevaraj.mat'))
    
    n=o['Freq'].shape[0]
    alphaCalculated=np.zeros(n)
    for i in range(0,n):
        
        alphaCalculated[i]=DevarajSteffesDuong(np.asarray(o['Freq'][i]).T,o['T'][i],o['P'][i],o['H2mr'][i],o['Hemr'][i],o['H2Omr'][i],o['NH3mr'][i])
    
    alphaCalculated=alphaCalculated.reshape(alphaCalculated.shape[0],1)
    alph=o['alphaCalculated'].reshape(o['alphaCalculated'].shape[1],1)
    pylab.plot(o['Freq'],100.0*(alphaCalculated-alph)/alph,'bx')
    pylab.xlabel('Frequency (GHz)')
    pylab.ylabel('% Difference mine vs. matlab')
    pylab.show()    

        
    
    