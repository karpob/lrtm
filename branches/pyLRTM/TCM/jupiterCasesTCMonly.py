# Go crazy! run it a whole bunch of times!
#
import matplotlib
matplotlib.rc('text', usetex = True)
matplotlib.use('Agg')
from DeBoerTCM import DeBoerTCM
from scipy.io import savemat
import numpy
import scipy.interpolate
from pyTCM.cloudPlot import cloudPlot
from pyTCM.tpResidual import tpResidual
from pyTCM.constituentPlot import constituentPlot
from pyTCM.tpProfile import tpProfile
import pylab
#olderr = numpy.seterr(all='ignore')
oblateness_factor=0.935 

#oblateness_factor=1 
#####################################################
#Planet verticies                                   #
#####################################################
ao=71492*100000.  # along x (in centimeters)
bo=ao        # along y
co=ao*oblateness_factor  # along z

#select_h2h2_model
#6=orton et al based upon interpolation
select_h2h2_model=6 

#select_ammonia_model
#Hanley Steffes Model

select_ammonia_model=1

#select_water_model
#1 corrected deboer water vapor model
#2 (to be implemented goodman 1969 water vapor model
#3 Karpowicz and Steffes 2011, Icarus absorption model.
select_water_model=3

#include cloud absorption?
#1=yes
#0=no
include_clouds=1

# refractivity_source
# Select the author you believe is right with regards to values for refractivity (used for raypath calculations)
#
# refractivity_source=0  # No bending due to refraction n=1.0D0
 
refractivity_source=1  # Karpowicz w/Clouds H2, He, CH4 etc.. using Essen, and other sources
# #####################################################
#load craft 

#spacecraft parameters
Rayorigin=numpy.array([ao+(4500.*100000.), 0., 0.]) 
Raydirection=numpy.array([-1, 0., 0.]) 

#Get Beam parameters
N_ring_one=8  #Number of beams in first phi ring
Nphi=3  #Number of rings in phi
BWHM=10  # Beamwidth Half-maximum

#Cassini/Measured Antenna Pattern
cassini_pattern=0 
cassini_data_path='none' 

# Juno bands
#f=[0.6,1.2] #,2.4,4.8,9.6,23]  #operating frequency in GHz
#f=[0.5:0.1:10,10:1:25]  #too fine try something else..
f=numpy.r_[numpy.arange(0.6,1,0.1),numpy.arange(1.5,9.5,0.5),numpy.arange(10.,22.,2.),23] 
 
nfreq=len(f)
Selected_Model='Mean_Lindal'

Model_Names=['Mean_Lindal','Mean_Seiff','Depleted_Ammonia', 'Enhanced_Ammonia',
             'Depleted_Water','Enhanced_Water','Hot_Spot']



# Temperature forcing
TP_list=['Seiff_Jupiter','Lindal_Jupiter','Lindal_Saturn','whatever_is_in_TCM_mex']

if( Selected_Model==Model_Names[1] or Selected_Model==Model_Names[6] ):
    TP_force='Seiff_Jupiter' 
else:
    TP_force='Lindal_Jupiter' 



# Select a Filename
if(include_clouds==1):cloudz='_Clouds' 
else:cloudz='_No_Clouds' 

output_filename=Selected_Model+cloudz+'.mat' 


#Abundances
XH2S_rel_H2=7.7e-5 
XNH3_rel_H2=[4.55e-4,4.55e-4,2.0e-4,7.1e-4,4.55e-4,4.55e-4,2.3e-4] 
XH2O_rel_H2=[6.3838e-3,6.3838e-3,6.3838e-3,6.3838e-3,2.5535e-3,1.0214e-2,6e-4] 
XCH4_rel_H2=2.1e-3 
XPH3_rel_H2=6e-7 
XHe_rel_H2=0.157 
XH2=numpy.zeros([len(XNH3_rel_H2)])
XH2S=numpy.zeros([len(XNH3_rel_H2)]) 
XNH3=numpy.zeros([len(XNH3_rel_H2)])
XH2O=numpy.zeros([len(XNH3_rel_H2)])
XCH4=numpy.zeros([len(XNH3_rel_H2)])
XPH3=numpy.zeros([len(XNH3_rel_H2)])
XHe=numpy.zeros([len(XNH3_rel_H2)])
for i in range(0,len(XNH3_rel_H2)):
    XH2[i]=1/(1+XHe_rel_H2+XH2S_rel_H2+XNH3_rel_H2[i]+XH2O_rel_H2[i]+ XCH4_rel_H2+ XPH3_rel_H2) 

    #Get Mole Fractions
    XH2S[i]=XH2S_rel_H2*XH2[i] 
    XNH3[i]=XNH3_rel_H2[i]*XH2[i] 
    XH2O[i]=XH2O_rel_H2[i]*XH2[i] 
    XCH4[i]=XCH4_rel_H2*XH2[i] 
    XPH3[i]=XPH3_rel_H2*XH2[i] 
    XHe[i]=XHe_rel_H2*XH2[i] 


#Misc DeBoer TCM inputs
#Guess for deep P,T to match 1 bar level.
inputPar={}
inputPar['P_temp']=1000.  
inputPar['T_temp']=1583.782470703125 
inputPar['g0']=2330.  #2417 
inputPar['R0']=ao 
inputPar['P0']=1.0 
inputPar['T_targ']=165. 
inputPar['P_targ']=1. 
inputPar['P_term']=0.001 
inputPar['use_lindal']='Y' 
inputPar['SuperSatSelf_H2S']=0 
inputPar['SuperSatSelf_NH3']=0 
inputPar['SuperSatSelf_PH3']=0 
inputPar['SuperSatSelf_H2O']=0
inputPar['SuperSatSelf1']=0.0
inputPar['SuperSatSelf2']=0.0
inputPar['SuperSatSelf3']=0.0
inputPar['SuperSatSelf4']=0.0 
inputPar['supersatNH3']=0 
inputPar['supersatH2S']=0 
inputPar['AutoStep']=False
inputPar['AutoStep_constant']=8 
inputPar['fp']=0.25 
inputPar['dz']=1 
inputPar['use_dz']=False 
inputPar['dP_init']=0.05 
inputPar['dP_fine']=0.05 
inputPar['P_fine_start']=13 
inputPar['P_fine_stop']=1 
inputPar['frain']=3. 
inputPar['select_ackerman']=0 
inputPar['XCO']=0.0
#table_output=[XH2 XHe (1e6)*XH2S (1e6)*XNH3 (1e6)*XH2O (1e6)*XCH4 (1e6)*XPH3] 
#to_dlm=transpose(table_output) 
#dlmwrite('mole_fractions_jupiter.dat',to_dlm,'delimiter','&','precision','#.4f')

for i in range(0,len(XNH3)):
    if( Selected_Model==Model_Names[i]):
        case_select=i 
inputPar['XH2S']=XH2S[case_select]
inputPar['XNH3']=XNH3[case_select]
inputPar['XH2O']=XH2O[case_select]
inputPar['XCH4']=XCH4[case_select]
inputPar['XPH3']=XPH3[case_select]
inputPar['XHe']=XHe[case_select]



me,tcme1,tcmp,DSOL_NH31=DeBoerTCM(TP_list,TP_force,oblateness_factor,inputPar)



inputPar['fp']=666
################################
# Thermo-chemical Modeling
################################
me,tcme2,tcmp,DSOL_NH32=DeBoerTCM(TP_list,TP_force,oblateness_factor,inputPar)

inputPar['fp']=0.25
inputPar['P_term']=1.342e-3
inputPar['T_targ']=166.10
inputPar['P_targ']=1.0
me,tcme3,tcmp,DSOL_NH33=DeBoerTCM(TP_list,'Seiff_Jupiter',oblateness_factor,inputPar)
                                

inputPar['fp']=666
################################
# Thermo-chemical Modeling
################################
me,tcme4,tcmp,DSOL_NH34=DeBoerTCM(TP_list,'Seiff_Jupiter',oblateness_factor,inputPar)

                                                                 
################################
cloudPlot(tcme2,tcme1,'Cloud Densities','CloudsLindal.pdf')
constituentPlot(tcme2,tcme1,'Constituent Profiles','ConstituentsLindal.pdf')

cloudPlot(tcme4,tcme3,'Cloud Densities','CloudsSeiff.pdf')
constituentPlot(tcme4,tcme3,'Constituent Profiles','ConstituentsSeiff.pdf')

tpProfile(tcme1,tcme2,tcme3,tcme4,'Temperature Profiles','TP.pdf')
tpResidual(tcme1,tcme2,tcme3,tcme4,'$\Delta T$=$T_{updated}$ - $T_{original}$','$\Delta P$=$P_{updated}$ - $P_{original}$','Residual.pdf')
#print 'T, P, nh3, h2o, cloud'

