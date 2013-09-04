# Go crazy! run it a whole bunch of times!
#
from TCM.DeBoerTCM import DeBoerTCM
from rtm.maintamone import maintamone
from scipy.io import savemat
import numpy
import scipy.interpolate
olderr = numpy.seterr(all='ignore')
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
inputPar['dP_init']=10.0 
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

################################
# Thermo-chemical Modeling
################################
me,tcme,tcmp,DSOL_NH31=DeBoerTCM(TP_list,TP_force,oblateness_factor,inputPar) 
################################
# Set Spacecraft Orientation
theta=0
Raydirection[0]=-1.0*numpy.cos(theta*(numpy.pi/180.)) 
Raydirection[1]=numpy.sin(theta*(numpy.pi/180.)) 

Tbeam_nadir=[]
zenith_nadir=[]
weighting_function_a_nadir=[]
refractive_index=[]
#Run Radiative Transfer model for all frequencies
for j in range(0,len(f)):
    no_ph3=0  
    [Tbeam_nadir_t,zenith_nadir_t,weighting_function_a_nadir_t,refractive_index_t,Tatma,intercepts_boresight,intercepts_b]= maintamone(Raydirection,Rayorigin,
                                                                                   tcme,tcmp,ao,bo,co,f[j],no_ph3,select_h2h2_model,select_ammonia_model,
                                                                                   select_water_model,include_clouds,N_ring_one,Nphi,BWHM,refractivity_source,
                                                                                   cassini_pattern,cassini_data_path) 
    Tbeam_nadir.append(Tbeam_nadir_t)
    zenith_nadir.append(zenith_nadir_t)
    weighting_function_a_nadir.append(weighting_function_a_nadir_t)
    refractive_index.append(refractive_index_t)
    print 'Tatma',Tatma
    print 'Tbeam_nadir', Tbeam_nadir
    """
    a={}
    a['Tbeam_nadir']=Tbeam_nadir
    a['tcme']=tcme
    a['weighting_function_a_nadir']=weighting_function_a_nadir
    savemat(output_filename,{'a':a}) 
    """


# Set Spacecraft Orientation

theta=54.5910 #along z
#theta=54.5912 #along y
#theta=60 
#Set Viewing Geometry

X_direction=-numpy.cos(theta*(numpy.pi/180)) 
Y_direction=0. 
Z_direction=numpy.sin(theta*(numpy.pi/180)) 
Raydirection=numpy.array([X_direction, Y_direction, Z_direction]) 

#Run Radiative Transfer Model For all Frequencies.
Tbeam_limb=[]
zenith_limb=[]
weighting_function_a_limb=[]

for j in range(0,len(f)):
    no_ph3=0 
    Tbeam_limb_t,zenith_limb_t,weighting_function_a_limb_t,refractive_index_t,Tatma,intercepts_boresight,intercepts_b= maintamone(Raydirection,Rayorigin,
                                                                                tcme,tcmp,ao,bo,co,f[j],no_ph3,select_h2h2_model,select_ammonia_model,
                                                                                select_water_model,include_clouds,N_ring_one,Nphi,BWHM,refractivity_source,
                                                                                cassini_pattern,cassini_data_path) 
                                                                                
    
    Tbeam_limb.append(Tbeam_limb_t)
    zenith_limb.append(zenith_limb_t)
    weighting_function_a_limb.append(weighting_function_a_limb_t)
    
    
    
    """
    a['Tbeam_limb']=Tbeam_limb    
    a['weighting_function_a_limb']=weighting_function_a_limb 
    
    print 'Tbeam_limb',Tbeam_limb
   
    savemat(output_filename,{'a':a})  
    """    

Tbeam_nadir=numpy.asarray(Tbeam_nadir)
zenith_nadir=numpy.asarray(zenith_nadir)
weighting_function_a_nadir=numpy.asarray(weighting_function_a_nadir)
refractive_index=numpy.asarray(refractive_index)
    
Tbeam_limb=numpy.asarray(Tbeam_limb)
zenith_limb=numpy.asarray(zenith_limb)
weighting_function_a_limb=numpy.asarray(weighting_function_a_limb)

a={}
R=100.*(Tbeam_nadir-Tbeam_limb)/Tbeam_nadir 
a['R']=R
savemat(output_filename,{'a':a}) 
numpy.seterr(*olderr) 
