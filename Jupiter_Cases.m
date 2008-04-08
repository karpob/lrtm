% Go crazy! run it a whole bunch of times!
%
clear all;
oblateness_factor=0.935;
%oblateness_factor=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Planet verticies                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ao=71492*100000; % along x (in centimeters)
bo=ao;       % along y
co=ao*oblateness_factor; % along z

%select_h2h2_model
%1=joiner
%2=goodman
%3=goodman by joiner
%4=borysow
%5=borysow with orton et al, 2007 modification
%6=orton et al based upon interpolation
select_h2h2_model=5;

%select_ammonia_model
%1 original hoffman coding of spilker
%2 toms code joiner-steffes
%3 " berge-gulkis
%4 " mohammed-steffes
%5 " spilker
%6 Hanley Steffes Model
% Note, 1 and 5 won't work for current Jupiter/adams data set
% spilker correction factor C goes negative, giving negative absorption
% coefficient

select_ammonia_model=6;

%select_water_model
%1 original deboer water vapor model
%2 corrected deboer water vapor model
%3 (to be implemented goodman 1969 water vapor model
select_water_model=2;

%include cloud absorption?
%1=yes
%0=no
include_clouds=0

% refractivity_source
% Select the author you believe is right with regards to values for refractivity (used for raypath calculations)
%
% refractivity_source=0; % No bending due to refraction n=1.0D0
% refractivity_source=1; % Original DeBoer/Hoffman H2/He refractivity 
% refractivity_source=2; % Karpowicz H2/He refractivity using original Essen data
% refractivity_source=3; % Karpowicz H2, He, CH4 etc.. using Essen, and other sources
 refractivity_source=4; % Karpowicz w/Clouds H2, He, CH4 etc.. using Essen, and other sources
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load craft;

%spacecraft parameters
Rayorigin=[ao+(4500*100000) 0 0];
Raydirection=[-1 0 0];

%Get Beam parameters
N_ring_one=8; %Number of beams in first phi ring
Nphi=3; %Number of rings in phi
BWHM=10; % Beamwidth Half-maximum

%Cassini/Measured Antenna Pattern
cassini_pattern=0;
cassini_data_path='none';

% Juno bands
%f=[0.6,1.2];%,2.4,4.8,9.6,23]; %operating frequency in GHz
f=[0.5:0.1:10,10:1:25];
%f=[0.5,0.6,1,25];
nfreq=length(f)
Selected_Model='Enhanced_Water'

Model_Names={'Mean_Lindal','Mean_Seiff','Depleted_Ammonia', 'Enhanced_Ammonia',...
             'Depleted_Water','Enhanced_Water','Hot_Spot'}



% Temperature forcing
TP_list={'Seiff_Jupiter','Lindal_Jupiter','Lindal_Saturn','whatever_is_in_TCM_mex'}
if(strcmp(Selected_Model,Model_Names(2))|strcmp(Selected_Model,Model_Names(7)))
    TP_force='Seiff_Jupiter';
else
    TP_force='Lindal_Jupiter';
end


% Select a Filename
if(include_clouds==1)
    cloudz='_Clouds';
else
    cloudz='_No_Clouds';
end
output_filename=strcat(Selected_Model,cloudz,'.mat');


%Abundances
XH2S_rel_H2=7.7e-5;
XNH3_rel_H2=[4.55e-4,4.55e-4,2.0e-4,7.1e-4,4.55e-4,4.55e-4,2.3e-4];
XH2O_rel_H2=[6.3838e-3,6.3838e-3,6.3838e-3,6.3838e-3,2.5535e-3,1.0214e-2,6e-4];
XCH4_rel_H2=2.1e-3;
XPH3_rel_H2=6e-7;
XHe_rel_H2=0.157;
for i=1:length(XNH3_rel_H2)
    XH2(i)=1/(1+XHe_rel_H2+XH2S_rel_H2+XNH3_rel_H2(i)+XH2O_rel_H2(i)+ XCH4_rel_H2+ XPH3_rel_H2);

    %Get Mole Fractions
    XH2S(i)=XH2S_rel_H2*XH2(i);
    XNH3(i)=XNH3_rel_H2(i)*XH2(i);
    XH2O(i)=XH2O_rel_H2(i)*XH2(i);
    XCH4(i)=XCH4_rel_H2*XH2(i);
    XPH3(i)=XPH3_rel_H2*XH2(i);
    XHe(i)=XHe_rel_H2*XH2(i);
end

%Misc DeBoer TCM inputs
%P_temp=1000;
%T_temp=1303.349;
P_temp=6000;
T_temp=2200;

g0_i=2330; %2417;
R0e_i=ao;
P0_i=1;
T_targ_i=166;
P_targ_i=1;
P_term_i=0.141;
use_lindal='Y';
SuperSatSelf_H2S=0;
SuperSatSelf_NH3=0;
SuperSatSelf_PH3=0;
SuperSatSelf_H2O=0;
supersatNH3=0;
supersatH2S=0;
AutoStep_constant=8;
fp=-1;
dz=1;
XCO=0;
use_dz=0;
dP_init=1;
dP_fine=0.1;
P_fine_start=10;
P_fine_stop=1;
frain=3;
select_ackerman=0;

 table_output=[XH2;XHe;(1e6)*XH2S;(1e6)*XNH3;(1e6)*XH2O;(1e6)*XCH4;(1e6)*XPH3];
    to_dlm=transpose(table_output);
    dlmwrite('mole_fractions_jupiter.dat',to_dlm,'delimiter','&','precision','%.4f')

for i=1:length(XNH3)
    if(strcmp(Selected_Model,Model_Names(i)))
        case_select=i;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thermo-chemical Modeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[me,tcme,tcmp,DSOL_NH3]=DeBoer_TCM(TP_list,TP_force,XH2S(case_select),XNH3(case_select),XH2O(case_select),XCH4(case_select),...
                                XPH3(case_select),XHe(case_select),XCO,P_temp,T_temp, g0_i,R0e_i,...
                                P0_i,T_targ_i,P_targ_i,P_term_i,...
                                use_lindal,SuperSatSelf_H2S,SuperSatSelf_NH3,...
                                SuperSatSelf_PH3,SuperSatSelf_H2O,supersatNH3,...
                                supersatH2S,AutoStep_constant,fp,dz,oblateness_factor,use_dz,dP_init,dP_fine,P_fine_start,P_fine_stop,frain,select_ackerman);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Spacecraft Orientation
theta=0
Raydirection(1)=-cos(theta*(pi/180));
Raydirection(2)=sin(theta*(pi/180));

%Run Radiative Transfer model for all frequencies
for j=1:length(f)
    no_ph3=0; 
    [Tbeam_nadir(j),zenith_nadir(j),weighting_function_a_nadir(:,:,j),refractive_index(:,j)]= maintamone(Raydirection,Rayorigin,...
                                                                                   tcme,tcmp,ao,bo,co,f(j),no_ph3,select_h2h2_model,select_ammonia_model,...
                                                                                   select_water_model,include_clouds,N_ring_one,Nphi,BWHM,refractivity_source,...
                                                                                   cassini_pattern,cassini_data_path);
    clear maintamone;
    Tbeam_nadir
end


% Set Spacecraft Orientation

theta=54.5910 %along z
%theta=54.5912 %along y
%theta=60;
%Set Viewing Geometry

X_direction=-cos(theta*(pi/180));
Y_direction=0;
Z_direction=sin(theta*(pi/180));
Raydirection=[X_direction Y_direction Z_direction];

%Run Radiative Transfer Model For all Frequencies.

for j=1:length(f)
    no_ph3=0;
    [Tbeam_limb(j),zenith_limb(j),weighting_function_a_limb(:,:,j)]= maintamone(Raydirection,Rayorigin,...
                                                                                tcme,tcmp,ao,bo,co,f(j),no_ph3,select_h2h2_model,select_ammonia_model,...
                                                                                select_water_model,include_clouds,N_ring_one,Nphi,BWHM,refractivity_source,...
                                                                                cassini_pattern,cassini_data_path);
    clear maintamone;
    Tbeam_limb
end

R=100*(Tbeam_nadir-Tbeam_limb)./Tbeam_nadir;

save(output_filename);
