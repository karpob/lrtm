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
co=ao; % along z

% Select a Filename
output_filename='give_me_a_cool_name.mat';

%select_ammonia_model
%1 original hoffman coding of spilker
%2 toms code joiner-steffes
%3 " berge-gulkis
%4 " mohammed-steffes
%5 " spilker
% Note, 1 and 5 won't work for current Jupiter/adams data set
% spilker correction factor C goes negative, giving negative absorption
% coefficient

select_ammonia_model=2;

%select_water_model
%1 original deboer water vapor model
%2 corrected deboer water vapor model
%3 (to be implemented goodman 1969 water vapor model
select_water_model=2;

%include cloud absorption?
%1=yes
%0=no
include_clouds=0;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load craft;

%spacecraft parameters
Rayorigin=[ao+(4500*100000) 0 0];
Sphereradius=ao;
Spherecenter=[0 0 0];
Raydirection=[-1 0 0];

%Get Beam parameters
N_ring_one=8; %Number of beams in first phi ring
BWHM=10; % Beamwidth Half-maximum

% Juno bands
f=[0.6,1.2,2.4,4.8,9.6,23]; %operating frequency in GHz

%data from Lena Adams (Atreya's group)

% For data as of Feb 6, 2007
%adams_data_directory='Adams_model';
%adams_suffix='jup_test';

%For data as of June 12, 2007 "3 times Solar H_2O, and 3 times Solar NH_3"
%adams_data_directory='adams_model_3x';
%adams_suffix='j08';

%For data as of June 12, 2007 "6 times Solar H_2O and 3 times Solar NH_3"
%adams_data_directory='adams_model_6x'
%adams_suffix='j09';
% [tcme,tcmp]=get_adams_data(oblateness_factor,adams_data_directory,adams_suffix);

 
%Use DeBoer TCM

%Water Enhancement
H2O_enhancement=3;
NH3_enhancement=3;

% Temperature forcing
TP_list={'Seiff_Jupiter','Lindal_Jupiter','Lindal_Saturn','whatever_is_in_TCM_mex'}
TP_force='Seiff_Jupiter'; %'Lindal'


%"Solar" Abundances
XH2S_rel_H2=3.1e-5;
XNH3_rel_H2=1.35e-4;
XH2O_rel_H2=1.26e-3;
XCH4_rel_H2=5.5e-4;
XPH3_rel_H2=5.14e-7;
XAr_rel_H2=3.4e-06;
XHe_rel_H2=0.194
XH2_i=1/(1+XHe_rel_H2+XH2S_rel_H2+XNH3_rel_H2+XH2O_rel_H2+ XCH4_rel_H2+ XPH3_rel_H2 + XAr_rel_H2);

%Get Mole Fractions
XH2S_i=3*XH2S_rel_H2*XH2_i
XNH3_i=NH3_enhancement*XNH3_rel_H2*XH2_i
XH2O_i=H2O_enhancement*XH2O_rel_H2*XH2_i
XCH4_i=3*XCH4_rel_H2*XH2_i
XPH3_i=3*XPH3_rel_H2*XH2_i
XAr_i=3*XAr_rel_H2*XH2_i
XHe_i=XHe_rel_H2*XH2_i
XCO=XAr_i;

%Misc DeBoer TCM inputs
P_temp=1000;
T_temp=1303.349;
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
fp=0.25;
dz=1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thermo-chemical Modeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[me,tcme,tcmp]=DeBoer_TCM(TP_list,TP_force,XH2S_i,XNH3_i,XH2O_i,XCH4_i,...
                                XPH3_i,XHe_i,XCO,P_temp,T_temp, g0_i,R0e_i,...
                                P0_i,T_targ_i,P_targ_i,P_term_i,...
                                use_lindal,SuperSatSelf_H2S,SuperSatSelf_NH3,...
                                SuperSatSelf_PH3,SuperSatSelf_H2O,supersatNH3,...
                                supersatH2S,AutoStep_constant,fp,dz,oblateness_factor);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set Spacecraft Orientation
theta=0
Raydirection(1)=-cos(theta*(pi/180));
Raydirection(2)=sin(theta*(pi/180));

%Run Radiative Transfer model for 6 frequencies
for j=1:6
    no_ph3=0; 
    [Tbeam_nadir(j),zenith_nadir(j),weighting_function_a_nadir(:,:,j)]= maintamone(Spherecenter,Sphereradius,Raydirection,Rayorigin,...
                                                                                   tcme,tcmp,ao,bo,co,f(j),no_ph3,select_ammonia_model,...
                                                                                   select_water_model,include_clouds,N_ring_one,BWHM);
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

%Run Radiative Transfer Model For 6 Frequencies.

for j=1:6
    no_ph3=0;
    [Tbeam_limb(j),zenith_limb(j),weighting_function_a_limb(:,:,j)]= maintamone(Spherecenter,Sphereradius,Raydirection,Rayorigin,...
                                                                                tcme,tcmp,ao,bo,co,f(j),no_ph3,select_ammonia_model,...
                                                                                select_water_model,include_clouds,N_ring_one,BWHM);
    clear maintamone;
    Tbeam_limb
end

R=100*(Tbeam_nadir-Tbeam_limb)./Tbeam_nadir;

save(output_filename);
