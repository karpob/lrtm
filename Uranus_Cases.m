% Go crazy! run it a whole bunch of times!
%
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Planet verticies
% Note all planet verticies are equal oblateness factor
% applied inside DeBoer_TCM function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%oblateness_factor=0.935; %Jupiter oblateness factor R_pole/R_equator
%oblateness_factor=0.902; %Saturn oblateness factor
oblateness_factor=0.977; %Uranus oblateness factor
%oblateness_factor=0.983; %Neptune oblateness factor
%ao=71492*100000; % along x (in centimeters) Jupiter
%ao=60268*100000; %along x (in centimeters) Saturn
ao=25559*100000; % along x (in centimeters) Uranus
%ao=24764*100000; % along x (in centimeters) Neptune
bo=ao;       % along y
co=ao*oblateness_factor; % along z

%select_h2h2_model
%1=joiner
%2=goodman
%3=goodman by joiner
%4=borysow
%5=borysow with orton modification
select_h2h2_model=4;

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
Rayorigin=[ao+2*ao 0 0];
Sphereradius=ao;
Spherecenter=[0 0 0];
Raydirection=[-1 0 0];

%Get Beam parameters
N_ring_one=8; %Number of beams in first phi ring
BWHM=10; % Beamwidth Half-maximum

% Juno bands
%f=[0.6,1.2];%,2.4,4.8,9.6,23]; %operating frequency in GHz
f=[10:10:225];
nfreq=length(f)
Selected_Model='Mean_Lindal'


% Temperature forcing
TP_list={'Seiff_Jupiter','Lindal_Jupiter','Lindal_Saturn','Lindal_Uranus','Lindal_Neptune','whatever_is_in_TCM_mex'}

TP_force='Lindal_Uranus';

% Select a Filename
if(include_clouds==1)
    cloudz='_Clouds';
else
    cloudz='_No_Clouds';
end
Selected_Model='Uranus';
output_filename=strcat(Selected_Model,cloudz,'.mat');


%Abundances
XH2S_rel_H2=1.18e-4;
XNH3_rel_H2=1.53e-4;
XH2O_rel_H2=1.213e-3;
XCH4_rel_H2=2.0e-2;
XPH3_rel_H2=4.297e-7;
XHe_rel_H2=0.15;
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
P_temp=500;
T_temp=573;
g0_i=870; %2417;
R0e_i=ao;
P0_i=1;
T_targ_i=100.9;
P_targ_i=2.309;
P_term_i=0.001;
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
XCO=0;
 table_output=[XH2;XHe;(1e6)*XH2S;(1e6)*XNH3;(1e6)*XH2O;(1e6)*XCH4;(1e6)*XPH3];
    to_dlm=transpose(table_output);
    dlmwrite('mole_fractions_uranus.dat',to_dlm,'delimiter','&','precision','%.4f')

%for i=1:length(XNH3)
%    if(strcmp(Selected_Model,Model_Names(i)))
%        case_select=i;
%    end
%end
case_select=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thermo-chemical Modeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[me,tcme,tcmp]=DeBoer_TCM(TP_list,TP_force,XH2S(case_select),XNH3(case_select),XH2O(case_select),XCH4(case_select),...
                                XPH3(case_select),XHe(case_select),XCO,P_temp,T_temp, g0_i,R0e_i,...
                                P0_i,T_targ_i,P_targ_i,P_term_i,...
                                use_lindal,SuperSatSelf_H2S,SuperSatSelf_NH3,...
                                SuperSatSelf_PH3,SuperSatSelf_H2O,supersatNH3,...
                                supersatH2S,AutoStep_constant,fp,dz,oblateness_factor);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Spacecraft Orientation
theta=0
Raydirection(1)=-cos(theta*(pi/180));
Raydirection(2)=sin(theta*(pi/180));

%Run Radiative Transfer model for all frequencies
for j=1:length(f)
    no_ph3=0; 
    [Tbeam_nadir(j),zenith_nadir(j),weighting_function_a_nadir(:,:,j)]= maintamone(Spherecenter,Sphereradius,Raydirection,Rayorigin,...
                                                                                   tcme,tcmp,ao,bo,co,f(j),no_ph3,select_h2h2_model,select_ammonia_model,...
                                                                                   select_water_model,include_clouds,N_ring_one,BWHM);
    clear maintamone;
    Tbeam_nadir
end


% Set Spacecraft Orientation

theta=17 %along z
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
    [Tbeam_limb(j),zenith_limb(j),weighting_function_a_limb(:,:,j)]= maintamone(Spherecenter,Sphereradius,Raydirection,Rayorigin,...
                                                                                tcme,tcmp,ao,bo,co,f(j),no_ph3,select_h2h2_model,select_ammonia_model,...
                                                                                select_water_model,include_clouds,N_ring_one,BWHM);
    clear maintamone;
    Tbeam_limb
end

R=100*(Tbeam_nadir-Tbeam_limb)./Tbeam_nadir;

save(output_filename);
