% Go crazy! run it a whole bunch of times!
%
clear all;
oblateness_factor=0.935;
output_filename='orton_test_jupiter.mat'
[tcme,tcmp,Re,gravity]= get_orton_data(oblateness_factor,'Jupiter');
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
%6=borysow with orton et al modification with a second modification below 3 GHz
select_h2h2_model=6;

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

% refractivity_source
% Select the author you believe is right with regards to values for refractivity (used for raypath calculations)
%
%refractivity_source=0; % No bending due to refraction n=1.0D0
 refractivity_source=1; % Original DeBoer/Hoffman H2/He refractivity 
% refractivity_source=2; % Karpowicz H2/He refractivity using original Essen data
% refractivity_source=3; % Karpowicz H2, He, CH4 etc.. using Essen, and other sources
%refractivity_source=4; % Karpowicz w/Clouds H2, He, CH4 etc.. using Essen, and other sources
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load craft;

%spacecraft parameters
Rayorigin=[ao+(4500*100000) 0 0];
%Rayorigin=[ao+1000000 0 0]
Raydirection=[-1 0 0];

%Get Beam parameters
N_ring_one=4; %Number of beams in first phi ring
Nphi=3; %Number of rings in phi
BWHM=10; % Beamwidth Half-maximum

% Juno bands
%f=[0.6,1.2];%,2.4,4.8,9.6,23]; %operating frequency in GHz
%f=[0.5:0.1:10,10:1:25];

f=[0.999308193333;1.498962290;4.99654096667;14.98962290;29.979245800;42.827494000;99.930819333;299.792459;999.308193333;2997.92458];
disp('poop2')
theta=0
Raydirection(1)=-cos(theta*(pi/180));
Raydirection(2)=sin(theta*(pi/180));

%Run Radiative Transfer model for all frequencies
for j=1:length(f)
    no_ph3=0; 
    [Tbeam_nadir(j),zenith_nadir(j),weighting_function_a_nadir(:,:,j),refractive_index(:,j),Tatma(j)]= maintamone(Raydirection,Rayorigin,...
                                                                                   tcme,tcmp,ao,bo,co,f(j),no_ph3,select_h2h2_model,select_ammonia_model,...
                                                                                   select_water_model,include_clouds,N_ring_one,Nphi,BWHM,refractivity_source);
    clear maintamone;
    Tatma
end


% Set Spacecraft Orientation

theta=61.74;%30;%54.5910 %along z
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
    [Tbeam_limb(j),zenith_limb(j),weighting_function_a_limb(:,:,j),refractive_index(:,j),Tatm_a_limb(j)]= maintamone(Raydirection,Rayorigin,...
                                                                                tcme,tcmp,ao,bo,co,f(j),no_ph3,select_h2h2_model,select_ammonia_model,...
                                                                       select_water_model,include_clouds,N_ring_one,Nphi,BWHM,refractivity_source);
    clear maintamone;
    disp('poop3') 
    Tatm_a_limb
    cos((pi/180)*zenith_limb(j))
end

R=100*(Tbeam_nadir-Tbeam_limb)./Tbeam_nadir;

save(output_filename);
