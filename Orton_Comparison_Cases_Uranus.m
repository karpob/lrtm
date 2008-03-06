% Go crazy! run it a whole bunch of times!
%
clear all;
%oblateness_factor=0.935; % Jupiter
%oblateness_factor=0.902; % Saturn
%oblateness_factor=1;
oblateness_factor=0.977; %Uranus
output_filename='orton_test_Uranus_borysowh2h2).mat'
[tcme,tcmp,Re,gravity]= get_orton_data(oblateness_factor,'Uranus');
%oblateness_factor=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Planet verticies                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ao_jupiter=71492*100000; % along x (in centimeters)
ao_saturn=60268*100000;
ao_uranus=25559*100000;

ao=ao_uranus;
bo=ao;       % along y
co=ao*oblateness_factor; % along z

%select_h2h2_model
%1=joiner
%2=goodman
%3=goodman by joiner
%4=borysow
%5=borysow with orton et al, 2007 modification
%6=orton et al based upon interpolation
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

% refractivity_source
% Select the author you believe is right with regards to values for refractivity (used for raypath calculations)
%
%refractivity_source=0; % No bending due to refraction n=1.0D0
% refractivity_source=1; % Original DeBoer/Hoffman H2/He refractivity 
 refractivity_source=2; % Karpowicz H2/He refractivity using original Essen data
% refractivity_source=3; % Karpowicz H2, He, CH4 etc.. using Essen, and other sources
%refractivity_source=4; % Karpowicz w/Clouds H2, He, CH4 etc.. using Essen, and other sources
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load craft;

%spacecraft parameters
Rayorigin=[ao+(45000*100000) 0 0];
%Rayorigin=[ao+1000000 0 0]
Raydirection=[-1 0 0];

%Get Beam parameters
N_ring_one=4; %Number of beams in first phi ring
Nphi=3; %Number of rings in phi
BWHM=10; % Beamwidth Half-maximum

% Juno bands
%f=[0.6,1.2];%,2.4,4.8,9.6,23]; %operating frequency in GHz
%f=[0.5:0.1:10,10:1:25];

f=[0.999308193333;1.498962290;4.99654096667;14.98962290;29.979245800;42.827494000;99.930819333;299.792459;999.308193333;2997.92458]
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

theta=19.625;%30;%54.5910 %along z
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

orton=load('orton_uranus_Tb.dat');
delta_nadir=100*(Tatma-transpose(orton(:,1)))./transpose(orton(:,1));
delta_mu=100*(Tatm_a_limb-transpose(orton(:,2)))./transpose(orton(:,2));
M=[Tatma;Tatm_a_limb;delta_nadir;delta_mu];
dlmwrite('Uranus.dat',transpose(M),'delimiter','\t','precision',6)

save(output_filename);
