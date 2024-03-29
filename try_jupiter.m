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
%select_h2h2_model
%1=joiner
%2=goodman
%3=goodman by joiner
%4=borysow
%5=borysow with orton modification
select_h2h2_model=5;
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
include_clouds=1;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load craft;

%spacecraft parameters
Rayorigin=[ao+(4500*100000) 0 0];
Sphereradius=ao;
Spherecenter=[0 0 0];
Raydirection=[-1 0 0];

%Get Beam parameters
N_ring_one=4; % Number of samples in first ring
BWHM=10; % Beamwidth Half-maximum


% Juno bands
f=[0.6,1.2,2.4,4.8,9.6,23]; %operating frequency in GHz

%data from Lena Adams (Atreya's group)

% For data as of Feb 6, 2007
%adams_data_directory='Adams_model';
%adams_suffix='jup_test';

%For data as of June 12, 2007 "3 times Solar H_2O, and 3 times Solar NH_3"
adams_data_directory='adams_model_3x';
adams_suffix='j08';

%For data as of June 12, 2007 "6 times Solar H_2O and 3 times Solar NH_3"
%adams_data_directory='adams_model_6x'
%adams_suffix='j09';
 [tcme,tcmp]=get_adams_data(oblateness_factor,adams_data_directory,adams_suffix);

  theta=0
  Raydirection(1)=-cos(theta*(pi/180));
  Raydirection(2)=sin(theta*(pi/180));
     for j=1:6
        no_ph3=0; 
         [Tbeam_nadir(j),zenith_nadir(j),weighting_function_a_nadir(:,:,j)]= maintamone(Spherecenter,Sphereradius,Raydirection,Rayorigin,tcme,tcmp,ao,bo,co,f(j),no_ph3,select_h2h2_model,select_ammonia_model,select_water_model,include_clouds,N_ring_one,BWHM);
          Tbeam_nadir
         clear maintamone;
     end

  theta=54.5910 %along z
 %theta=54.5912 %along y
 %theta=60;
 %Set Viewing Geometry
 
 X_direction=-cos(theta*(pi/180));
 Y_direction=0;
 Z_direction=sin(theta*(pi/180));
 
  Raydirection=[X_direction Y_direction Z_direction];
  no_ph3=0;
  for j=1:6
      no_ph3=0;
      [Tbeam_limb(j),zenith_limb(j),weighting_function_a_limb(:,:,j)]= maintamone(Spherecenter,Sphereradius,Raydirection,Rayorigin,tcme,tcmp,ao,bo,co,f(j),no_ph3,select_h2h2_model,select_ammonia_model,select_water_model,include_clouds,N_ring_one,BWHM);
      clear maintamone;
      Tbeam_limb
  end
  R=100*(Tbeam_nadir-Tbeam_limb)./Tbeam_nadir;

save(output_filename);


