% Go crazy! run it a whole bunch of times!
%
clear all;
oblateness_factor=0.935;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Planet verticies                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ao=71492*100000; % along x (in centimeters)
bo=ao;       % along y
co=oblateness_factor*ao; % along z

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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load craft;

%spacecraft parameters
Rayorigin=[ao+(4500*100000) 0 0];
Sphereradius=ao;
Spherecenter=[0 0 0];
Raydirection=[-1 0 0];

% Juno bands
f=[0.6,1.2,2.4,4.8,9.6,23]; %operating frequency in GHz

%data from Lena Adams (Atreya's group)

[tcme,tcmp]=get_adams_data(oblateness_factor);

% theta=0
% Raydirection(1)=-cos(theta*(pi/180))
 %Raydirection(2)=sin(theta*(pi/180))
  for j=1:6
     no_ph3=0; 
      [Tbeam_nadir(j),zenith_nadir(j),weighting_function_a_nadir(:,:,j),weighting_function_b_nadir(:,:,j),weighting_function_c_nadir(:,:,j),weighting_function_d_nadir(:,:,j)]= maintamone(Spherecenter,Sphereradius,Raydirection,Rayorigin,tcme,tcmp,ao,bo,co,f(j),no_ph3,select_ammonia_model,select_water_model);
       Tbeam_nadir
      clear maintam;
  end

 % theta=53.9233 %along z
 theta=54.5912 %along y
 %theta=60;
  Raydirection(1)=-cos(theta*(pi/180))
  Raydirection(2)=sin(theta*(pi/180))
  no_ph3=0;
  for j=1:6
      no_ph3=0;
      [Tbeam_limb(j),zenith_limb(j),weighting_function_a_limb(:,:,j),weighting_function_b_limb(:,:,j),weighting_function_c_limb(:,:,j),weighting_function_d_limb(:,:,j)]= maintamone(Spherecenter,Sphereradius,Raydirection,Rayorigin,tcme,tcmp,ao,bo,co,f(j),no_ph3,select_ammonia_model,select_water_model);
  end
  R=100*(Tbeam_nadir-Tbeam_limb)./Tbeam_nadir;
  
save weighting_functions.mat
 %plot_replicate_hoffman_results
% jims_zenith_seventy_five
