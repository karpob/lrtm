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

%  theta=0
%  Raydirection(1)=-cos(theta*(pi/180));
%  Raydirection(2)=sin(theta*(pi/180));
%     for j=1:6
%        no_ph3=0; 
%         [Tbeam_nadir(j),zenith_nadir(j),weighting_function_a_nadir(:,:,j)]= maintamone(Spherecenter,Sphereradius,Raydirection,Rayorigin,tcme,tcmp,ao,bo,co,f(j),no_ph3,select_ammonia_model,select_water_model);
%          Tbeam_nadir
%         clear maintamone;
%     end

  theta=54.59115 %along z
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
      [Tbeam_limb(j),zenith_limb(j),weighting_function_a_limb(:,:,j),intercepts_boresight(:,:,j),intercepts_b(:,:,:,j),intercepts_c(:,:,:,j),intercepts_d(:,:,:,j)]= maintamone(Spherecenter,Sphereradius,Raydirection,Rayorigin,tcme,tcmp,ao,bo,co,f(j),no_ph3,select_ammonia_model,select_water_model);
      clear maintamone;
      Tbeam_limb
  end
 % R=100*(Tbeam_nadir-Tbeam_limb)./Tbeam_nadir;

save weighting_functions.dat
 ray_antenna_to_planet=[Rayorigin(1) 0 0; intercepts_boresight(1,:,1)]
 plot3(ray_antenna_to_planet(:,1),ray_antenna_to_planet(:,2),ray_antenna_to_planet(:,3),'color','g');
 hold on;
 
 x_axis=[-Rayorigin(1)*1 0 0;Rayorigin(1)*1 0 0];
 plot3(x_axis(:,1),x_axis(:,2),x_axis(:,3),'color','k');
 
 y_axis=[0 -Rayorigin(1)*1  0;0 Rayorigin(1)*1  0]
 plot3(y_axis(:,1),y_axis(:,2),y_axis(:,3),'color','k');
 
 z_axis=[0 0 -Rayorigin(1);0 0 Rayorigin(1)*1 ];
 plot3(z_axis(:,1),z_axis(:,2),z_axis(:,3),'color','k');
 plot3(intercepts_boresight(:,1),intercepts_boresight(:,2),intercepts_boresight(:,3))
 
 [m,n,o,p]=size(intercepts_b);
 for i=1:m
    plot3(intercepts_b(i,:,1,1),intercepts_b(i,:,2,1),intercepts_b(i,:,3,1),'color','r')
 end
 
[m,n,o,p]=size(intercepts_c);
for i=1:m
    plot3(intercepts_c(i,:,1,1),intercepts_c(i,:,2,1),intercepts_c(i,:,3,1),'color','m')
end
[m,n,o,p]=size(intercepts_d);
for i=1:m
    plot3(intercepts_d(i,:,1,1),intercepts_d(i,:,2,1),intercepts_d(i,:,3,1),'color','c')
end

