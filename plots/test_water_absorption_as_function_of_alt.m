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
P=tcme(:,1);
T=tcme(:,2);
XH2=tcme(:,4);
XHe=tcme(:,5);
P_He=P.*XHe;
P_H2=P.*XH2;
XH2O=tcme(:,8);
P_H2O=P.*XH2O;

cd h2o
for i=1:length(P)
    alpha_watervapor_deboer(i)=fwatervaporalpha_deboer_corrected(f(1),T(i),P_H2(i),P_He(i),P_H2O(i));
    alpha_watervapor_goodman(i)=fwatervaporalpha_goodman(f(1),P(i),T(i),XH2O(i),XH2(i),XHe(i));
    alpha_watervapor_deboer_2x(i)=fwatervaporalpha_deboer_corrected(f(1),T(i),P_H2(i),P_He(i),2*P_H2O(i));
    alpha_watervapor_goodman_2x(i)=fwatervaporalpha_goodman(f(1),P(i),T(i),2*XH2O(i),XH2(i),XHe(i));
end
cd ..

figure(1)
hold on
plot(alpha_watervapor_deboer,P,'k')
plot(alpha_watervapor_deboer_2x,P,'k--')
plot(alpha_watervapor_goodman,P,'r')
plot(alpha_watervapor_goodman_2x,P,'r--')

legend('Deboer Water Vapor','Deboer Water Vapor 2x','Goodman Water Vapor','Goodman Water Vapor 2x');

% % theta=0
% % Raydirection(1)=-cos(theta*(pi/180))
%  %Raydirection(2)=sin(theta*(pi/180))
%   for j=1:6
%      no_ph3=0; 
%       [Tbeam_nadir(j),zenith_angle_nadir(j)]= maintamone(Spherecenter,Sphereradius,Raydirection,Rayorigin,tcme,tcmp,ao,bo,co,f(j),no_ph3,select_ammonia_model,select_water_model);
%        Tbeam_nadir
%       clear maintam;
%   end
% 
%  % theta=53.9233 %along z
%  theta=54.5912 %along y
%  %theta=60;
%   Raydirection(1)=-cos(theta*(pi/180))
%   Raydirection(2)=sin(theta*(pi/180))
%   no_ph3=0;
%   for j=1:6
%       no_ph3=0;
%       [Tbeam_limb(j),zenith_limb(j)]= maintamone(Spherecenter,Sphereradius,Raydirection,Rayorigin,tcme,tcmp,ao,bo,co,f(j),no_ph3,select_ammonia_model,select_water_model);
%   end
%   R=100*(Tbeam_nadir-Tbeam_limb)./Tbeam_nadir;
%   
% save Jupiter_joiner_steffes_wv_Goodman_2times_h2o.mat
 %plot_replicate_hoffman_results
% jims_zenith_seventy_five