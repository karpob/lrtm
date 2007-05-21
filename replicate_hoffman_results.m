% Go crazy! run it a whole bunch of times!
%
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Planet verticies at 1 bar pressure                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%He=54.5;
%Hp=41.0;
ao=6.0268e9; % along x
bo=ao;       % along y
co=5.4364e9; % along z


%ao=ao+He*(1/0.5);
%bo=ao+He*(1/0.5);
%co=bo+Hp*(1/0.5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
He_range=[0.01:0.01:0.1]
[n,mm]=size(He_range)

zstep=2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                   %
% Input variables from Hoffman's Thesis p.97        % 
%  table 4.3:                                       %
%                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Model_names={'A','A2','A3','A4','B','B2','B3','C','C2','C3','C4','D','E','F','G'};
Tbeam_thesis=[146.68,146.54,146.28,146.68,145.89,146.13,146.27,146.43,146.13,146.30,146.56,146.18,145.90,145.40,144.75];

Tbeam_thesis_seventy_five=[139.58,139.24,138.62,139.58,138.45,139.06,139.42,138.91,138.89,138.60,139.23,138.90,138.63,138.04,137.51];

dz=[zstep,zstep,zstep,zstep,zstep,zstep,zstep,zstep,zstep,zstep,zstep,zstep,zstep,zstep,zstep]; % Set dz to zero for autostep
%dz=zstep:0.01:5;

XH2S_i=[410e-6,410e-6,410e-6,410e-6,30.8e-6,30.8e-6,30.8e-6,308e-6,154e-6,308e-6,308e-6,30.8e-6,0.0,30.8e-6,30.8e-6];

XNH3_i=[520e-6,520e-6,520e-6,520e-6,468e-6,468e-6,468e-6,412.5e-6,412.5e-6,412.5e-6,412.5e-6,281e-6,468e-6,468e-6,468e-6];

XH2O_i=(1e-6)*[1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0];

XCH4_i=[4.2e-3,4.2e-3,4.2e-3,4.2e-3,4.2e-3,4.2e-3,4.2e-3,4.2e-3,4.2e-3,4.2e-3,4.2e-3,4.2e-3,4.2e-3,4.2e-3,4.2e-3];

XPH3_i=[0,6.2e-6,12.4e-6,3.1e-6,12.4e-6,6.2e-6,3.1e-6,9.3e-6,9.3e-6,12.4e-6,6.2e-6,9.3e-6,6.2e-6,12.4e-6,12.4e-6];

XCO=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];

P_temp=[70.0,70.0,70.0,70.0,70.0,70.0,70.0,70.0,70.0,70.0,70.0,70.0,70.0,70.0,70.0];

T_temp=[489.5,489.5,489.5,489.5,489.5,489.5,489.5,489.5,489.5,489.5,489.5,489.5,489.5,489.5,489.5];

g0_i=[900.0,900.0,900.0,900.0,900.0,900.0,900.0,900.0,900.0,900.0,900.0,900.0,900.0,900.0,900.0];

g0_p=[900.0,900.0,900.0,900.0,900.0,900.0,900.0,900.0,900.0,900.0,900.0,900.0,900.0,900.0,900.0];

%g0_p=[1200.0,1200.0,1200.0,1200.0,1200.0,1200.0,1200.0,1200.0,1200.0,1200.0,1200.0,1200.0,1200.0,1200.0,1200.0];

R0e_i=[ao,ao,ao,ao,ao,ao,ao,ao,ao,ao,ao,ao,ao,ao,ao];

R0p_i=[co,co,co,co,co,co,co,co,co,co,co,co,co,co,co];

P0_i=[1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0];

T_targ_i=[146.2,146.2,146.2,146.2,146.2,146.2,146.2,146.2,146.2,146.2,146.2,146.2,146.2,146.2,146.2];

P_targ_i=[1.29848,1.29848,1.29848,1.29848,1.29848,1.29848,1.29848,1.29848,1.29848,1.29848,1.29848,1.29848,1.29848,1.29848,1.29848];

P_term_i=(0.5/0.0631)*[0.0631,0.0631,0.0631,0.0631,0.0631,0.0631,0.0631,0.0631,0.0631,0.0631,0.0631,0.0631,0.0631,0.0631,0.0631];

use_lindal={'Y','Y','Y','Y','Y','Y','Y','Y','Y','Y','Y','Y','Y','Y','Y'};

n_lindal_pts=[27,27,27,27,27,27,27,27,27,27,27,27,27,27,27];

SuperSatSelf_H2S=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];

SuperSatSelf_NH3=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];

SuperSatSelf_PH3=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];

SuperSatSelf_H2O=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];

supersatNH3=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.1,0.2];

supersatH2S=[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0];

AutoStep_constant_range=8*[1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0];
%AutoStep_constant=8;
He=0.03;
XHe_i=[He,He,He,He,He,He,He,He,He,He,He,He,He,He,He];
%AutoStep_range=5.6:0.01:7;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Run TCM, and arrange variables for maintam      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for kk=1:1%length(AutoStep_range)
AutoStep_constant=8;
for k=1:1
    cd TCM_mex
    [P,T,XH2,XHe,XH2S,XNH3,XH20,XCH4,XPH3,clouds,DNH4SH,DH2S,DNH3,DH2O,DCH4,DPH3,DSOL,g,mu,ref_w_o,ref_w,z]=TCM(dz(k),XHe_i(k),XH2S_i(k),XNH3_i(k),XH2O_i(k),XCH4_i(k),XPH3_i(k),XCO(k),P_temp(k),T_temp(k),g0_i(k),R0e_i(k),P0_i(k),T_targ_i(k),P_targ_i(k),P_term_i(k),1,27,SuperSatSelf_H2S(k),SuperSatSelf_NH3(k),SuperSatSelf_PH3(k),SuperSatSelf_H2O(k),supersatNH3(k),supersatH2S(k),AutoStep_constant);
    clear TCM;
    
   
    tcme(:,:,k)=load('tcm.out');
%    tcmp(1:me(k),1:22,k)=[tcme(1:me(k),1:2,k),(co/ao).*tcme(1:me(k),3,k),tcme(1:me(k),4:22,k)];
%    tcme(:,:,k)=tcmel(1:slicefun,1:22,k);
    tcmp(:,1:22,k)=[tcme(:,1:2,k),(co/ao).*tcme(:,3,k),tcme(:,4:22,k)];
    cd ..
  %  tcme(1:me(k),1:22,k)=[P(1:me(k)),T(1:me(k)),z(1:me(k)),XH2(1:me(k)),XHe(1:me(k)),XH2S(1:me(k)),XNH3(1:me(k)),XH20(1:me(k)),XCH4(1:me(k)),XPH3(1:me(k)),clouds(1:me(k)),DNH4SH(1:me(k)),DH2S(1:me(k)),DNH3(1:me(k)),DH2O(1:me(k)),DCH4(1:me(k)),DPH3(1:me(k)),DSOL(1:me(k)),g(1:me(k)),mu(1:me(k)),ref_w_o(1:me(k)),ref_w(1:me(k))];
end
%  for k=1:1
%      cd TCM_mex    
%      [P,T,XH2,XHe,XH2S,XNH3,XH20,XCH4,XPH3,clouds,DNH4SH,DH2S,DNH3,DH2O,DCH4,DPH3,DSOL,g,mu,ref_w_o,ref_w,z]=TCM(dz(k),XHe_i(k),XH2S_i(k),XNH3_i(k),XH2O_i(k),XCH4_i(k),XPH3_i(k),XCO(k),P_temp(k),T_temp(k),g0_p(k),R0p_i(k),P0_i(k),T_targ_i(k),P_targ_i(k),P_term_i(k),1,27,SuperSatSelf_H2S(k),SuperSatSelf_NH3(k),SuperSatSelf_PH3(k),SuperSatSelf_H2O(k),supersatNH3(k),supersatH2S(k),AutoStep_constant);
%  %    cd ..
%      clear TCM;
%      
%      [mp(k),n]=size(P)
%      tcmp(:,:,k)=load('tcm.out');
%%  %    tcmp(:,:,k)=tcmpl(1:slicefun,1:22,k);
%      cd ..
%%    %tcmp(1:mp(k),1:22,k)=[P(1:mp(k)),T(1:mp(k)),z(1:mp(k)),XH2(1:mp(k)),XHe(1:mp(k)),XH2S(1:mp(k)),XNH3(1:mp(k)),XH20(1:mp(k)),XCH4(1:mp(k)),XPH3(1:mp(k)),clouds(1:mp(k)),DNH4SH(1:mp(k)),DH2S(1:mp(k)),DNH3(1:mp(k)),DH2O(1:mp(k)),DCH4(1:mp(k)),DPH3(1:mp(k)),DSOL(1:mp(k)),g(1:mp(k)),mu(1:mp(k)),ref_w_o(1:mp(k)),ref_w(1:mp(k))];
%  end
 % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load craft;
Raydirection
Rayorigin

f=13.78; %operating frequency in GHz

%select_ammonia_model
%1 original hoffman coding of spilker
%2 toms code joiner-steffes
%3 " berge-gulkis
%4 " mohammed-steffes
%5 " spilker
% Note, 1 and 5 won't work for current Jupiter/adams data set
% spilker correction factor C goes negative, giving negative absorption
% coefficient

select_ammonia_model=1;

%select_water_model
%1 original deboer water vapor model
%2 corrected deboer water vapor model
%3 (to be implemented goodman 1969 water vapor model
select_water_model=1;
 
%Rayorigin=[4.2e10 0 0]
%Sphereradius=ao;
% theta=7.889
theta=0;
 Raydirection(1)=-cos(theta*(pi/180))
 Raydirection(2)=sin(theta*(pi/180))

% ao=6.0268e9; % along x
%bo=ao;       % along y
%co=5.4364e9; % along z
%tcmp(1:me(1),3,1)
for j=1:15
    if(j==1)
        no_ph3=0;
    else
        no_ph3=1;
    end
%    me(k)
    [Tbeam(j),jims_zenith_seventy_five(j)]= maintamone(Spherecenter,Sphereradius,Raydirection,Rayorigin,tcme(:,:,j),tcmp(:,:,j),ao,bo,co,f,no_ph3,select_ammonia_model,select_water_model);
    Tbeam(j)
    Model_names(j)
    residual(j)=Tbeam(j)-Tbeam_thesis(j)
    clear maintam;
end
end 
% theta=7.5
% Raydirection(1)=-cos(theta*(pi/180))
% Raydirection(2)=sin(theta*(pi/180))
% 
% for j=1:15
%     [Tbeam_seventy_five(j),jims_zenith_seventy_five(j)]= maintam(Spherecenter,Sphereradius,Raydirection,Rayorigin,tcme(1:me(j),:,j),tcmp(1:mp(j),:,j),ao,bo,co,f);
% end
% %plot_replicate_hoffman_results
% jims_zenith_seventy_five