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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Change Number of radii from 2-8
number_of_saturn_radii=2 
include_clouds=0;
theta=17

% Name the file
output_filename='Saturn_Rsat_2_No_Clouds_theta_17.mat'
%
% Cases for Sept 2007 Modeling Study
%
Model_names={'Mean','En_Ammonia','Dep_Ammonia','En_Water','Dep_Water',...
              'En_Phosphine','Dep_Phosphine','Mean_p_warming','Dep_Ammonia_p_Warming','Absent_Ammonia'};

%Fractions relative to H2

qHe=0.135;
qH2S=2.5874e-04;

qNH3=[3.5e-4, 6e-4, 1e-4, 3.5e-4, 3.5e-4...
      3.5e-4, 3.5e-4, 3.5e-4, 1e-4, 0.0 ];
  
qH2O=[3.6e-4, 3.6e-4, 3.6e-4, 5.5e-4 1.7e-4...
      3.6e-4, 3.6e-4, 3.6e-4, 3.6e-5 3.6e-4 ];

qCH4=5.1e-3;

qPH3=[ 6.4e-6, 6.4e-6, 6.4e-6, 6.4e-6, 6.4e-6,...
       7.4e-6, 3.9e-6, 6.4e-6, 6.4e-6, 6.4e-6 ];

%Calculate Actual Mole Fraction of each species
for i=1:length(qNH3)
    Model_names(i);
    xH2(i)=1/(1+qHe+qH2S+qNH3(i)+qCH4+qPH3(i));
    xHe(i)=qHe*xH2(i);
    xH2S(i)=qH2S*xH2(i);
    xNH3(i)=qNH3(i)*xH2(i);
    xH2O(i)=qH2O(i)*xH2(i);
    xCH4(i)=qCH4*xH2(i);
    xPH3(i)=qPH3(i)*xH2(i);    
end

% Write the table of values to LaTeX (sort of) with minor species in ppm
table_output=[xH2;xHe;(1e6)*xH2S;(1e6)*xNH3;(1e6)*xH2O;(1e6)*xCH4;(1e6)*xPH3];
to_dlm=transpose(table_output);
dlmwrite('mole_fractions_saturn.dat',to_dlm,'delimiter','&','precision','%.4f')


   
%Misc DeBoer TCM params
fp=-1;
xCO=0;
P_temp=70;
T_temp=489.5;
g0=900.0;
R0=ao;
P0=1.0
T_targ=146.2;
P_targ=1.29848
P_term=0.0631;
use_lindal=1;
n_lindal_pts=27;
SuperSatSelf_H2S=0.0;
SuperSatSelf_NH3=0.0;
SuperSatSelf_PH3=0.0;
SuperSatSelf_H2O=0.0;
supersatNH3=0.0;
supersatH2S=0.0;
AutoStep_constant=8;
dz=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Run TCM, and arrange variables for maintam      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AutoStep_constant=8;

for k=1:length(xHe)
    % If cases includes a 10% warming, do it, otherwise copy lindal 1985
    % profile
    if(strcmp(Model_names(k),'Mean_p_warming')|strcmp(Model_names(k),'Dep_Ammonia_p_Warming'))
        TP_loaded=load('TCM_mex/TP.SAT'); % load lindal profile
        TP_modified=[TP_loaded(:,1),TP_loaded(:,2)*1.1];
        dlmwrite('TCM_mex/TP.high',TP_modified,'delimiter','\t','precision','%.4f');
        system('cp TCM_mex/TP.high TCM_mex/TP.TCM')
    else
        system('cp TCM_mex/TP.SAT TCM_mex/TP.TCM')
    end
    
    cd TCM_mex
    [P,T,XH2,XHe,XH2S,XNH3,XH20,XCH4,XPH3,...
     clouds,DNH4SH_i,DH2S_i,DNH3_i,DH2O_i,DCH4_i,DPH3_i,DSOL_i,...
     g,mu,ref_w_o,ref_w,z]=TCM(dz,xHe(k),xH2S(k),xNH3(k),xH2O(k),...
                               xCH4(k),xPH3(k),xCO,...
                               P_temp,T_temp,g0,...
                               R0,P0,T_targ,...
                               P_targ,P_term,use_lindal,n_lindal_pts,...
                               SuperSatSelf_H2S,SuperSatSelf_NH3,...
                               SuperSatSelf_PH3,SuperSatSelf_H2O,...
                               supersatNH3,supersatH2S,...
                               AutoStep_constant,fp);
    clear TCM;
    [me(k),n]=size(P);
    cd ..
    
    [DNH4SH,DH2S,DNH3,DH2O,DCH4,DPH3,DSOL]=filter_clouds(clouds,DNH4SH_i,DH2S_i,DNH3_i,DH2O_i,DCH4_i,DPH3_i,DSOL_i);
   
    tcme(1:me(k),1:22,k)=[P(1:me(k)),T(1:me(k)),z(1:me(k)),XH2(1:me(k)),...
                          XHe(1:me(k)),XH2S(1:me(k)),XNH3(1:me(k)),...
                          XH20(1:me(k)),XCH4(1:me(k)),XPH3(1:me(k)),...
                          clouds(1:me(k)),DNH4SH(1:me(k)),DH2S(1:me(k)),...
                          DNH3(1:me(k)),DH2O(1:me(k)),DCH4(1:me(k)),DPH3(1:me(k)),...
                          DSOL(1:me(k)),g(1:me(k)),mu(1:me(k)),ref_w_o(1:me(k)),ref_w(1:me(k))];
                      
    tcmp(1:me(k),1:22,k)=[tcme(1:me(k),1:2,k),(co/ao).*tcme(1:me(k),3,k),tcme(1:me(k),4:22,k)];
end

%Filter values for clouds, make sure that the "clouds" variable is
%consistent with cloud bulk density


%
%


%Beam pattern parameters
N_ring_one=8;
BWHM=0.37;

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

select_ammonia_model=2;

%select_water_model
%1 original deboer water vapor model
%2 corrected deboer water vapor model
%3 (to be implemented goodman 1969 water vapor model
select_water_model=2;



%Spacecraft Parameters
Rayorigin=[(number_of_saturn_radii+1)*ao 0 0]
Sphereradius=ao;
Spherecenter=[0 0 0];
Raydirection=[-1 0 0];

no_ph3=1; %use PH3 decay
for j=1:length(xHe)
    [Tbeam_nadir(j),zenith_nadir(j),weighting_function_a_nadir(1:me(j)+1,j)]= maintamone(Spherecenter,Sphereradius,Raydirection,Rayorigin,...
                           tcme(1:me(j),1:22,j),tcmp(1:me(j),1:22,j),ao,bo,co,...
                           f,no_ph3,select_ammonia_model,select_water_model,...
                           include_clouds,N_ring_one,BWHM);
end

X_direction=-cos(theta*(pi/180));
Y_direction=0;
Z_direction=sin(theta*(pi/180));
Raydirection=[X_direction Y_direction Z_direction];

for j=1:length(xHe)       
    [Tbeam_limb(j),zenith_limb(j),weighting_function_a_limb(1:me(j)+1,j)]= maintamone(Spherecenter,Sphereradius,Raydirection,Rayorigin,...
                                tcme(1:me(j),1:22,j),tcmp(1:me(j),1:22,j),ao,bo,co,...
                                f,no_ph3,select_ammonia_model,select_water_model,...
                                include_clouds,N_ring_one,BWHM);
end

R=100*(Tbeam_nadir-Tbeam_limb)./Tbeam_nadir;
save(output_filename);
