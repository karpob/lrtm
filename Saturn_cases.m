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
oblateness_factor=(co/ao);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_filename='Saturn_Rsat_2.mat'

%Change Number of radii from 2-8
number_of_saturn_radii=2 

for cases=1:4
    if(cases<=4)
        theta=17.15;
        BWHM=0.37;
    %elseif(cases<=8)
    %    theta=16.52;
    %    BWHM=1.0;
    %elseif(cases<=12)
    %    theta=17.15;
    %    BWHM=0.37;
    %elseif(cases<=16)
    %    theta=16.52;
    %    BWHM=1.0;
    end
    
    if((cases==1)|(cases==3))%((cases==1)|(cases==2)|(cases==5)|(cases==6)|(cases==9)|(cases==10)|(cases==13)|(cases==14))
        select_ammonia_model=2;
    else
        select_ammonia_model=3;
    end
    %if%((cases==1)|(cases==3)|(cases==9)|(cases==11)|(cases==13)|(cases==15))
        select_water_model=2;
    %else
    %    select_water_model=3;
    %end
    if(cases<3)
        include_clouds=0;
    else
        include_clouds=1;
    end

%Beam pattern parameters
    N_ring_one=8;% First ring has 8 rays
    Nphi=3; % 3 rings


    f=13.78; %operating frequency in GHz

%select_h2h2_model
%1=joiner
%2=goodman
%3=goodman by joiner
%4=borysow
%5=borysow with orton modification
select_h2h2_model=1;

%select_ammonia_model
%1 original hoffman coding of spilker
%2 toms code joiner-steffes
%3 " berge-gulkis
%4 " mohammed-steffes
%5 " spilker
% Note, 1 and 5 won't work for current Jupiter/adams data set
% spilker correction factor C goes negative, giving negative absorption
% coefficient

%select_ammonia_model=2;

%select_water_model
%1 original deboer water vapor model
%2 corrected deboer water vapor model
%3 goodman 1969 water vapor model
%select_water_model=2;

% refractivity_source
% Select the author you believe is right with regards to values for refractivity (used for raypath calculations)
%
% refractivity_source=0; % No bending due to refraction n=1.0D0
 refractivity_source=1; % Original DeBoer/Hoffman H2/He refractivity 
% refractivity_source=2; % Karpowicz H2/He refractivity using original Essen data
% refractivity_source=3; % Karpowicz H2, He, CH4 etc.. using Essen, and other sources
%refractivity_source=4; % Karpowicz w/Clouds H2, He, CH4 etc.. using Essen, and other sources
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Cases for Sept 2007 Modeling Study
%
    Model_names={'En_Ammonia','Mean','40%1.4bar','40%2bar','Mean_p_warming','10per','Depleted'};

%Fractions relative to H2

    qHe=0.135;
    qH2S=2.5874e-04;

    qNH3=[6e-4,4.3e-4, 4.3e-4,4.3e-4, 4.3e-4,4.3e-4,2.6e-4];
  
    qH2O=[3.6e-3, 3.6e-3, 3.6e-3, 3.6e-3,3.6e-3,3.6e-3,3.6e-3 ];

    qCH4=5.1e-3;

    qPH3=[6.4e-6, 6.4e-6, 6.4e-6, 6.4e-6, 6.4e-6,6.4e-6,6.4e-6 ];

%Calculate Actual Mole Fraction of each species
    for i=1:length(qNH3)
        Model_names(i);
        xH2(i)=1/(1+qHe+qH2S+qNH3(i)+qH2O(i)+qCH4+qPH3(i));
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
    g0_i=900.0;
    R0e_i=ao;
    P0_i=1.0
    T_targ_i=146.2;
    P_targ_i=1.29848
    P_term_i=0.0631;
    use_lindal='Y';
    n_lindal_pts=27;
    SuperSatSelf_H2S=0.0;
    SuperSatSelf_NH3=0.0;
    SuperSatSelf_PH3=0.0;
    SuperSatSelf_H2O=0.0;
    supersatNH3=0.0;
    supersatH2S=0.0;
    AutoStep_constant=8;
    dz=1;
    TP_list={'Seiff_Jupiter','Lindal_Jupiter','Lindal_Saturn','whatever_is_in_TCM_mex'}
    AutoStep_constant=8;
    use_dz=1;
    dP_init=1;
    dP_fine=0.0001;
    P_fine_start=10;
    P_fine_stop=1;
    frain=0;
    select_ackerman=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Run TCM, and arrange variables for maintam      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for k=1:length(xHe)
    % [me(k),n]=size(P);
    % If cases includes a 10% warming, do it, otherwise copy lindal 1985
    % profile
        if(strcmp(Model_names(k),'Mean_p_warming')|strcmp(Model_names(k),'Dep_Ammonia_p_Warming'))
            TP_loaded=load('TCM_mex/TP.SAT'); % load lindal profile
            TP_modified=[TP_loaded(:,1),TP_loaded(:,2)*1.1];
            dlmwrite('TCM_mex/TP.high',TP_modified,'delimiter','\t','precision','%.4f');
            system('cp TCM_mex/TP.high TCM_mex/TP.TCM')
            TP_force='whatever_is_in_TCM_mex';
        else
            TP_force='Lindal_Saturn';
        end
    
    [me(k),tcme(1:me(k),:,k),tcmp(1:me(k),:,k)]=DeBoer_TCM(TP_list,TP_force,xH2S(k),xNH3(k),xH2O(k),xCH4(k),...
                                xPH3(k),xHe(k),xCO,P_temp,T_temp, g0_i,R0e_i,...
                                P0_i,T_targ_i,P_targ_i,P_term_i,...
                                use_lindal,SuperSatSelf_H2S,SuperSatSelf_NH3,...
                                SuperSatSelf_PH3,SuperSatSelf_H2O,supersatNH3,...
                                supersatH2S,AutoStep_constant,fp,dz,oblateness_factor,use_dz,dP_init,dP_fine,P_fine_start,P_fine_stop,frain,select_ackerman);
        if (k==3)
            for jj=1:me(k)
                if(tcme(jj,1,k)<1.4)
                    tcme(jj,7,k)=0.40*tcme(jj,7,k);
                end
            end
        end
        if (k==4)
            for jj=1:me(k)
                if(tcme(jj,1,k)<3)
                    tcme(jj,7,k)=0.40*tcme(jj,7,k);
                end
            end
        end
        if (k==6)
            for jj=1:me(k)
                if(tcme(jj,1,k)<3)
                    tcme(jj,7,k)=0.10*tcme(jj,7,k);
                end
            end
        end
    end


%Spacecraft Parameters
    Rayorigin=[(number_of_saturn_radii+1)*ao 0 0];
    Raydirection=[-1 0 0];

    no_ph3=1; %use PH3 decay

    for j=1:length(xHe)
        [Tbeam_nadir(j,cases),zenith_nadir(j,cases),weighting_function_a_nadir(1:me(j)+1,j,cases)]= maintamone(Raydirection,Rayorigin,...
                               tcme(1:me(j),1:22,j),tcmp(1:me(j),1:22,j),ao,bo,co,...
                               f,no_ph3,select_h2h2_model,select_ammonia_model,select_water_model,...
                               include_clouds,N_ring_one,Nphi,BWHM,refractivity_source);
    end

    X_direction=-cos(theta*(pi/180));
    Y_direction=0;
    Z_direction=sin(theta*(pi/180));
    Raydirection=[X_direction Y_direction Z_direction];

    for j=1:length(xHe)       
        [Tbeam_limb(j,cases),zenith_limb(j,cases),weighting_function_a_limb(1:me(j)+1,j,cases)]= maintamone(Raydirection,Rayorigin,...
                                    tcme(1:me(j),1:22,j),tcmp(1:me(j),1:22,j),ao,bo,co,...
                                    f,no_ph3,select_h2h2_model,select_ammonia_model,select_water_model,...
                                    include_clouds,N_ring_one,Nphi,BWHM,refractivity_source);
                               clear ph3decay;
                               
                             % Sorry! ugly piece to flip around PH3 so you
                             % can plot it nicely against P.
                               xxxPH3=tcme(1:me(j),10,j);   % Get the value of PH3 from tcm matrix
                               ph3_size=size(xxxPH3,1);     % get the size of the PH3 vector
                               kk=(0:ph3_size-1);           % create a reverse index of kk
                               xxPH3=[0;xxxPH3(ph3_size-kk)]; % flip the phosphine vector and add a zero (hangover from hoffman code)
                               P=tcme(1:me(j),1,j);           % Get the pressure vector
                               Pp=[0;P(ph3_size-kk)];         % add zero and flip Pressure vector
                               decay=ph3decay(Pp,xxPH3);      % run the phosphine decay routine
                               close_ph3=decay(ph3_size-kk);  % flip the phosphine back
                               ph3_flipped=fliplr(close_ph3);  % flip it again
                               tcme(1:me(j),10,j)=ph3_flipped; % save it as a vector
    end
     
    R(:,cases)=100*(Tbeam_nadir(:,cases)-Tbeam_limb(:,cases))./Tbeam_nadir(:,cases);
end
save(output_filename);
