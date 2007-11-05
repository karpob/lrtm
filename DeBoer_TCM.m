function [me,tcme,tcmp,DSOL_NH3]=DeBoer_TCM(TP_list,TP_force,XH2S_i,XNH3_i,XH2O_i,XCH4_i,...
                                XPH3_i,XHe_i,XCO,P_temp,T_temp, g0_i,R0e_i,...
                                P0_i,T_targ_i,P_targ_i,P_term_i,...
                                use_lindal,SuperSatSelf_H2S,SuperSatSelf_NH3,...
                                SuperSatSelf_PH3,SuperSatSelf_H2O,supersatNH3,...
                                supersatH2S,AutoStep_constant,fp,dz,oblateness_factor)
if(strcmp(use_lindal,'Y'))
    lindal_profile_switch=1;
else
    lindal_profile_switch=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copy over desired Temp Pressure profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(strcmp(TP_list(1),TP_force))
    system('cp TCM_mex/TP_Seiff_Jupiter.JUP TCM_mex/TP.TCM');
    TP_in_directory=load('TCM_mex/TP.TCM');
    n_lindal=length(TP_in_directory);    

elseif(strcmp(TP_list(2),TP_force))
    system('cp TCM_mex/TP.JUP TCM_mex/TP.TCM');
    TP_in_directory=load('TCM_mex/TP.TCM');
    n_lindal=length(TP_in_directory);

elseif(strcmp(TP_list(3),TP_force))
    system('cp TCM_mex/TP.SAT TCM_mex/TP.TCM');
    TP_in_directory=load('TCM_mex/TP.TCM');
    n_lindal=length(TP_in_directory);

elseif(strcmp(TP_list(4),TP_force))
    system('cp TCM_mex/TP.URN TCM_mex/TP.TCM');
    TP_in_directory=load('TCM_mex/TP.TCM');
    n_lindal=length(TP_in_directory);

elseif(strcmp(TP_list(5),TP_force))
    system('cp TCM_mex/TP.NEP TCM_mex/TP.TCM');
    TP_in_directory=load('TCM_mex/TP.TCM');
    n_lindal=length(TP_in_directory);
elseif(strcmp(TP_list(6),TP_force))
    TP_in_directory=load('TCM_mex/TP.TCM');
    n_lindal=length(TP_in_directory);
else
    disp('I do not have a TP profile to force. You enjoy failure.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run DeBoer TCM (must be compiled first using:
%  >>mex *.C -o TCM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd TCM_mex
[P,T,XH2,XHe,XH2S,XNH3,XH20,XCH4,XPH3,...
    clouds,DNH4SH_i,DH2S_i,DNH3_i,DH2O_i,DCH4_i,DPH3_i,DSOL_i,...
    g,mu,ref_w_o,ref_w,z,DSOL_NH3]=TCM(dz,XHe_i,XH2S_i,XNH3_i,XH2O_i,XCH4_i,XPH3_i,XCO,...
                              P_temp,T_temp,g0_i,R0e_i,P0_i,T_targ_i,P_targ_i,P_term_i,1,n_lindal,...
                              SuperSatSelf_H2S,SuperSatSelf_NH3,SuperSatSelf_PH3,...
                              SuperSatSelf_H2O,supersatNH3,supersatH2S,...
                              AutoStep_constant,fp);
clear TCM;
[me,n]=size(P);
cd ..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove fictional clouds from DeBoer TCM using built-in filter 'clouds'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[DNH4SH,DH2S,DNH3,DH2O,DCH4,DPH3,DSOL]=filter_clouds(clouds,DNH4SH_i,DH2S_i,DNH3_i,DH2O_i,DCH4_i,DPH3_i,DSOL_i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Scale Polar profile according to oblateness factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Equatorial
tcme(1:me,1:22)=[P(1:me),T(1:me),z(1:me),XH2(1:me),XHe(1:me),XH2S(1:me),XNH3(1:me),XH20(1:me),XCH4(1:me),XPH3(1:me),...
                 clouds(1:me),DNH4SH(1:me),DH2S(1:me),DNH3(1:me),DH2O(1:me),DCH4(1:me),DPH3(1:me),DSOL(1:me),...
                 g(1:me),mu(1:me),ref_w_o(1:me),ref_w(1:me)];
%Polar
tcmp(1:me,1:22)=[tcme(1:me,1:2),oblateness_factor.*tcme(1:me,3),tcme(1:me,4:22)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
