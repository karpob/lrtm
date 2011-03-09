clear all;
oblateness_factor=0.935;
%oblateness_factor=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Planet verticies                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ao=71492*100000; % along x (in centimeters)
bo=ao;       % along y
co=ao*oblateness_factor; % along z

Model_Names={'Mean_Lindal','Mean_Seiff','Depleted_Ammonia', 'Enhanced_Ammonia',...
             'Depleted_Water','Enhanced_Water'}



% Temperature forcing
%TP_list={'Seiff_Jupiter','Lindal_Jupiter','Lindal_Saturn','whatever_is_in_TCM_mex'}

TP_force='Lindal_Jupiter';

TP_list={'Seiff_Jupiter','Lindal_Jupiter','Lindal_Saturn','whatever_is_in_TCM_mex'}

output_filename='generatedProfiles.mat'

%Abundances
XH2S_rel_H2=7.7e-5;
XNH3_rel_H2=[4.55e-4,4.55e-4,2.0e-4,7.1e-4,4.55e-4,4.55e-4];
XH2O_rel_H2=[6.3838e-3,6.3838e-3,6.3838e-3,6.3838e-3,2.5535e-3,1.0214e-2];
XCH4_rel_H2=2.1e-3;
XPH3_rel_H2=6e-7;
XHe_rel_H2=0.157;
for i=1:length(XNH3_rel_H2)
    XH2(i)=1/(1+XHe_rel_H2+XH2S_rel_H2+XNH3_rel_H2(i)+XH2O_rel_H2(i)+ XCH4_rel_H2+ XPH3_rel_H2);

    %Get Mole Fractions
    XH2S(i)=XH2S_rel_H2*XH2(i);
    XNH3(i)=XNH3_rel_H2(i)*XH2(i);
    XH2O(i)=XH2O_rel_H2(i)*XH2(i);
    XCH4(i)=XCH4_rel_H2*XH2(i);
    XPH3(i)=XPH3_rel_H2*XH2(i);
    XHe(i)=XHe_rel_H2*XH2(i);
end

%Misc DeBoer TCM inputs
%Guess for deep P,T to match 1 bar level.
P_temp=1000;
%T_temp=1.555678710937500e+03;
T_temp=1583.782470703125;
%P_temp=6000;
%T_temp=2200;

g0_i=2330; %2417;
R0e_i=ao;
P0_i=1;
T_targ_i=166;
P_targ_i=1;
P_term_i=0.141;
use_lindal='Y';
SuperSatSelf_H2S=0;
SuperSatSelf_NH3=0;
SuperSatSelf_PH3=0;
SuperSatSelf_H2O=0;
supersatNH3=0;
supersatH2S=0;
AutoStep_constant=8;
fp=666;
dz=1;
XCO=0;
use_dz=0;
dP_init=10;
dP_fine=0.25;
P_fine_start=13;
P_fine_stop=1;
frain=3;
select_ackerman=0;

 table_output=[XH2;XHe;(1e6)*XH2S;(1e6)*XNH3;(1e6)*XH2O;(1e6)*XCH4;(1e6)*XPH3];
    to_dlm=transpose(table_output);
    dlmwrite('mole_fractions_jupiter.dat',to_dlm,'delimiter','&','precision','%.4f')

for i=1:length(XNH3)

    
    case_select=i; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thermo-chemical Modeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [me,tcme,tcmp,DSOL_NH3]=DeBoer_TCM(TP_list,TP_force,XH2S(case_select),XNH3(case_select),XH2O(case_select),XCH4(case_select),...
                                XPH3(case_select),XHe(case_select),XCO,P_temp,T_temp, g0_i,R0e_i,...
                                P0_i,T_targ_i,P_targ_i,P_term_i,...
                                use_lindal,SuperSatSelf_H2S,SuperSatSelf_NH3,...
                                SuperSatSelf_PH3,SuperSatSelf_H2O,supersatNH3,...
                                supersatH2S,AutoStep_constant,fp,dz,oblateness_factor,use_dz,dP_init,dP_fine,P_fine_start,P_fine_stop,frain,select_ackerman);

    Data(i,:)={Model_Names(i),tcme};                        
end
save(output_filename);
