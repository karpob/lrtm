function kappa=findkappa(f,T,P,P_H2,P_He,P_NH3,P_H2O,P_CH4,P_PH3,P_H2S,XH2,XHe,XNH3,XH2O,select_ammonia_model,select_water_model)

% function kappa=findkappa(f,T,P_H2,P_He,P_NH3,P_H2O,P_CH4,P_PH3)
% Finds the atmospheric absorption along the raypath
% File structure of TCM.out
% P T XH2 XHe XH2S XNH3 XH2O XCH4 XPH3  clouds DNH4SH DH2S DNH3 DH2O DCH4 DPH3 DSOL g mu refr_w/o refr_ w/
% Passes in should be in T(Kelvin), P(bars)
% Some models require different units, but will be interpreted internally
% Output shoud be in optical depths per centimeter
% This puts out optical depths per kilometer
% THIS VERSION USED FOR IMAGING DOES ALL PRESSURE LEVELS-TO BE USED MULTIPLE TIMES FOR
% THE DIFFERENT LOOK ANGLES

% Includes the opcacity due to Ammonia, Phosphine, Hydrogen Sulfide, Hydrogen, Water Vapor
% inlcuding H2 collisions with CH4
%maxmaster=max(masterindex)			% masterindex may go up and down-only need to calc repeaters once
OpticaldepthstodB=434294.5;
stopindex=size(T,1);

% Call ammonia
cd NH3
for k=1:stopindex
    if(select_ammonia_model==1)
        alphanh3(k)=falphanh3(f,T(k),P_H2(k),P_He(k),P_NH3(k));			% Will need to call for each P(z)
    end
    if(select_ammonia_model==2)
        alphanh3(k)=nh3jsmodel(f,T(k),P(k),XH2(k),XHe(k),XNH3(k))/OpticaldepthstodB;
    end
    if(select_ammonia_model==3)
        alphanh3(k)=nh3bgmodel(f,T(k),P(k),XH2(k),XHe(k),XNH3(k))/OpticaldepthstodB;
    end
    if(select_ammonia_model==4)
        alphanh3(k)=nh3msmodel_Ka_band(f,T(k),P(k),XH2(k),XHe(k),XNH3(k))/OpticaldepthstodB;
    end
    if(select_ammonia_model==5)
        alphanh3(k)=nh3spilkermodel(f,T(k),P(k),XH2(k),XHe(k),XNH3(k))/OpticaldepthstodB;
    end
end
cd ..

% Call phosphine
cd ph3
for k=1:stopindex
   alphaph3(k)=falphaph3(f,T(k),P_H2(k),P_He(k),P_PH3(k));	% Vector in f
end   
cd ..
% Call Hydrogen Sulfide
cd h2s
for k=1:stopindex
   alphah2s(k)=falphah2s(f,T(k),P_H2(k),P_He(k),P_H2S(k));	% Vector in f
end   
cd ..
% Call Hydrogen
cd h2
for k=1:stopindex
   alphah2(k)=falphah2(f,T(k),P_H2(k),P_He(k),P_CH4(k));
end   
cd ..
% Call Water Vapor
cd h2o
for k=1:stopindex
   if(select_water_model==1)
        alphah2o(k)=fwatervaporalpha_deboer_original(f,T(k),P_H2(k),P_He(k),P_H2O(k));
   end
   if(select_water_model==2)
        alphah2o(k)=fwatervaporalpha_deboer_corrected(f,T(k),P_H2(k),P_He(k),P_H2O(k));
   end
   if(select_water_model==3)
       alphah2o(k)=fwatervaporalpha_goodman(f,P(k),T(k),XH2O(k),XH2(k),XHe(k));
   end
end
cd ..
% Sum up and get total kappa
for k=1:stopindex
    kappa_1(k)=alphanh3(k)+alphaph3(k)+alphah2(k)+alphah2o(k)+alphah2s(k);
end

save alfs alphanh3 alphaph3 alphah2s alphah2 alphah2o
kappa=[kappa_1'];	% (in 1/cm) if limb-have passed same pressure twice-masterindex keeps track
kappa_old=kappa;
save old_kappa kappa_old
f_old=f;
T_old=T;
P_old=P;
P_H2_old=P_H2;
P_He_old=P_He;
P_NH3_old=P_NH3;
P_H2O_old=P_H2O;
P_CH4_old=P_CH4;
P_PH3_old=P_PH3;
P_H2S_old=P_H2S;
XH2_old=XH2;
XHe_old=XHe;
XNH3_old=XNH3;
XH2O_old=XH2O;
select_ammonia_model_old=select_ammonia_model;
select_water_model_old=select_water_model;
save old_kappa_conditions f_old T_old P_old P_H2_old P_He_old P_NH3_old P_H2O_old P_CH4_old P_PH3_old P_H2S_old XH2_old XHe_old XNH3_old XH2O_old select_ammonia_model_old select_water_model_old;

