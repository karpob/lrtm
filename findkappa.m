function kappa=findkappa(f,T,P,P_H2,P_He,P_NH3,P_H2O,P_CH4,P_PH3,P_H2S,XH2,XHe,XNH3,XH2O,DNH4SH,DH2S,DNH3,DH2O,DSOL,select_h2h2_model,select_ammonia_model,select_water_model,include_clouds)

% function kappa=findkappa(f,T,P_H2,P_He,P_NH3,P_H2O,P_CH4,P_PH3)
% Finds the atmospheric absorption at each pressure level.
% 
%
%
%            VARIABLE DEFINITIONS:
%
%                        -->INPUT:
%                               ->f: Frequency in GHz
%                               ->T; Temperature Profile (K)
%                               ->P: Pressure Profile (bars)
%                               ->P_H2: Partial pressure of hydrogen/H2 (bars)
%                               ->P_He: Partial pressure of helium (bars)
%                               ->P_NH3: Partial pressure of ammonia (bars)
%                               ->P_H2S: Partial pressure of hydrogen sulfide (bars)
%                               ->XH2: absolute Mole fraction of hydrogen/H2 
%                               ->XHe: absolute Mole fraction of He
%                               ->XNH3: absolute Mole fraction ammonia
%                               ->XH2O: absolute mole fraction of water
%                               ->DNH4SH: cloud density of ammonium hydrosulfide cloud (g/cm^3)
%                               ->DNH3: cloud density of ammonia ice (g/cm^3)
%                               ->DH2O: cloud density of water ice (g/cm^3)
%                               ->DSOL: cloud density of water ammonia solution (g/cm^3)
%                               ->select_h2h2_model: select hydrogen continuum absorption model (1=Joiner Model, 2=Goodman, 3=Goodman by Joiner, 4=Borysow, 5=Borysow Orton modification
%                               ->select_ammonia_model: select ammonia absorption model (1=Spilker original Hoffman coding,2=Joiner Steffes 3=Berge Gulkis 4=Mohammed Steffes KA_band 5=Spilker Model
%                               ->select_water_model: select water absorption model (1=DeBoer Original, 2=DeBoer Corrected, 3=Goodman)
%                               ->include_clouds: select whether or not to include cloud absorption (1=Yes, 0=No)
%                        <--OUTPUT:
%                               <-kappa:layer(s) absorption coefficient in cm^-1
% Output shoud be in optical depths per centimeter

% THIS VERSION USED FOR IMAGING DOES ALL PRESSURE LEVELS-TO BE USED MULTIPLE TIMES FOR
% THE DIFFERENT LOOK ANGLES

% Includes the opcacity due to Ammonia, Phosphine, Hydrogen Sulfide, Hydrogen, Water Vapor
% inlcuding H2 collisions with CH4
%maxmaster=max(masterindex)			% masterindex may go up and down-only
%need to calc repeaters once

OpticaldepthstodB=434294.5;
stopindex=size(T,1);
stopindex=stopindex-1;
 % can't include bottom P pressure..no appropriate depth associated with it.
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
    if(select_ammonia_model==6)
       alphanh3(k)=nh3hsmodel(f,T(k),P(k),XH2(k),XHe(k),XNH3(k))/OpticaldepthstodB;
    end
    if(select_ammonia_model==7)
       alphanh3(k)=nh3hsmodel2(f,T(k),P(k),XH2(k),XHe(k),XNH3(k))/OpticaldepthstodB;
    end
   if(select_ammonia_model==8)
      alphanh3(k)=devaraj_nh3_model_0(f,T(k),P(k),XH2(k),XHe(k),XNH3(k))/OpticaldepthstodB;
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
   if(select_h2h2_model==1)
       alphah2(k)=falphah2(f,T(k),P_H2(k),P_He(k),P_CH4(k));
   end
   if(select_h2h2_model==2)
       alphah2(k)=falphah2_goodman(f,T(k),P_H2(k),P_He(k));
   end
   if(select_h2h2_model==3)
       alphah2(k)=falphah2_goodman_by_joiner(f,T(k),P_H2(k),P_He(k));
   end
   if(select_h2h2_model==4)
       alphah2(k)=falpha_borysow(f,T(k),P_H2(k));
   end
   if(select_h2h2_model==5)
       alphah2(k)=falpha_orton(f,T(k),P_H2(k));
   end
   if(select_h2h2_model==6)
       alphah2(k)=falpha_orton_quantum_h2h2(f,T(k),P_H2(k),P_He(k),P_CH4(k),0);
       % by default we use equilibruim H2, last arg set to zero you want something different...(ie. normal, intermediate, etc change this)
   end
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
   if(select_water_model==4)
      density_h2_g_m3=(P_H2(k)*2.01594)/((8.314472e-5)*T(k));
      density_he_g_m3=(P_He(k)*4.0026)/((8.314310e-5)*T(k));
      density_h2o_g_m3=(P_H2O(k)*(8.314472/0.46141805))/((8.314472e-5)*T(k));
      alphah2o(k)=karpowicz_h2o_model(f,density_h2o_g_m3,density_h2_g_m3,density_he_g_m3,T(k));
   end
end
cd ..
if(include_clouds==1)
    cd clouds
    for k=1:stopindex
        complex_dielectric_SOL(k)=get_complex_dielectric_constant_water(f,T(k)); %Use equations in Ulaby, Fung, Moore for water
        complex_dielectric_H2S(k)=1+sqrt(-1)*1; %Don't care no H2S clouds here.
        complex_dielectric_NH3(k)=(((1.48)^2)-((8.73e-4)^2)) + sqrt(-1)*(2*1.48*8.73e-4); %Use Howett, Carlson, Irwin,Calcutt ref index nu=1300 cm^-1
        complex_dielectric_H2O(k)=3.15 + sqrt(-1)*1e-3;                                   %Wost Case at 33 GHz Matsuoka,Fujita,Mae 1996
        complex_dielectric_NH4SH(k)=(((2.72)^2)-((7.83e-4)^2)) + sqrt(-1)*(2*2.72*7.83e-4);%Use Howett, Carlson, Irwin,Calcutt ref index nu=1300 cm^-1
    
        alpha_NH4SH(k)=rayleigh_absorption(f,DNH4SH(k),1,complex_dielectric_NH4SH(k));
        alpha_H2S(k)=rayleigh_absorption(f,DH2S(k),1,complex_dielectric_H2S(k));
        alpha_NH3(k)=rayleigh_absorption(f,DNH3(k),1,complex_dielectric_NH3(k));
        alpha_H2O(k)=rayleigh_absorption(f,DH2O(k),1,complex_dielectric_H2O(k));
        alpha_SOL(k)=rayleigh_absorption(f,DSOL(k),1,complex_dielectric_SOL(k));
    end
    cd .. 
end

    
% Sum up and get total kappa
if(include_clouds==1)
    for k=1:stopindex
        kappa_1(k)=alphanh3(k)+alphaph3(k)+alphah2(k)+alphah2o(k)+alphah2s(k)+alpha_NH4SH(k)+alpha_H2S(k)+alpha_NH3(k)+alpha_H2O(k)+alpha_SOL(k);
    end
end
if(include_clouds==0)
    for k=1:stopindex
        kappa_1(k)=alphanh3(k)+alphaph3(k)+alphah2(k)+alphah2o(k)+alphah2s(k);
    end
end

%save alfs
kappa=[kappa_1'];	% (in 1/cm) if limb-have passed same pressure twice-masterindex keeps track
save kappa
%kappa_old=kappa;
%save old_kappa kappa_old
%f_old=f;
%T_old=T;
%P_old=P;
%P_H2_old=P_H2;
%P_He_old=P_He;
%P_NH3_old=P_NH3;
%P_H2O_old=P_H2O;
%P_CH4_old=P_CH4;
%P_PH3_old=P_PH3;
%P_H2S_old=P_H2S;
%XH2_old=XH2;
%XHe_old=XHe;
%XNH3_old=XNH3;
%XH2O_old=XH2O;
%select_ammonia_model_old=select_ammonia_model;
%select_water_model_old=select_water_model;
%save old_kappa_conditions f_old T_old P_old P_H2_old P_He_old P_NH3_old P_H2O_old P_CH4_old P_PH3_old P_H2S_old XH2_old XHe_old XNH3_old XH2O_old select_ammonia_model_old select_water_model_old;

