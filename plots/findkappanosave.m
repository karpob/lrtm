function kappa=findkappanosave(f,T,P_H2,P_He,P_NH3,P_H2O,P_CH4,P_PH3,P_H2S)

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

stopindex=size(T,1);
pp=cd;
for k=1:stopindex
   % Call ammonia
   %alphanh3(k)=falphanh3(f,T(k),P_H2(k),P_He(k),P_NH3(k));		% Vector in f
   cd c:\models\rtm\imaging\nh3
   %krap=falphabg;
   %alphanh3=feval(krap,f,T,P_H2,P_He,P_NH3);			% Will need to call for each P(z)
   alphanh3(k)=falphanh3(f,T(k),P_H2(k),P_He(k),P_NH3(k));			% Will need to call for each P(z)
   % Call phosphine
   cd c:\models\rtm\imaging\ph3
   alphaph3(k)=falphaph3(f,T(k),P_H2(k),P_He(k),P_PH3(k));	% Vector in f
   
   
   % Call Hydrogen Sulfide
   cd c:\models\rtm\imaging\h2s
   alphah2s(k)=falphah2s(f,T(k),P_H2(k),P_He(k),P_H2S(k));	% Vector in f
   
   
   % Call Hydrogen
   cd c:\models\rtm\imaging\h2
   alphah2(k)=falphah2(f,T(k),P_H2(k),P_He(k),P_CH4(k));
   
   
   % Call Water Vapor
   cd c:\models\rtm\imaging\h2o
   alphah2o(k)=fwatervaporalpha(f,T(k),P_H2(k),P_He(k),P_H2O(k));
      % Sum up
   
   kappa_1(k)=alphanh3(k)+alphaph3(k)+alphah2(k)+alphah2o(k)+alphah2s(k);
end
%kappa_2=(kappa_1*1e5)';		% Put into optical depths per km
cd (pp)
%save alfs alphanh3 alphaph3 alphah2s alphah2 alphah2o
kappa=[kappa_1'];	% (in 1/cm) if limb-have passed same pressure twice-masterindex keeps track
