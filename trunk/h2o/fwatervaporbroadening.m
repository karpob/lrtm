function gammawater=fwatervaporbroadening(P_H2O,P_H2,P_He,T)

% function gammawater=fwatervaporbroadening(P_H2O,P_H2,P_He,T)
% P's should be row vectors
% Unit is lower case gamma
% Pressure Broadening term adapted for outer planet use
% From DeBeor thesis, see main water attenuation program for cite
% Checks with "manual" calculations in Excel, see file 'checkh20 (tab waterbroad)

% Finds gamma


% There are contributions to each resonant line (j) from each broadening agent(i)
% Broadening agents are mainly H2, He, and H2O

To=300;				% Reference Temperature 300 Kelvin
Tdiv=To/T;			% Ratio of Ref Temp to actual temp, used often


fo=[22.23515, 183.31012, 323, 325.1538, 380.1968, 390, 436, 438, 442, 448.0008];
Ep=[644, 196, 1850, 454, 306, 2199, 1507, 1070, 1507, 412];
A=[1.0, 41.9, 334.4, 115.7, 651.8, 127.0, 191.4, 697.6, 590.2, 973.1];
g_H2=[2.395, 2.4000, 2.395, 2.395, 2.390, 2.395, 2.395, 2.395, 2.395, 2.395];
g_He=[0.67, 0.71, 0.67, 0.67, 0.63, 0.67, 0.67, 0.67, 0.67, 0.67];
g_H2O=[10.67, 11.64, 9.59, 11.99, 12.42, 9.16, 6.32, 8.34, 6.52, 11.57];

zeta_H2=[0.900, 0.950, 0.900, 0.900, 0.850, 0.900, 0.900, 0.900, 0.900, 0.900];
zeta_He=[0.515, 0.490, 0.515, 0.490, 0.540, 0.515, 0.515, 0.515, 0.515, 0.515];
zeta_H2O=[0.626, 0.649, 0.420, 0.619, 0.630, 0.330, 0.290, 0.360, 0.332, 0.510];

Tdep=[Tdiv.^(zeta_H2);Tdiv.^(zeta_He);Tdiv.^(zeta_H2O)]; % temp dep
gamma_v=[g_H2*P_H2.*Tdep(1,:); g_He*P_He.*Tdep(2,:); g_H2O*P_H2O.*Tdep(3,:)];
gammawater=sum(gamma_v,1);		% The ten line broad params-check to see if works



