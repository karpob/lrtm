function [Np,Npp]=refractivity_solution_cloud_89(T,D_sol,rho,f)
% Calculates the real (Np) and imaginary (Npp) parts of refractvitiy
%
%     Input Variables:
%                  -->T: Temperature in Kelvins
%                  -->D_Sol: Cloud Density in g/cm^3;
%                  -->rho: Density of Cloud Material in g/cm^3
%                  -->f: frequency in GHz
%
%
%     Output Variables:
%                  <--Np: Nprime or real part of refractivity
%                  <--Npp: N double prime or imaginary part of refractivity

D_sol_g_m3=(D_sol*1e6)/rho; % convert to g/m^3 This is necessary, since N', N'' are small parts of total refractive index
% Use Refractivity from Liebe, 1989 (most closely matches his notation in
% paper), but 1993 will work too

V=300/T;
fd=20.09-142*(V-1)+294*(V-1)^2;
Nu=(f/fd);
fs=590-1500*(V-1);
Ns=f/fs;
Epinf=5.48;
Eopt=3.51;
Eps=103.3*(V-1)+77.66;
Epp=Eopt+(Eps-Epinf)/(1+Nu^2)+(Epinf-Eopt)/(1+Ns^2);
Eppp=((Eps-Epinf)*Nu)/(1+Nu^2)+((Epinf-Eopt)*Ns)/(1+Ns^2);
Ep=(2+Epp)/Eppp;
Nwpp=4.5*D_sol_g_m3/(Eppp*(1+Ep^2));
Dw=-4.5*D_sol_g_m3*(Ep/(Eppp*(1+Ep*Ep)))+4.5*D_sol_g_m3/(Eps+2);
Np=Dw;
Npp=Nwpp;