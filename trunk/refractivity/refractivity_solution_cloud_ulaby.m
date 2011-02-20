function [Np,Npp]=refractivity_solution_cloud_ulaby(T,D_sol,rho,f)
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


% Use dielectric constant from Ulably

Eps=88.045-0.4147.*T+((6.295e-4).*T.*T)+((1.075e-5)*T.*T.*T); %Static dielectric constant of pure H2O
cd ../clouds
C=get_complex_dielectric_constant_water(f,T);
cd ..
Epp=real(C);
Eppp=imag(C);
cd refractivity
D_sol_g_m3=(D_sol*1e6)/rho; % convert to g/m^3 This is necessary, since N', N'' are small parts of total refractive index
% Use Refractivity from Liebe, 1989 (most closely matches his notation in
% paper), but 1993 will work too

Ep=(2+Epp)/Eppp; % What he calls y in his 1989 paper
Nwpp=(4.5*D_sol_g_m3)/(Eppp*(1+Ep^2)); % imaginary part of refractivity
Dw=-4.5*D_sol_g_m3*(Ep/(Eppp*(1+Ep*Ep)))+4.5*D_sol_g_m3/(Eps+2); % real part of refractivity

Np=Dw;
Npp=Nwpp;
