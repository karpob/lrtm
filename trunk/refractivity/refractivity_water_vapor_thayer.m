function N=refractivity_water_vapor_thayer(T,P)
% This function gives the value for refractivity (N) for a given Temperature and Pressure
%     P--> is in Bars
%     T--> is in Kelvins
%     N <-- N-units (unitless)
K3=3.776e5; % Value given by Thayer, 1974
t=T-273.15; % get value in deg C
Pmb=P.*1000; % convert to millibars
Zw=1+1650*(Pmb./(T.*T.*T)).*(1-0.01317*t+1.74e-4*t.*t+1.44e-6*t.*t.*t); % inverse compressibility factor
N=K3.*(Pmb./(T.*T)).*Zw; %Refractivity value
