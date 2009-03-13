function[alpha,alpha_H2_prime,alpha_He_prime,alpha_CH4_prime]=...
        falpha_orton_quantum_h2h2(f,T,P_H2,P_He,P_CH4,normal_or_equilibrium)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function uses data supplied by Glenn Orton from his quantum scattering
% code for H2-H2, H2-He, and H2-CH4. The data is interpolated using interp2
% with spline option (along frequency, and temperature grid).
%     Input
%     f     --> Frequency in GHz
%     T     --> Temperature in K
%     P_H2  --> Patial Pressure from H2 (bars)
%     P_He  --> Partial Pressure from He (bars)
%     P_CH4 --> Partial Pressure from CH4 (bars)
%     normal_or_equilibrium --> select normal (1) or eq(0) hydrogen
%
%     Output
%     alpha <-- total CIA absorption from H2
%     alpha_H2_prime  <-- CIA absorption from H2-H2
%     alpha_He_prime  <-- CIA absorption from H2-He
%     alpha_CH4_prime <-- CIA absorption from H2-CH4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert stuff to wavenumbers and amagats...let the fun begin.%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GHztoHz=1e9;
f=f*GHztoHz;
c=2.99792458e10; %cm/sec
nu=f./c;

Lo=2.68719e19; %Loschmidt number molecules/cm^3 at stp

Pascal_per_bar=100000;
cm3_per_m3=1e6;
grams_per_kg=1000;
molecules_per_mole=6.0221415e23;

RH2=4124.18; %Nm/Kg/K
RHe=2077;
RCH4=518.3;

rho_he=(Pascal_per_bar*P_He)/(RHe*T);
rho_he=rho_he/cm3_per_m3;
rho_he=rho_he*grams_per_kg;

rho_h2=(Pascal_per_bar*P_H2)/(RH2*T);
rho_h2=rho_h2/cm3_per_m3;
rho_h2=rho_h2*grams_per_kg;

rho_ch4=(Pascal_per_bar*P_CH4)/(RCH4*T);
rho_ch4=rho_ch4/cm3_per_m3;
rho_ch4=rho_ch4*grams_per_kg;

grams_per_mole_h2=2*1.007822;
grams_per_mole_He=4.003;
grams_per_mole_CH4=16.04;

rho_h2=rho_h2/grams_per_mole_h2;
rho_h2=rho_h2*molecules_per_mole;
amagat_h2=rho_h2/Lo;
amagat_sq_h2=amagat_h2^2;

rho_he=rho_he/grams_per_mole_He;
rho_he=rho_he*molecules_per_mole;
amagat_he=rho_he/Lo;
amagat_sq_he=amagat_he^2;

rho_ch4=rho_ch4/grams_per_mole_CH4;
rho_ch4=rho_ch4*molecules_per_mole;
amagat_ch4=rho_ch4/Lo;
amagat_sq_ch4=amagat_ch4^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enough having fun with amagats, time to read in the data             %
% Read in data supplied by Glenn Orton from His quantum scattering code%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name='cia.longrev08.data';%'orton_h2h2_grid.txt';
file_handle=fopen(name,'r');
line1=fscanf(file_handle,'%f %f %f',[3,1]);
line2=fscanf(file_handle,'%f',[1,1]);
line_length=line2;
line_format=strcat('%',sprintf('%d', line_length),'f');
wave_numbers=fscanf(file_handle,line_format,[line2,1]);
alpha_eH2=fscanf(file_handle,line_format,[line1(1),line2]);
alpha_nH2=fscanf(file_handle,line_format,[line1(1),line2]);
alpha_eH2_He=fscanf(file_handle,line_format,[line1(1),line2]);
alpha_nH2_He=fscanf(file_handle,line_format,[line1(1),line2]);
alpha_eH2_CH4=fscanf(file_handle,line_format,[line1(1),line2]);
alpha_nH2_CH4=fscanf(file_handle,line_format,[line1(1),line2]);
fclose(file_handle);

normal_or_equilibrium==1;
for i=1:length(wave_numbers)
   if(wave_numbers<0.1)
       scale_factor=100*(wave_numnbers(i)*wave_numbers(i));
   else
       scale_factor=1;
   end
end

% Close the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pick which set of absorption coefficients you want. (normal, or eq)
%
% Absorption coefficient of things with normal hydrogen=1, 
% equilibrium hydrogen=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if(normal_or_equilibrium==0)
%    alpha_H2=alpha_eH2;
%    alpha_He=alpha_eH2_He;
%    alpha_CH4=alpha_eH2_CH4;
%end
%if(normal_or_equilibrium==1)
%    alpha_H2=alpha_nH2;
%    alpha_He=alpha_nH2_He;
%    alpha_CH4=alpha_nH2_CH4;
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Start interpolation scheme                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% "evenly spaced T in log space used for calculation 
%   (ie. 10 temperatures in Orton's file)

%Ti=[40,51.667,66.724,86.177,111.302,143.753,185.664,239.794,309.705,400];
logTi=linspace(log(40),log(400),10);
% interpolate values for given wavenumber (nu), and while we're at it multiply
% by amagat^2 of each constituent
no_div_by_three=0;
np=0;
theta_x=-85.4;
k=0;
fp=0;
no_actual=0;

for i=1:201
    k=i-1;
    no_div_by_three=(2*(2*k+1)+1)*exp(theta_x*(2*k+1)*((2*k+2)/T))+no_div_by_three;
    np=(2*2*k+1)*exp(theta_x*(2*k+1)*((2*k)/T))+np;
end

no_actual=no_div_by_three*3;
np_actual=np;

    alpha_H2_prime_1=amagat_sq_h2*exp(interp2(wave_numbers,logTi,alpha_eH2,nu,log(T),'spline'));
    alpha_He_prime_1=amagat_he*amagat_h2*exp(interp2(wave_numbers,logTi,alpha_eH2_He,nu,log(T),'spline'));
    alpha_CH4_prime_1=amagat_h2*amagat_ch4*exp(interp2(wave_numbers,logTi,alpha_eH2_CH4,nu,log(T),'spline'));
%   alpha_H2_prime_2=amagat_sq_h2*exp(interp2(wave_numbers,logTi,alpha_nH2,nu,log(T),'spline'));
%    alpha_He_prime_2=amagat_he*amagat_h2*exp(interp2(wave_numbers,logTi,alpha_nH2_He,nu,log(T),'spline'));
%    alpha_CH4_prime_2=amagat_h2*amagat_ch4*exp(interp2(wave_numbers,logTi,alpha_nH2_CH4,nu,log(T),'spline'));

%
%    alpha_H2_prime_2=amagat_sq_h2*exp(interp2(wave_numbers,logTi,alpha_nH2,nu,log(T),'spline'));
%    alpha_He_prime_2=amagat_sq_he*exp(interp2(wave_numbers,logTi,alpha_nH2_He,nu,log(T),'spline'));
%    alpha_CH4_prime_2=amagat_sq_ch4*exp(interp2(wave_numbers,logTi,alpha_nH2_CH4,nu,log(T),'spline'));

% Add up alpha's of H2-H2, H2-He, H2-CH4
alpha_1=alpha_H2_prime_1+alpha_He_prime_1+alpha_CH4_prime_1;
%alpha_2=alpha_H2_prime_2+alpha_He_prime_2+alpha_CH4_prime_2;
fp=(np_actual)/(np_actual+no_actual);
%fo=(no_actual)/(np_actual+no_actual);
%if(T>400)
%    alpha=alpha_2;%+fo*alpha_2;%*(np_actual)/(np_actual+no_actual); % ((np_actual)/(np_actual+no_actual))*alpha_1+((no_actual)/(np_actual+no_actual))*alpha_2;
%else
    alpha=alpha_1;
%end
% If nu<0.1 apply correction factor.
%if(nu<0.1)
%  alpha=100*(nu^2)*alpha;
%end
