function alphanh3=nh3hsmodel(f,T,P,H2mr,Hemr,NH3mr) % Hanley-Steffes model 

%input f as a vector in GHz, T in Kelvin, P in Bars, H2mr, Hemr and NH3mr in mole
%fraction, so for a 1% mixture of NH3, NH3mr = 0.01, etc.  f is the only
%variable that can be a vector.  Opacity (alphanh3) is returned as a vector matching f
%in dB/km.

load nh3lincat190Latest;	
%Loads the Poynter Pickett (JPL) nh3 line catalog.  The workspace contains the arrays for Freq (fo) in GHz, intensity (Io) in inverse cm,
%lower state energy (Eo) in inverse cm and self broadening due to NH3 (gammaNH3o) in GHz/bar.
%Constants
GHztoinv_cm=1/29.9792458;	            %for converting GHz to inverse cm
OpticaldepthstodB=434294.5;				%converts from cm^-1 to dB/km
torrperatm=760;
bartoatm=0.987;
GHztoMHz=1000;
hc=19.858252418E-24;			%planks (J.s) light (cm/s)
k=1.38*10^-23;					%boltzmann's in J/K or N.m/K
No=6.02297e23;					%Avogadros Number [mole^-1]
R=8.31432e7;					%Rydberg's [erg/mole-K]
To=300;			                %Ref temp for P/P Catalogue
dynesperbar=1e6;				%dyne=bar/1e6;
coef=dynesperbar*No/R;          %See Appendix D: Using the Poyter-Pickett Catalogs

PH2=P*H2mr; %partial pressures
PHe=P*Hemr;
PNH3=P*NH3mr;
%1.495;%-0.34;%0.296;
%deltH2=0.0643;
%calculate vector linewidth
xi1=0.7725; %0.7723;%0.7723;%0.7893;%0.934;%2/3;
xi2=2/3;
xi3=1;
xi12=0.7755; %0.7693;%0.8161;%0.825;
xi22=0.7761; %0.7699;%2/3;%0.76;
xi32=1.5521; %1.7271;%1;
Tdiv=To/T;
gnu1=1.7225; %1.7077;%1.7119;%1.647;
gnu2=0.5803; %0.7117;%0.76;
gnu3=0.8329; %0.8433;%0.823;%0.8665;%0.767;%0.823;
gH2=gnu1*PH2;
gHe=gnu2*PHe;
% for inc=1:length(J)
% gammaNH3o(inc)=(C1*(K(inc)^2/(J(inc)*(J(inc)+1)))^Exp1+C2*(K(inc)^2/(J(inc)*(J(inc)+1)))^Exp2)*Kscale(inc)^Exp3;
% gNH3(inc)=gnu3*PNH3*gammaNH3o(inc);
% gamma(inc)=(gH2)*((Tdiv)^(xi1))+(gHe)*((Tdiv)^(xi2))+gNH3(inc)*(300/T)^(xi3);
%gammaNH3o=(C1*(K(inc)^2/(J(inc)*(J(inc)+1)))^Exp1+C2*(K(inc)^2/(J(inc)*(J(inc)+1)))^Exp2)*Kscale(inc)^Exp3;
gNH3=gnu3*PNH3*gammaNH3o;
gamma=(gH2)*((Tdiv)^(xi1))+(gHe)*((Tdiv)^(xi2))+gNH3*(300/T)^(xi3);

delt=0.0048*gamma; %0.0041*gamma;
znu1=1.0941; %1.1561;%1.3256;%1.377;
znu2=1.8449; %1.4414;%0.32;
znu3=0.5468; %0.5463;%0.5737;%0.6131;%0.5736;
zH2=znu1*PH2;
zHe=znu2*PHe;
zNH3=znu3*PNH3*gammaNH3o;
zeta=(zH2)*((Tdiv)^(xi12))+(zHe)*((Tdiv)^(xi22))+zNH3*(300/T)^(xi32);

zetasize=size(fo,1);
pst=delt*PNH3;%*ones(zetasize,1);							% answer in GHz
%Coupling element, pressure shift and dnu or gamma are in GHz, need to convert brlineshape to inverse cm which is done below

n=size(f,2);  %returns the number of columns in f
m=size(fo,1); %returns the number of rows in fo
% f1 f2 f3 f4 ....fn  n times where n is the number of frequency steps
% f1 f2 f3 f4 ....fn				in the observation range                            
% ...
% f1 f2 f3 f4 ....fn 
% m times where m is the number of spectral lines

nones=ones(1,n);
mones=ones(m,1); 
f_matrix=mones*f;
fo_matrix=fo*nones;

% The 10^6 allows use of P(bar) for P(dynes/cm^2)

eta=3/2;			% for symmetric top molecule
expo=-(1/T-1/To)*Eo*hc/k;
ST=Io.*exp(expo);	% S(T) =S(To)converted for temperature
Con=0.9404;%(0.9349+PH2/T*0.5388-(PH2/T)^2*0.102);%(0.9349+PH2/T*0.535+(PH2/T)^2*0.0312);%(0.9404);
alpha_noshape=Con*coef*(PNH3/To)*((To/T)^(eta+2)).*ST;%0.9387  
%Alpha Max Found

%Ben Reuven lineshape calculated by the brlineshape function gives the answer in GHz
%Here we change from GHz to inverse cm.
lineshape=0;
dnu_matrix=gamma*nones;    
ce_matrix=zeta*nones;
pst_matrix=pst*nones;
if lineshape<0.5
Aa=(2/pi)*((f_matrix./fo_matrix).^2);			

Bb=(dnu_matrix-ce_matrix).*(f_matrix.^2);
Cc=dnu_matrix+ce_matrix;
Dd=((fo_matrix+pst_matrix).^2) + (dnu_matrix.^2)-(ce_matrix.^2);
Ee=f_matrix.^2;
Jj=(fo_matrix+pst_matrix).^2;
Gg=dnu_matrix.^2;
Hh=ce_matrix.^2;
Ii=4*(f_matrix.^2).*(dnu_matrix.^2);
Ff=Aa.*(Bb+Cc.*Dd)./(((Ee-Jj-Gg+Hh).^2)+Ii);	
else
Aa=(2/pi)*2*(f_matrix.^2).*dnu_matrix;			

Bb=(fo_matrix.^2-f_matrix.^2).^2;
Cc=4*dnu_matrix.^2.*f_matrix.^2;
Ff=Aa./(Bb+Cc);  
end

Fbr=(1/GHztoinv_cm).*Ff;

alpha_noshape_matrix=alpha_noshape*nones;
%con1=1.7197;
%con2=0.0619;
%con3=5.4015;
%con4=0.2660;
%D=con1.*(f_matrix.^con2)./((con3+abs(f_matrix-fo_matrix)).^con4);

br_alpha_matrix=alpha_noshape_matrix.*Fbr;%.*D;

% foRot=[140141.8974 572498.0678 1168451.6877 1214858.6009 1215245.1833 1763525.3868 1763601.8638 1763821.3719 1808935.5500 1810377.7915 2357210.3572 2357726.7497 2358563.2307 2400017.6324 2400578.3942 2402264.8766 2405121.2992 2948410.6478 2948669.3987 2949480.4272 2950814.5336]'/1000;
% IoRot=[2.9612788E-24 1.0606676E-20 2.9273813E-20 8.2904728E-20 3.1664986E-20 2.0686132E-19 9.3579193E-20 6.1713375E-20 9.8487088E-20 6.5024841E-20 1.6891931E-19 1.4258679E-19 1.8191952E-19 3.6718024E-19 1.7525858E-19 1.4807417E-19 1.8922534E-19 4.5287599E-19 2.2129892E-19 2.0430514E-19 3.4039196E-19]';
% EoRot=[984.1708 0.3967 16.5667 19.4932 15.7763 60.0163 56.3125 45.1906 55.542 44.3993 115.8816 104.7871 86.2612 118.8411 115.1399 104.0254 85.465 198.8972 195.2146 184.1564 165.6913]';
% expoRot=-(1/T-1/To)*EoRot*hc/k;
% STRot=IoRot.*exp(expoRot);	
% gammaNH3oRot=[25 0 16.3 0 16.3 0 12.5 20.45 12.5 20.45 11.2 17.3 22.45 0 11.2 17.3 22.45 0 9.8 14.75 19.35]';
% 
% 
% gNH3rot=gnu3*PNH3*gammaNH3oRot;
% 
% gammaRot=(gH2)*((Tdiv)^(xi1))+(gHe)*((Tdiv)^(xi2))+gNH3rot*(300/T)^(xi3);
% 
% n=size(f,2);  %returns the number of columns in f
% m=size(foRot,1);
% nones=ones(1,n);
% mones=ones(m,1); 
% f_matrix=mones*f;
% fo_matrix=foRot*nones;
% dnu_matrix=gammaRot*nones;
% Aa=4*(f_matrix.^2).*dnu_matrix;			
% Bb=(fo_matrix.^2-f_matrix.^2).^2;
% Cc=4*f_matrix.^2.*dnu_matrix.^2;
% FRot=Aa./(Bb+Cc);
% FbrRot=(1/GHztoinv_cm).*FRot;
% 
% alpha_rot=1.9387*coef*(PNH3/To)*((To/T)^(eta+2))*STRot*nones.*FbrRot;  
alpha_opdep=sum(br_alpha_matrix,1);%+sum(alpha_rot,1);

%sums up the element in the matrix to calculate the alpha in optical depths or inverse cm
%alpha_opdep=sum(br_alpha_matrix,1);

%C=1;
alphanh3=alpha_opdep*434294.5;%*0.907;
%answer in inverse cm converted to dB/km


