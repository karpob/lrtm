function alphanh3=nh3jsmodel(f,T,P,H2mr,Hemr,NH3mr) % Joiner-Steffes model

%input f as a vector in GHz, T in Kelvin, P in Bars, H2mr, Hemr and NH3mr in mole
%fraction, so for a 1% mixture of NH3, NH3mr = 0.01, etc.  f is the only
%variable that can be a vector.  Opacity (alphanh3) is returned as a vector matching f
%in dB/km.

load nh3lincatPK16;	
%Loads the nh3 line catalog.  The workspace contains the arrays for Freq (fo) in GHz, intensity (Io) in inverse cm,
%lower state energy (Eo) in inverse cm and self broadening due to NH3 (gammaNH3o) in GHz/bar.
%Constants
fo=Fo/1000; %convert MHz to GHz
OpticaldepthstodB=434294.5;				%converts from cm^-1 to dB/km
To=300;

PH2=P*H2mr;
PHe=P*Hemr;
PNH3=P*NH3mr;

%calculate vector linewidth
xi=2/3;
xi2=1;
Tdiv=To/T;

gammaNH3o=Width*760/1000;
gnu1=1.69;
gnu2=0.75;
gnu3=0.6;
gH2=gnu1*PH2;
gHe=gnu2*PHe;
gNH3=gnu3*PNH3*gammaNH3o;
gamma=((gH2+gHe)*((Tdiv)^xi))+gNH3*(300/T)^(xi2);

znu1=1.35;
znu2=0.30;
znu3=0.20;
zH2=znu1*PH2;
zHe=znu2*PHe;
zNH3=znu3*PNH3*gammaNH3o;
zeta=((zH2+zHe)*((Tdiv)^xi))+zNH3*(300/T)^(xi2);

zetasize=size(fo,1);
pst=-0.45*PNH3*ones(zetasize,1);							% answer in GHz

alpha_noshape=1214*((2*J+1).*K.^2./(J.*(J+1))).*fo.^2.*S.*(PNH3)./(T^(7/2)).*exp(4.8/T*(1.09*K.^2-2.98*J.*(J+1)));

n=size(f,2);  %returns the number of columns in f
m=size(fo,1); %returns the number of rows in fo

nones=ones(1,n);
mones=ones(m,1); 
f_matrix=mones*f;
fo_matrix=fo*nones;
dnu_matrix=gamma*nones;    
ce_matrix=zeta*nones;
pst_matrix=pst*nones;


Aa=2*((f_matrix./fo_matrix).^2);			

Bb=(dnu_matrix-ce_matrix).*(f_matrix.^2);
Cc=dnu_matrix+ce_matrix;
Dd=((fo_matrix+pst_matrix).^2) + (dnu_matrix.^2)-(ce_matrix.^2);
Ee=f_matrix.^2;
Jj=(fo_matrix+pst_matrix).^2;
Gg=dnu_matrix.^2;
Hh=ce_matrix.^2;
Ii=4*(f_matrix.^2).*(dnu_matrix.^2);
Fbr=Aa.*(Bb+Cc.*Dd)./(((Ee-Jj-Gg+Hh).^2)+Ii);						% m x n matrix

alpha_noshape_matrix=alpha_noshape*nones;

br_alpha_matrix=alpha_noshape_matrix.*Fbr;

%sums up the element in the matrix to calculate the alpha in optical depths or inverse cm

B00=0.07*T*(1-exp(-28.6/T));
B10=0.075*T*(1-exp(-57.2/T))*exp(-28.6/T);
B11=0.053*T*(1-exp(-57.2/T))*exp(-23.3/T);
fo00=595.987; %=19.88*(J+1) cm^-1
fo10=1191.975; %=19.88*(J+1) cm^-1
fo11=1191.975; %=19.88*(J+1) cm^-1

% JS parameters
GH2=1.69; 
GHe=0.75;
GNH3=0.6;

% % BG parameters
% GH2=2.318; 
% GHe=0.79;
% GNH3=0.75;

gamma00=GH2*PH2*(300/T)^(2/3)+GHe*PHe*(300/T)^(2/3)+GNH3*PNH3*(300/T)*14*750.06/1000;
gamma10=GH2*PH2*(300/T)^(2/3)+GHe*PHe*(300/T)^(2/3)+GNH3*PNH3*(300/T)*24*750.06/1000;
gamma11=GH2*PH2*(300/T)^(2/3)+GHe*PHe*(300/T)^(2/3)+GNH3*PNH3*(300/T)*24*750.06/1000;

fo_matrix=fo00*nones;
dnu_matrix=gamma00*nones;
Aa=4*(f.^2).*dnu_matrix;			
Bb=(fo_matrix.^2-f.^2).^2;
Cc=4*f.^2.*dnu_matrix.^2;
F00=Aa./(Bb+Cc);

fo_matrix=fo10*nones;
dnu_matrix=gamma10*nones;
Aa=4*(f.^2).*dnu_matrix;			
Bb=(fo_matrix.^2-f.^2).^2;
Cc=4*f.^2.*dnu_matrix.^2;
F10=Aa./(Bb+Cc);

fo_matrix=fo11*nones;
dnu_matrix=gamma11*nones;
Aa=4*(f.^2).*dnu_matrix;			
Bb=(fo_matrix.^2-f.^2).^2;
Cc=4*f.^2.*dnu_matrix.^2;
F11=Aa./(Bb+Cc);

A00=1.826e3*fo00^2*PNH3/(T^(7/2))*B00;
A10=1.826e3*fo10^2*PNH3/(T^(7/2))*B10;
A11=1.826e3*fo11^2*PNH3/(T^(7/2))*B11;
alpha_opdep=sum(br_alpha_matrix,1)+F00.*A00+F10.*A10+F11.*A11;

alphanh3=alpha_opdep*OpticaldepthstodB;
%answer in dB/km

