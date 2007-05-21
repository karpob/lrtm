function alphanh3=nh3bgmodel(f,T,P,H2mr,Hemr,NH3mr) % Berge-Gulkis model

%input f as a vector in GHz, T in Kelvin, P in Bars, H2mr, Hemr and NH3mr in mole
%fraction, so for a 1% mixture of NH3, NH3mr = 0.01, etc.  f is the only
%variable that can be a vector.  Opacity (alphanh3) is returned as a vector matching f
%in dB/km.

Joiner=0;
load nh3lincatPK16;	
%Loads the nh3 line catalog.  The workspace contains the arrays for Freq (fo) in GHz, intensity (Io) in inverse cm,
%lower state energy (Eo) in inverse cm and self broadening due to NH3 (gammaNH3o) in GHz/bar.
%Constants
Fo=Fo/1000; %convert MHz to GHz
OpticaldepthstodB=434294.5;				%converts from cm^-1 to dB/km

bartoatm=0.986923;

To=300;

P_H2=P*H2mr; % calculate partial pressures
P_He=P*Hemr;
P_NH3=P*NH3mr;
if Joiner==1  % Attempt at recreating Joiner's incorrect interpretation of the Berge & Gulkis model
    PH2=P_H2;
    PHe=P_He;
    PNH3=P_NH3;
    Ampli=1214;
    C=1.0075+0.0304*(PNH3/T)+0.0537*(PH2/T).^2;
    delt=-0.45;
else
    PH2=P_H2*bartoatm;
    PHe=P_He*bartoatm;
    PNH3=P_NH3*bartoatm;
    Ampli=1230;
    C=1.0075+(0.0308+0.0552*(PH2/T)).*(PH2/T);
    delt=-0.45;
end

%calculate vector linewidth
xi=2/3;
xi2=1;
Tdiv=To/T;

%Loads the nh3 line catalog Width parameter in MHz/torr.  

gammaNH3o=Width;%*760/1000; % MHz/torr to GHz/atm (not necessarily)
gnu1=2.318;
gnu2=0.79;
gnu3=0.75;
gH2=gnu1*PH2;
gHe=gnu2*PHe;
gNH3=gnu3*PNH3*gammaNH3o;
gamma=((gH2+gHe)*((Tdiv)^xi))+gNH3*(300/T)^(xi2);

znu1=1.92;
znu2=0.30;
znu3=0.49;
zH2=znu1*PH2;
zHe=znu2*PHe;
zNH3=znu3*PNH3*gammaNH3o;
zeta=((zH2+zHe)*((Tdiv)^xi))+zNH3*(300/T)^(xi2);

zetasize=size(zeta,1);
pst=delt*PNH3*ones(zetasize,1);							

alpha_noshape=Ampli*((2*J+1).*K.^2./(J.*(J+1))).*Fo.^2.*S.*(PNH3)./(T^(7/2)).*exp(4.8/T*(1.09*K.^2-2.98*J.*(J+1)));  

n=size(f,2);  %returns the number of columns in f
m=size(Fo,1); %returns the number of rows in Fo

nones=ones(1,n);
mones=ones(m,1); 
f_matrix=mones*f;
Fo_matrix=Fo*nones;
dnu_matrix=gamma*nones;    
ce_matrix=zeta*nones;
pst_matrix=pst*nones;

Aa=2*((f_matrix./Fo_matrix).^2);			% no gamma needed here as they cancel

Bb=(dnu_matrix-ce_matrix).*(f_matrix.^2);
Cc=dnu_matrix+ce_matrix;
Dd=((Fo_matrix+pst_matrix).^2) + (dnu_matrix.^2)-(ce_matrix.^2);
Ee=f_matrix.^2;
Jj=(Fo_matrix+pst_matrix).^2;
Gg=dnu_matrix.^2;
Hh=ce_matrix.^2;
Ii=4*(f_matrix.^2).*(dnu_matrix.^2);
Fbr=Aa.*(Bb+Cc.*Dd)./(((Ee-Jj-Gg+Hh).^2)+Ii);						% m x n matrix

alpha_noshape_matrix=alpha_noshape*nones;

br_alpha_matrix=alpha_noshape_matrix.*Fbr;

%sums up the element in the matrix to calculate the alpha in optical depths or inverse cm
alpha_opdep=sum(br_alpha_matrix,1); % cm^-1

alphanh3=alpha_opdep.*C*OpticaldepthstodB;
%answer in dB/km

