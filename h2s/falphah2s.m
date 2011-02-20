function alphah2s=falphah2s(f,T,P_H2,P_He,P_H2S)

% 12/12/MM
% function alphaH2S=falphaH2S(f,T,P,xH2,xHe,xH2S)
% T in kelvin, P in bars, xH2,xHe,xH2S = fraction (x/100) of each 

if P_H2==0
   alphah2s=0;
   return
end


load h2slin;	% use for testing-NEED H2S.lin
% ??? is it Converted P/P catalogue with GHz; nm^2MHz; 1/cm (Freqline Int Eo)

%Constants
GHztoinv_cm=1/29.9792458;	% for converting GHz to inverse cm.
OpticaldepthstodB=434294.5;				%/* converts from cm^-1 to dB/km */
torrperatm=760;
atmperbar=0.987;
GHzperMHz=1/1000;
hc=19.858252418E-24;			%	planks (J.s) light (cm/s)
K=1.38*10^-23;					%	boltzmann's in J/K or N.m/K
No=6.02297e23;					% Avogadros Number [mole^-1]
R=8.31432e7;					% Rydberg's [erg/mole-K]
To=300;			%Ref temp for P/P Catalogue



%Convert pressure percents to partial pressures
PH2=P_H2;						% Equations in bar ???
PHe=P_He;						% Equations in bar ???
PH2S=P_H2S;						% Equations in bar ???


%calculate vector linewidth
[gammaH2S,psiH2S] = brparametersh2s(To,T,PH2,PHe,PH2S);
psisize=size(psiH2S,2);
pst_p=pst*PH2S;							% answer in GHz

%convert to inverse cm
d_nu=gammaH2S*GHztoinv_cm;			%Correct for units of d_nu in GHz to (cm^-1)
ce_cm=psiH2S.*GHztoinv_cm;
pst_cm=pst.*GHztoinv_cm;

% Convert Intensity at 300 Kelvin to Intensity at T kelvin:

numberdensity=(10^6)*No*PH2S/(R*T);	% [molec/cm^3]
												% The 10^6 allows use of P(bar) for P(dynes/cm^2)

eta=3/2;										% for symmetric top molecule
expo=-(1/T-1/To)*Ep*hc/K;
ST=(2.997847*10^18)*Io.*exp(expo).*(To/T)^(eta+1);	% S(T) =S(To)converted for temperature

alpha_max=ST.*((d_nu).^-1)*((pi*2.9979E18)^-1)*numberdensity; % 
%%%% Alpha Max Found


%%%% Get Lineshape
% Lineshapes calculated using inverse centimeters
f_cm=f.*GHztoinv_cm;					% Convert observation freq to cm^-1
fo_cm=fo.*GHztoinv_cm;

% DONT NEED VAN VLECK
%Fvvw=vvwlineshape(f_cm,fo_cm,d_nu);	% Matrix of form f(1)*fo(1)  f(2)*fo(1) f(3)*fo(1)
%												  f(1)*fo(2)  f(2)*fo(2) f(3)*fo(2)

Fbr=brlineshape(f_cm,fo_cm,d_nu,ce_cm,pst_cm);
% Fbr=brlineshape(f_cm,fo_cm,d_nu,zed,zed); %Becomes VVW lineshape

% Put alpha_max(Eo) into a matrix m x n like ie. Fvvw
n=size(f,1);
nones=ones(n,1);

alpha_max_matrix=nones*alpha_max;
% Multiply element by element alpha(1)*F(f,lf(1)) etc...
%vvw_alpha_matrix=Fvvw.*alpha_max_matrix;
br_alpha_matrix=Fbr.*alpha_max_matrix;

% Sum each column, which corresponds to one observation frequency.
%ie if 1-10Ghz, first column are the contributions of all line freq at 1 Ghz
%vvw_alpha_neper=sum(vvw_alpha_matrix,1).*(pi);
%vvw_alpha=vvw_alpha_neper.*OpticaldepthstodB;

br_alpha_opdep=sum(br_alpha_matrix,2).*(pi);
%C=-0.337+(T/110.4)-((T^2)/70600);		% Correction factor

alphah2s=br_alpha_opdep;
alphaH2Sdb=br_alpha_opdep.*OpticaldepthstodB;	% Want in optical depths per cm
%alphaH2S=br_alpha_neper.*OpticaldepthstodB*C;
