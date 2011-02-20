function alphaph3=falphaph3(f,T,P_H2,P_He,P_PH3);% Revised 11/29/MM


% function alphaph3=falphaph3(f,T,P_H2,P_He,P_PH3);
% A CHANGE, the pressures passed to this are already in partial pressures


load catalog.dat
load a1
load a2


%Constants
GHztoinv_cm=1/29.9792458;	% for converting GHz to inverse cm.
OpticaldepthstodB=434294.5; 				%/* converts from cm^-1 to dB/km */
torrperatm=760;
atmperbar=0.987;
GHzperMHz=1/1000;
hc=19.858252418E-24;			%	planks (J.s) light (cm/s)
K=1.38*10^-23;					%	boltzmann's in J/K or N.m/K
No=6.02297e23;					% Avogadros Number [mole^-1]
R=8.31432e7;					% Rydberg's [erg/mole-K]
To=300;			%Ref temp for P/P Catalogue



% Parameters
%Enter user inputs
%fprintf(1,'Using alphashape \n')
%fprintf(1,'frequency is %6.2f GHz \n',f)
%fprintf(1,'Hydrogen mixing ratio is %6.2f \n',xH2)
%fprintf(1,'Helium mixing ratio is %6.2f \n',xHe)
%fprintf(1,'Phosphine mixing ratio is %6.2f \n',xPH3)






%Convert pressure percents to partial pressures
%% This version assumes already converted tp P_H2 etc...
dP_H2=P_H2;
dP_He=P_He;
dP_PH3=P_PH3;


%calculate vector linewidth
d_nu_GHz=lnwidth(To,T,'H2',dP_H2,'He',dP_He,'PH3',dP_PH3);
d_nu=d_nu_GHz*GHztoinv_cm;			%Correct for units of d_nu in GHz to (cm^-1)



%Parse data file
fo=catalog(:,1);			%Line freq (GHz)
Spp=catalog(:,2);			%intensity nm^2MHz
ELO=catalog(:,3);			%Energy of lower state 1/cm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%padsize=size(fo,1)-size(A,1);
%padding=ones(padsize,1);
%AA=[A', padding']';
%Spp(1:40)=Spp(1:40).*A;
%a1=[15 17 19 21 23 25 26 28 30 32 33 34 35 36 37 38 39 40 41];
%a2=[2:14 16 18 20 22 24 27 29 31];
%a1=a1-1;
%a2=a2-1;
%a1=[1 4 7 9 11 14 16 18 20 22 24 25 27 29 31:40];
%a2=[2 3 5 6 8 10 12 13 15 17 19 21 23 26 28 30]; 

%Spp(a1)=Spp(a1).*A(1);
%Spp(a2)=Spp(a2).*A(2);
%Spp(2:41)=Spp(2:41).*A(2);
%Spp(41:105)=Spp(41:105).*A(3);

Spp(a1)=Spp(a1).*(36.65);
Spp(a2)=Spp(a2).*(2.76);

% Convert Intensity at 300 Kelvin to Intensity at T kelvin:

numberdensity=(10^6)*No*dP_PH3/(R*T);	% [molec/cm^3]
												% The 10^6 allows use of P(bar) for P(dynes/cm^2)

eta=3/2;										% for symmetric top molecule
expo=-(1/T-1/To)*ELO*hc/K;
ST=Spp.*exp(expo).*(To/T)^(eta+1);	% S(T) =S(To)converted for temperature

alpha_max=ST.*((d_nu).^-1)*((pi*2.9979E18)^-1)*numberdensity; % 


%%%% Alpha Max Found


%%%% Get Lineshape
% Lineshapes calculated using inverse centimeters
f_cm=f.*GHztoinv_cm;					% Convert observation freq to cm^-1
fo_cm=fo.*GHztoinv_cm;

Fvvw=vvwlineshape(f_cm,fo_cm,d_nu);	% Matrix of form f(1)*fo(1)  f(2)*fo(1) f(3)*fo(1)
%												  f(1)*fo(2)  f(2)*fo(2) f(3)*fo(2)


%ce=couplingelement(A,fo,d_nu_GHz,dP_H2,dP_He,dP_PH3,T);		% answer in GHz
%pst=pshiftterm(fo,dP_H2,dP_He,dP_PH3);							% answer in GHz

%convert to inverse cm
%ce_cm=ce.*GHztoinv_cm;
%pst_cm=pst.*GHztoinv_cm;
%zed=zeros(size(fo,1),1);

%Fbr=brlineshape(f_cm,fo_cm,d_nu,ce_cm,pst_cm);
% Fbr=brlineshape(f_cm,fo_cm,d_nu,zed,zed); %Becomes VVW lineshape

% Put alpha_max(ELO) into a matrix m x n like ie. Fvvw


n=size(f,2);
nones=ones(1,n);
alpha_max_matrix=alpha_max*nones;

%Multiply element by element alpha(1)*F(f,lf(1)) etc...
vvw_alpha_matrix=Fvvw.*alpha_max_matrix;
%br_alpha_matrix=Fbr.*alpha_max_matrix;

%Sum each column, which corresponds to one observation frequency.
%ie if 1-10Ghz, first column are the contributions of all line freq at 1 Ghz
vvw_alpha_neper=sum(vvw_alpha_matrix,1).*(pi);
alphaph3=vvw_alpha_neper.*OpticaldepthstodB;

%br_alpha_neper=sum(br_alpha_matrix,1).*(pi);
%br_alpha=br_alpha_neper.*OpticaldepthstodB;
